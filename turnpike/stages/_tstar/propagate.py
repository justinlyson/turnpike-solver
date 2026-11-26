"""Exact 0/1 bound propagation (compiled with numba).

The rule (Lemma S of the theory): for a system {Aeq z = beq, Aub z <= bub, z in {0,1}^nv} with integer data, value t of a
variable z_k is allowed by row r iff  a t  lies in [b - maxrest, b - minrest]  (equality row) or  a t <= b - minrest
(<= row), where [minrest, maxrest] is the range of the rest of row r under the current bounds.  A disallowed value is
removed; a variable with no allowed value proves that the system has no 0/1 point.  Every 0/1 point survives every step,
so a contradiction is a PROOF.  The rules are monotone narrowings, so every evaluation order reaches the same greatest
fixpoint; the event-driven loops below re-check a row only when one of its variables changed.

Two kernels:

* `propagate_face`  - plain propagation on T*'s face system (the propagation inside T*'s node test).
* `propagate_alias` - the same rule with ALIASING (RULE+ item 5): when a point column x_p is fixed to 1, every pair
  column y(p, q) is merged into x_q (McCormick forces y = x_q); rows are evaluated with columns mapped to their
  representatives and coefficients merged.  Used by `closure.PropNode` (the propagation inside close_k).

Integer arithmetic only (int64; coefficients and right-hand sides are counts).
"""
from __future__ import annotations

import numpy as np
import scipy.sparse as sp
from numba import njit


# ------------------------------------------------------------------------------------------------ plain face propagation
@njit(cache=True)
def _cval(indptr, indices, data, r, k):
    """A[r, k] for a CSR matrix with sorted row indices (binary search)."""
    lo_ = indptr[r]
    hi_ = indptr[r + 1]
    while lo_ < hi_:
        mid = (lo_ + hi_) >> 1
        if indices[mid] < k:
            lo_ = mid + 1
        else:
            hi_ = mid
    return data[lo_]


@njit(cache=True)
def _prop_face(indptr, indices, data, cptr, crow, b, meq, nv):
    m = indptr.shape[0] - 1
    lo = np.zeros(nv, np.int64)
    hi = np.ones(nv, np.int64)
    minact = np.zeros(m, np.int64)
    maxact = np.zeros(m, np.int64)
    for r in range(m):
        s0 = 0
        s1 = 0
        for p in range(indptr[r], indptr[r + 1]):
            a = data[p]
            if a > 0:
                s1 += a
            else:
                s0 += a
        minact[r] = s0
        maxact[r] = s1
    # maxfree[r] >= max |a| over the still-free variables of row r.  A row whose slacks are all >= maxfree[r] cannot
    # remove any value, so it is skipped (stale-high values are safe).
    maxfree = np.zeros(m, np.int64)
    for r in range(m):
        mx = 0
        for p in range(indptr[r], indptr[r + 1]):
            v = data[p] if data[p] > 0 else -data[p]
            if v > mx:
                mx = v
        maxfree[r] = mx
    queue = np.empty(m, np.int64)
    inq = np.ones(m, np.bool_)
    for r in range(m):
        queue[r] = r
    head = 0
    size = m
    while size > 0:
        r = queue[head]
        head += 1
        if head == m:
            head = 0
        size -= 1
        inq[r] = False
        br = b[r]
        eq = r < meq
        sl = br - minact[r]
        if eq:
            su = maxact[r] - br
            if su < sl:
                sl = su
        if sl >= maxfree[r]:
            continue
        newmax = 0
        for p in range(indptr[r], indptr[r + 1]):
            k = indices[p]
            a = data[p]
            al = a * lo[k]
            ah = a * hi[k]
            cmin = al if al < ah else ah
            cmax = al if al > ah else ah
            minrest = minact[r] - cmin
            maxrest = maxact[r] - cmax
            if eq:
                ok0 = (br - maxrest <= 0) and (0 <= br - minrest)
                ok1 = (br - maxrest <= a) and (a <= br - minrest)
            else:
                ok0 = 0 <= br - minrest
                ok1 = a <= br - minrest
            can0 = lo[k] == 0 and ok0
            can1 = hi[k] == 1 and ok1
            if not can0 and not can1:
                return False
            nlo = lo[k]
            nhi = hi[k]
            if not can0:
                nlo = 1
            if not can1:
                nhi = 0
            if nlo != lo[k] or nhi != hi[k]:
                for q in range(cptr[k], cptr[k + 1]):
                    rr = crow[q]
                    aa = _cval(indptr, indices, data, rr, k)
                    ol = aa * lo[k]
                    oh = aa * hi[k]
                    nl = aa * nlo
                    nh = aa * nhi
                    minact[rr] += (nl if nl < nh else nh) - (ol if ol < oh else oh)
                    maxact[rr] += (nl if nl > nh else nh) - (ol if ol > oh else oh)
                    if not inq[rr]:
                        inq[rr] = True
                        tail = head + size
                        if tail >= m:
                            tail -= m
                        queue[tail] = rr
                        size += 1
                lo[k] = nlo
                hi[k] = nhi
            if lo[k] != hi[k]:
                v = a if a > 0 else -a
                if v > newmax:
                    newmax = v
        maxfree[r] = newmax
    return True


def propagate_face(Aeq, beq, Aub, bub) -> bool:
    """Propagation to the fixpoint on T*'s face system.  Returns False iff a contradiction is derived (= proof)."""
    nv = Aeq.shape[1]
    A = sp.vstack([Aeq, Aub]).tocsr() if Aub.shape[0] else Aeq.tocsr()
    A.sum_duplicates()
    A.sort_indices()                    # explicit zeros are kept: they are checked like any entry (a row-range test)
    b = (np.concatenate([np.asarray(beq, np.int64), np.asarray(bub, np.int64)]) if Aub.shape[0]
         else np.asarray(beq, np.int64))
    C = A.tocsc()
    C.sort_indices()
    return bool(_prop_face(A.indptr.astype(np.int64), A.indices.astype(np.int64), A.data.astype(np.int64),
                           C.indptr.astype(np.int64), C.indices.astype(np.int64), b, Aeq.shape[0], nv))


# ------------------------------------------------------------------------------------------------ propagation with aliasing
@njit(cache=True)
def _push(r, queue, inq, head, size, m):
    if not inq[r]:
        inq[r] = True
        t = head + size
        if t >= m:
            t -= m
        queue[t] = r
        size += 1
    return size


@njit(cache=True)
def _narrow(t, nlo, nhi, lo, hi, rep, cptr, crow, incp, inc, ispt, pa, pb, queue, inq, head, size, m):
    """Set the bounds of representative t; queue its rows and the rows of every column aliased to it; if t is a point
    column that became 1, alias each free pair column y(t, o) to o (when o is still free)."""
    lo[t] = nlo
    hi[t] = nhi
    for q in range(cptr[t], cptr[t + 1]):
        size = _push(crow[q], queue, inq, head, size, m)
    if ispt[t]:
        for s in range(incp[t], incp[t + 1]):
            y = inc[s]
            if rep[y] == t:
                for q in range(cptr[y], cptr[y + 1]):
                    size = _push(crow[q], queue, inq, head, size, m)
        if nlo == 1:
            for s in range(incp[t], incp[t + 1]):
                y = inc[s]
                if rep[y] != y or lo[y] == hi[y]:
                    continue
                o = pb[y] if pa[y] == t else pa[y]
                if lo[o] == hi[o]:
                    continue                    # other end fixed: the McCormick rows decide y
                rep[y] = o
                for q in range(cptr[y], cptr[y + 1]):
                    size = _push(crow[q], queue, inq, head, size, m)
                for q in range(cptr[o], cptr[o + 1]):
                    size = _push(crow[q], queue, inq, head, size, m)
    return size


@njit(cache=True)
def propagate_alias(indptr, indices, data, cptr, crow, b, meq, lo, hi, rep, incp, inc, ispt, pa, pb,
                    start, nstart, scratch, mark, lst):
    """Propagation with aliasing, in place on the state (lo, hi, rep).  Rows `start[:nstart]` are queued first
    (nstart < 0: all rows).  Returns False iff a contradiction is derived."""
    m = indptr.shape[0] - 1
    queue = np.empty(m, np.int64)
    inq = np.zeros(m, np.bool_)
    head = 0
    size = 0
    if nstart < 0:
        for r in range(m):
            size = _push(r, queue, inq, head, size, m)
    else:
        for i in range(nstart):
            size = _push(start[i], queue, inq, head, size, m)
    while size > 0:
        r = queue[head]
        head += 1
        if head == m:
            head = 0
        size -= 1
        inq[r] = False
        const = 0
        nl = 0
        for p in range(indptr[r], indptr[r + 1]):
            t = rep[indices[p]]
            a = data[p]
            if lo[t] == hi[t]:
                const += a * lo[t]
            else:
                if not mark[t]:
                    mark[t] = True
                    scratch[t] = 0
                    lst[nl] = t
                    nl += 1
                scratch[t] += a
        mn = const
        mx = const
        for i in range(nl):
            a = scratch[lst[i]]
            if a > 0:
                mx += a
            else:
                mn += a
        br = b[r]
        eq = r < meq
        if mn > br or (eq and mx < br):
            for i in range(nl):
                mark[lst[i]] = False
            return False
        for i in range(nl):
            t = lst[i]
            a = scratch[t]
            if lo[t] == hi[t] or a == 0:
                continue
            cmin = a if a < 0 else 0
            cmax = a if a > 0 else 0
            minrest = mn - cmin
            maxrest = mx - cmax
            if eq:
                ok0 = (br - maxrest <= 0) and (0 <= br - minrest)
                ok1 = (br - maxrest <= a) and (a <= br - minrest)
            else:
                ok0 = 0 <= br - minrest
                ok1 = a <= br - minrest
            if not ok0 and not ok1:
                for j in range(nl):
                    mark[lst[j]] = False
                return False
            if ok0 != ok1:
                v = 1 if ok1 else 0
                size = _narrow(t, v, v, lo, hi, rep, cptr, crow, incp, inc, ispt, pa, pb, queue, inq, head, size, m)
                break                           # row r is re-queued by _narrow (t is in r): activities changed
        for i in range(nl):
            mark[lst[i]] = False
    return True
