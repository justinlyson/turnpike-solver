"""The k-level child test close_k^LP (theory: PROOF_STPLUS sections 8.5, 8.12, 8.16; Lemma S2^LP).

A child that passed T*'s node test has a face system (tree.FaceData).  To it we add the Q rows (integer combinations of
equality rows, below), and work with STATES s = sets of 0/1 fixes of the columns.

Base test  B(s):  the propagation with aliasing from s reaches a contradiction, OR the exact face LP with the fixes of s
                  is infeasible (decided only through an exactly checked Farkas certificate, lp.FaceLP).
close_0(s) refutes iff B(s).  For d >= 1, close_d(s): repeat over the free point columns x -
    if B(current): REFUTED;  a leg x = v is refuted iff (d = 1) B(current + {x = v} propagated), or (d > 1)
    close_{d-1} refutes current + {x = v} (propagated);  both legs refuted: REFUTED;  one leg refuted: fix x to the
    other value and continue;  ALIVE at a fixpoint.
T*_k  =  T*'s node test, then KILL iff close_k refutes the child's (propagated) root state.

Soundness (Lemma S2^LP): propagation and the certified LP are sound; by induction on d a refuted leg has no 0/1
completion, so a refuted state has none, and T*_k never refuses a child that extends to a realization.
Monotonicity: B is monotone in s, so close_d is a closure under monotone inflationary operators: the verdict does not
depend on the column order.  Theorem ST (8.16): close_k refutes s iff s has a split refutation tree of Strahler rank
<= k.  Cost: at most (2c(c+1)+1)^k propagations and exact LP decisions per child test (c = number of point columns),
polynomial for fixed k.

Q rows (RULE+ item 2).  Let A_p be the decided points below W/2 and B_p = {W - p : p decided above W/2}.  For a decided
side P, c = 1 / P(z) as a power series (P(z) = sum_{p in P} z^p, constant term 1); Q(k) = sum_j c_j L(k - j) over the
face's top-window count rows L(k) (distance W - k).  Any integer combination of equality rows is valid, so the Q rows
are sound; the series is truncated before its first coefficient above 2^40 and rows whose absolute row sum could exceed
2^59 are dropped (both only weaken the test, and keep the int64 arithmetic exact).
"""
from __future__ import annotations

import numpy as np
import scipy.sparse as sp

from .. import TimeUp
from .lp import FaceLP
from .propagate import propagate_alias

Q_CAP = 2 ** 40              # series coefficients above this are dropped (truncation = weaker, sound)
Q_JMAX = 2 * 10 ** 6         # longest series computed


# ------------------------------------------------------------------------------------------------ labels and Q rows
def face_labels(T, face):
    """Row labels of the equality rows ('L', W - v) / ('v', v) / ('count',) / ('deg', p) and column labels
    (point p, or ('y', p, q) for a pair column), as the face system was built in tree.ClassTree._face_refuted."""
    W = T.W
    fl = [int(p) for p in face.free]

    def vlab(v):
        v = int(v)
        return ("L", W - v) if 2 * v > W else ("v", v)

    eql = [vlab(T.vals[i]) for i in face.rows_used] + [("count",)] + [("deg", p) for p in fl]
    coll = fl + [("y", fl[i], fl[j]) for i, j in zip(face.EI, face.EJ)]
    return eql, coll


def series_inverse(P, m):
    """c = 1 / P(z) mod z^(J+1) for a set P of positive exponents (P(z) = 1 + sum z^p), int64, J <= min(m, Q_JMAX),
    stopped before the first |c_t| > Q_CAP.  Computed in blocks of length min(P): a block depends only on earlier ones."""
    P = np.array(sorted(set(int(i) for i in P if i >= 1)), np.int64)
    if len(P) == 0:
        return np.ones(1, np.int64)
    m = min(int(m), Q_JMAX)
    p1 = int(P[0])
    c = np.zeros(m + 1, np.int64)
    c[0] = 1
    t = 1
    while t <= m:
        e = min(t + p1, m + 1)
        blk = np.zeros(e - t, np.int64)
        for i in P:
            if i > e - 1:
                break
            lo_ = t - i
            hi_ = e - i
            if hi_ <= 0:
                continue
            if lo_ < 0:
                blk[-lo_:] -= c[0:hi_]
            else:
                blk -= c[lo_:hi_]
        big = np.nonzero(np.abs(blk) > Q_CAP)[0]
        if len(big):
            J = t + int(big[0]) - 1
            c[t:J + 1] = blk[:J + 1 - t]
            return c[:J + 1]
        c[t:e] = blk
        t = e
    return c


def q_rows(Aeq, beq, eql, Pset):
    """Q(k) = sum_{k2 <= k} c_{k - k2} L(k2) over the face's L rows, as one exact sparse product.  Returns (rows as
    dicts column -> coefficient, right-hand sides)."""
    Lidx = {l[1]: r for r, l in enumerate(eql) if l[0] == "L"}
    if not Lidx:
        return [], []
    ks = np.array(sorted(Lidx), np.int64)
    rows_of = np.array([Lidx[int(k)] for k in ks], np.int64)
    c = series_inverse(Pset, int(ks[-1] - ks[0]))
    J = len(c) - 1
    ri, ci, vv = [], [], []
    for a, k in enumerate(ks):
        k2 = ks[:a + 1]
        d = k - k2
        msk = d <= J
        cf = c[d[msk]]
        nzm = cf != 0
        if nzm.sum() <= 1:
            continue                                   # a single term is just an L row
        ri.append(np.full(int(nzm.sum()), a))
        ci.append(rows_of[:a + 1][msk][nzm])
        vv.append(cf[nzm])
    if not ri:
        return [], []
    ri = np.concatenate(ri)
    ci = np.concatenate(ci)
    vv = np.concatenate(vv)
    C = sp.csr_matrix((vv, (ri, ci)), shape=(len(ks), Aeq.shape[0]), dtype=np.int64)
    Af = sp.csr_matrix(Aeq, dtype=np.int64)
    bound = ((abs(C).astype(np.float64) @ abs(Af).astype(np.float64)).sum(axis=1).A.ravel()
             + np.abs(C.astype(np.float64) @ np.asarray(beq, np.float64)))
    keep = np.nonzero((np.diff(C.indptr) > 1) & (bound < 2.0 ** 59))[0]
    C = C[keep]
    Q = (C @ Af).tocsr()
    Q.eliminate_zeros()
    qb = C @ np.asarray(beq, np.int64)
    R, b = [], []
    for i in range(Q.shape[0]):
        if Q.indptr[i + 1] == Q.indptr[i]:
            continue
        sl = slice(Q.indptr[i], Q.indptr[i + 1])
        R.append({int(k2): int(v) for k2, v in zip(Q.indices[sl], Q.data[sl])})
        b.append(int(qb[i]))
    return R, b


def rule_system(T, face):
    """The child's system for close_k: face equality rows + Q rows (Ae, be), the <= rows (Au, bu), column labels and the
    number of face equality rows (the LP uses the face rows only)."""
    W = T.W
    Aeq, beq, Aub, bub = face.Aeq, face.beq, face.Aub, face.bub
    eql, coll = face_labels(T, face)
    dec = face.decided
    Ap = [p for p in dec if 2 * p < W]
    Bp = [W - p for p in dec if 2 * p > W]
    QR, Qb = [], []
    for P in (Ap, Bp):
        r, b = q_rows(Aeq, beq, eql, P)
        QR += r
        Qb += b
    ncol = Aeq.shape[1]
    if QR:
        rr, cc, vv = [], [], []
        for i, row in enumerate(QR):
            for k2, v in row.items():
                rr.append(i)
                cc.append(k2)
                vv.append(v)
        Q = sp.csr_matrix((np.array(vv, np.int64), (rr, cc)), shape=(len(QR), ncol), dtype=np.int64)
        Ae = sp.vstack([Aeq.astype(np.int64), Q]).tocsr()
        be = np.concatenate([beq, np.array(Qb, np.int64)])
    else:
        Ae = Aeq.astype(np.int64).tocsr()
        be = beq
    return Ae, be, Aub.astype(np.int64).tocsr(), bub, coll, Aeq.shape[0]


# ------------------------------------------------------------------------------------------------ propagation states
class PropNode:
    """Propagation with aliasing on one child's system.  A state is (lo, hi, rep): bounds per column and the alias
    representative of each column.  `fix` copies the state, fixes one column and propagates from the touched rows."""

    def __init__(self, Ae, be, Au, bu, coll):
        A = sp.vstack([sp.csr_matrix(Ae, dtype=np.int64), sp.csr_matrix(Au, dtype=np.int64)]).tocsr()
        A.sum_duplicates()
        A.sort_indices()
        A.eliminate_zeros()
        self.indptr = A.indptr.astype(np.int64)
        self.indices = A.indices.astype(np.int64)
        self.data = A.data.astype(np.int64)
        C = A.tocsc()
        C.sort_indices()
        self.cptr = C.indptr.astype(np.int64)
        self.crow = C.indices.astype(np.int64)
        self.b = np.concatenate([np.asarray(be, np.int64), np.asarray(bu, np.int64)])
        self.meq = Ae.shape[0]
        nv = A.shape[1]
        self.nv = nv
        self.coll = coll
        pos = {l: i for i, l in enumerate(coll) if not isinstance(l, tuple)}
        self.ispt = np.array([not isinstance(l, tuple) for l in coll], np.bool_)
        pa = np.arange(nv, dtype=np.int64)
        pb = np.arange(nv, dtype=np.int64)
        incl = [[] for _ in range(nv)]
        for j, l in enumerate(coll):
            if isinstance(l, tuple):
                pa[j] = pos[l[1]]
                pb[j] = pos[l[2]]
                incl[pa[j]].append(j)
                incl[pb[j]].append(j)
        self.pa, self.pb = pa, pb
        self.incp = np.zeros(nv + 1, np.int64)
        self.incp[1:] = np.cumsum([len(x) for x in incl])
        self.inc = np.array([y for x in incl for y in x], np.int64) if self.incp[-1] else np.zeros(0, np.int64)
        self.scratch = np.zeros(nv, np.int64)
        self.mark = np.zeros(nv, np.bool_)
        self.lst = np.zeros(nv, np.int64)
        self.pts = [i for i in range(nv) if self.ispt[i]]
        self.calls = 0

    def _run(self, st, start, nstart):
        self.calls += 1
        lo, hi, rep = st
        return propagate_alias(self.indptr, self.indices, self.data, self.cptr, self.crow, self.b, self.meq, lo, hi, rep,
                               self.incp, self.inc, self.ispt, self.pa, self.pb, start, nstart,
                               self.scratch, self.mark, self.lst)

    def root(self):
        """Propagate from the empty state (all rows queued).  Returns (ok, state)."""
        st = (np.zeros(self.nv, np.int64), np.ones(self.nv, np.int64), np.arange(self.nv, dtype=np.int64))
        ok = self._run(st, np.zeros(1, np.int64), -1)
        return ok, st

    def fix(self, st, k, v):
        """Copy of st with column k fixed to v, propagated.  Returns (ok, state); ok False = contradiction (proof)."""
        lo, hi, rep = (st[0].copy(), st[1].copy(), st[2].copy())
        t = rep[k]
        if lo[t] == hi[t]:
            return (lo[t] == v), (lo, hi, rep)
        lo[t] = v
        hi[t] = v
        start = list(self.crow[self.cptr[t]:self.cptr[t + 1]])
        if self.ispt[t]:
            for y in self.inc[self.incp[t]:self.incp[t + 1]]:
                if rep[y] == t:
                    start += list(self.crow[self.cptr[y]:self.cptr[y + 1]])
            if v == 1:
                for y in self.inc[self.incp[t]:self.incp[t + 1]]:
                    if rep[y] != y or lo[y] == hi[y]:
                        continue
                    o = self.pb[y] if self.pa[y] == t else self.pa[y]
                    if lo[o] == hi[o]:
                        continue
                    rep[y] = o
                    start += list(self.crow[self.cptr[y]:self.cptr[y + 1]]) + list(self.crow[self.cptr[o]:self.cptr[o + 1]])
        s = np.asarray(sorted(set(start)) if start else [0], np.int64)
        ok = self._run((lo, hi, rep), s, len(s))
        return ok, (lo, hi, rep)


# ------------------------------------------------------------------------------------------------ close_k^LP
class CloseK:
    """close_k^LP on one child (see the module docstring).  refutes(k) -> True iff close_k refutes the child."""

    def __init__(self, T, face, deadline=None):
        W = T.W
        Ae, be, Au, bu, coll, nq_face = rule_system(T, face)
        self.node = PropNode(Ae, be, Au, bu, coll)
        self.lp = FaceLP(Ae, be, Au, bu, nq_face)
        self.order = sorted(self.node.pts, key=lambda i: min(coll[i], W - coll[i]))   # near-end columns first
        self.deadline = deadline

    def _base(self, st) -> bool:
        """B(s) for an already propagated state: the exact face LP (None = unknown = not refuted)."""
        return self.lp.dead(st) is True

    def _close(self, st, d):
        """Returns (alive, state at the fixpoint)."""
        if self._base(st):
            return False, st
        N = self.node
        changed = True
        while changed:
            changed = False
            for k in self.order:
                if self.deadline is not None and self.deadline.expired():
                    raise TimeUp()
                lo, hi, rep = st
                t = rep[k]
                if lo[t] == hi[t]:
                    continue
                alive = []
                for v in (0, 1):
                    ok, s2 = N.fix(st, k, v)
                    if not ok:
                        continue
                    if d == 1:
                        if self._base(s2):
                            continue
                        alive.append(s2)
                    elif self._close(s2, d - 1)[0]:
                        alive.append(s2)
                if not alive:
                    return False, st
                if len(alive) == 1:
                    st = alive[0]
                    changed = True
        return True, st

    def refutes(self, k: int) -> bool:
        ok, st = self.node.root()
        if not ok:
            return True
        if k <= 0:
            return self._base(st)
        return not self._close(st, k)[0]
