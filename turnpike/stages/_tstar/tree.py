"""The T* class tree: candidate classes, forced closure, and T*'s exact node test.

Setting.  D is the distance multiset, W = max D, and the end points 0 and W are placed first.  Every other point p of a
realization X lies in the candidate set {p : 0 < p < W, p in D, W - p in D}; candidates come in CLASSES c = min(p, W - p),
processed in increasing order.  At class c the number of still unused copies of the distance W - c,

    kappa(c) = D[W - c] - (copies already used by placed pairs),

says how many of {c, W - c} are in X (the pair (0, W - c) or (c, W) uses one copy each): kappa = 0 (neither),
kappa = 2 (both: forced), kappa = 1 (exactly one: a BRANCHING class, two children "c" and "W - c"), anything else: dead.
The forced classes up to the next branching class form the CLOSURE of a node (`_closure`); distances in the gap between
consecutive classes must be used up exactly when the gap closes (they can no longer be produced).  The first branching
class takes only the child "c" (mirror symmetry: X and W - X have the same distances).

T*'s node test (`test`), every part exact:
  (a) counting: enough candidates remain (raw and after the oracle);
  (b) the combinatorial oracle: a free candidate at a non-distance from a placed point, or over-using a distance with the
      placed points, is not in X; every distance still unused must be reachable by a placed-free or free-free pair;
  (c) degree peeling: a point of X among the free candidates has rem - 1 partners among them (the (rem-1)-core);
  (d) the SA_2 face system Z(v) (point columns x_p, pair columns y_pq; value rows, count row, degree rows; McCormick and
      forbidden-pair rows; box) - exact 0/1 propagation (`propagate.propagate_face`), then
  (e) its LP relaxation, refuted only through an exact integer Farkas certificate (`lp.certify`): first the equality rows
      (sub-systems of short rows first), then the full face.  HiGHS gets an ITERATION limit (no clock), so the same face
      always gets the same answer.  An LP that stops at its limit or whose ray fails the exact check only prunes less.

These are the THEORY settings of the research code (no caps, no LP back-off, LP on every child, propagation to the
fixpoint); the face system of the last tested node is kept (`last_face`) for the k-level test of `closure.py`.
"""
from __future__ import annotations

import bisect
from collections import Counter

import numpy as np
import scipy.sparse as sp

from ..  import TimeUp
from . import lp as LP
from .propagate import propagate_face

LP_ITER_FACTOR = 20          # HiGHS iteration limit = factor * (rows + cols) + floor
LP_ITER_FLOOR = 1000
CERT_BITS = 40               # scale of the rounded Farkas multipliers in `certify`
SUB_RAY_K = (4, 12)          # equality-stage sub-systems: rows with at most K nonzeros


class FaceData:
    """The face system Z(v) of one tested node, with what `closure.py` needs to label rows and columns."""

    __slots__ = ("Aeq", "beq", "Aub", "bub", "rows_used", "free", "EI", "EJ", "FI", "FJ", "decided")

    def __init__(self, **kw):
        for k, v in kw.items():
            setattr(self, k, v)


class ClassTree:
    """State of the T* tree on one list: placed points, used distance counts, and the node test.

    D: Counter distance -> multiplicity; n: number of points; deadline: object with .expired() or None.
    """

    def __init__(self, D, n, deadline=None, stats=None):
        self.D = D
        self.n = n
        self.deadline = deadline
        W = max(D)
        self.W = W
        self.pdt = np.int64 if W < 2 ** 62 else object
        self.vals = sorted(D)
        self.vid = {v: i for i, v in enumerate(self.vals)}
        self.V = len(self.vals)
        self.vals_np = np.array(self.vals, dtype=self.pdt)
        self.m_np = np.array([D[v] for v in self.vals], dtype=np.int64)
        self.cands = sorted({0, W} | {v for v in D if 0 < v < W and (W - v) in D})
        self.cands_np = np.array(self.cands, dtype=self.pdt)
        self.classes = sorted({min(p, W - p) for p in self.cands if 0 < p < W})
        # gap[ci] = index range of the distances in (W - classes[ci], W - classes[ci-1]): final once class ci is reached
        self.gap = []
        prev = 0
        for c in self.classes:
            i = bisect.bisect_right(self.vals, W - c)
            j = bisect.bisect_left(self.vals, W - prev)
            self.gap.append((i, j))
            prev = c
        self.used = np.zeros(self.V, dtype=np.int64)
        self.present = []
        self.tick = 0
        self.keep_face = False
        self.last_face = None
        self.st = stats if stats is not None else Counter()

    # ------------------------------------------------------------------ bookkeeping
    def _check_time(self):
        if self.deadline is not None and self.deadline.expired():
            raise TimeUp()

    def place(self, p):
        """Add point p if its distances to the placed points are still available; returns the undo record or None."""
        add = Counter()
        for q in self.present:
            i = self.vid.get(p - q if p > q else q - p)
            if i is None:
                return None
            add[i] += 1
        for i, k in add.items():
            if self.used[i] + k > self.m_np[i]:
                return None
        for i, k in add.items():
            self.used[i] += k
        self.present.append(p)
        return add

    def _unplace(self, p, add):
        q = self.present.pop()
        assert q == p
        for i, k in add.items():
            self.used[i] -= k

    def undo(self, trail):
        for p, a in reversed(trail):
            self._unplace(p, a)

    def closure(self, ci, trail):
        """Place every forced class from class index ci on.  Returns ('complete' | 'dead' | 'branch', class index)."""
        W, n, D, cl = self.W, self.n, self.D, self.classes
        while True:
            self.tick += 1
            if self.tick % 512 == 0:
                self._check_time()
            if ci == len(cl):
                return "complete", ci
            gi, gj = self.gap[ci]
            if gj > gi and not np.array_equal(self.used[gi:gj], self.m_np[gi:gj]):
                return "dead", ci
            c = cl[ci]
            if 2 * c == W:                                   # the middle class: one candidate point
                kap = n - len(self.present)
                if kap == 0:
                    ci += 1
                    continue
                if kap == 1:
                    a = self.place(c)
                    if a is None:
                        return "dead", ci
                    trail.append((c, a))
                    ci += 1
                    continue
                return "dead", ci
            j = self.vid.get(W - c)
            kap = D[W - c] - (int(self.used[j]) if j is not None else 0)
            if kap == 0:
                ci += 1
                continue
            if kap == 2:
                for p in (c, W - c):
                    a = self.place(p)
                    if a is None:
                        return "dead", ci
                    trail.append((p, a))
                ci += 1
                continue
            if kap == 1:
                return "branch", ci
            return "dead", ci

    def apply(self, ci, side):
        """Child `side` (0: point c, 1: point W - c) of the branching class ci, with its closure."""
        c = self.classes[ci]
        p = c if side == 0 else self.W - c
        trail = []
        a = self.place(p)
        if a is None:
            return "dead", ci, trail
        trail.append((p, a))
        s, ci2 = self.closure(ci + 1, trail)
        return s, ci2, trail

    def leaf(self):
        """The placed set if it uses D exactly (a realization), else None."""
        X = sorted(int(p) for p in self.present)
        if len(X) != self.n or np.any(self.used != self.m_np):
            return None
        return X

    # ------------------------------------------------------------------ the node test
    def _free(self, ci):
        """The oracle's free candidates at class ci and their (distance, candidate, count) incidences with placed points."""
        c = self.classes[ci]
        W = self.W
        lo = bisect.bisect_left(self.cands, c)
        hi = bisect.bisect_right(self.cands, W - c)
        free = self.cands_np[lo:hi]
        Z = np.array(self.present, dtype=self.pdt)
        F = len(free)
        if F == 0:
            return free, (np.zeros(0, np.int64),) * 3
        dz = np.abs(free[None, :] - Z[:, None])
        vi = np.searchsorted(self.vals_np, dz).astype(np.int64)
        vic = np.minimum(vi, self.V - 1)
        ok = (self.vals_np[vic] == dz).astype(bool)
        kill = ~ok.all(axis=0)                           # a distance to a placed point is not in D
        zz, pp = np.nonzero(ok)
        key = vic[zz, pp] * F + pp
        uk, cnt = np.unique(key, return_counts=True)
        kv = uk // F
        kp = uk % F
        over = self.used[kv] + cnt > self.m_np[kv]        # would over-use a distance together with the placed points
        kill[kp[over]] = True
        keep = ~kill
        newidx = -np.ones(F, dtype=np.int64)
        newidx[keep] = np.arange(int(keep.sum()))
        sel = keep[kp]
        return free[keep], (kv[sel], newidx[kp[sel]], cnt[sel])

    def test(self, ci) -> bool:
        """T*'s node test at the node whose next branching class is ci.  False = refuted (exactly)."""
        st = self.st
        self._check_time()
        self.last_face = None
        c = self.classes[ci]
        Z = len(self.present)
        rem = self.n - Z
        Fraw = bisect.bisect_right(self.cands, self.W - c) - bisect.bisect_left(self.cands, c)
        if Fraw < rem or rem < 0:
            return False
        free, (cv, cp, cc) = self._free(ci)
        F = len(free)
        if F < rem:
            return False
        resid = self.m_np - self.used
        reach = np.zeros(self.V, dtype=bool)
        reach[cv] = True
        if F >= 2:
            iu, ju = np.triu_indices(F, 1)
            d = free[ju] - free[iu]
            vi = np.searchsorted(self.vals_np, d).astype(np.int64)
            vic = np.minimum(vi, self.V - 1)
            okp = (self.vals_np[vic] == d).astype(bool)
            reach[vic[okp]] = True
        else:
            iu = ju = vic = np.zeros(0, dtype=np.int64)
            okp = np.zeros(0, dtype=bool)
        if np.any((resid != 0) & ~reach):
            return False
        if F == 0:
            return rem == 0 and not np.any(resid)
        if F >= 2 and rem >= 2:
            # degree peeling: the (rem - 1)-core of the graph of free candidates joined by distance pairs
            need = rem - 1
            alive = np.ones(F, dtype=bool)
            EI0, EJ0 = iu[okp], ju[okp]
            while True:
                m = alive[EI0] & alive[EJ0]
                deg = np.bincount(EI0[m], minlength=F) + np.bincount(EJ0[m], minlength=F)
                bad = alive & (deg < need)
                if not bad.any():
                    break
                alive[bad] = False
            na = int(alive.sum())
            if na < F:
                if na < rem:
                    return False
                newidx = -np.ones(F, dtype=np.int64)
                newidx[alive] = np.arange(na)
                selc = alive[cp]
                cv, cp, cc = cv[selc], newidx[cp[selc]], cc[selc]
                selp = alive[iu] & alive[ju]
                iu, ju, vic, okp = newidx[iu[selp]], newidx[ju[selp]], vic[selp], okp[selp]
                free = free[alive]
                F = na
                reach = np.zeros(self.V, dtype=bool)
                reach[cv] = True
                reach[vic[okp]] = True
                if np.any((resid != 0) & ~reach):
                    return False
        return not self._face_refuted(free, cv, cp, cc, iu, ju, vic, okp, resid, Z, rem)

    # ------------------------------------------------------------------ the SA_2 face: propagation, then exact LP
    def _face_refuted(self, free, cv, cp, cc, iu, ju, vic, okp, resid, Z, rem) -> bool:
        """Build the face system Z(v); True only if exact propagation or an exactly certified LP refutes it."""
        st = self.st
        self._check_time()
        F = len(free)
        n = self.n
        EI, EJ, EV = iu[okp], ju[okp], vic[okp]
        FI, FJ = iu[~okp], ju[~okp]
        E = len(EI)
        nv = F + E
        ar = np.arange(E)
        af = np.arange(F)
        # equality rows: one per distance value in use (count of remaining copies), the count row, one degree row per point
        rows_used, inv = np.unique(np.concatenate([cv, EV]), return_inverse=True)
        R = len(rows_used)
        r_ = [inv[:len(cv)], inv[len(cv):], np.full(F, R), R + 1 + af, R + 1 + EI, R + 1 + EJ]
        c_ = [cp, F + ar, af, af, F + ar, F + ar]
        a_ = [cc.astype(np.int64), np.ones(E, np.int64), np.ones(F, np.int64), np.full(F, Z - (n - 1), dtype=np.int64),
              np.ones(E, np.int64), np.ones(E, np.int64)]
        meq = R + 1 + F
        Aeq = sp.csr_matrix((np.concatenate(a_), (np.concatenate(r_), np.concatenate(c_))), shape=(meq, nv), dtype=np.int64)
        beq = np.concatenate([resid[rows_used].astype(np.int64), [rem], np.zeros(F, np.int64)])
        # <= rows: McCormick y <= x_i, y <= x_j, x_i + x_j - y <= 1, and x_i + x_j <= 1 for pairs at a non-distance
        nf = len(FI)
        r_ = [ar, ar, E + ar, E + ar, 2 * E + ar, 2 * E + ar, 2 * E + ar, 3 * E + np.arange(nf), 3 * E + np.arange(nf)]
        c_ = [F + ar, EI, F + ar, EJ, EI, EJ, F + ar, FI, FJ]
        one, mone = np.ones(E, np.int64), -np.ones(E, np.int64)
        a_ = [one, mone, one, mone, one, one, mone, np.ones(nf, np.int64), np.ones(nf, np.int64)]
        mub = 3 * E + nf
        Aub = sp.csr_matrix((np.concatenate(a_), (np.concatenate(r_), np.concatenate(c_))), shape=(mub, nv), dtype=np.int64)
        bub = np.concatenate([np.zeros(2 * E, np.int64), np.ones(E, np.int64), np.ones(nf, np.int64)])
        if self.keep_face:
            self.last_face = FaceData(Aeq=Aeq, beq=beq, Aub=Aub, bub=bub, rows_used=rows_used, free=free, EI=EI, EJ=EJ,
                                      FI=FI, FJ=FJ, decided=sorted(int(p) for p in self.present))
        st["propagations"] += 1
        if not propagate_face(Aeq, beq, Aub, bub):
            st["prop_refuted"] += 1
            return True
        return self._lp_refuted(Aeq, beq, Aub, bub)

    def _iters(self, A):
        return LP.iteration_limit(A, LP_ITER_FACTOR, LP_ITER_FLOOR)

    def _lp_refuted(self, Aeq, beq, Aub, bub) -> bool:
        """Two LP stages (equality rows, then the full face), each by dual simplex with an iteration limit; a stage cuts
        only through `certify` (exact)."""
        st = self.st
        st["lp_calls"] += 1
        E0 = sp.csr_matrix((0, Aeq.shape[1]), dtype=np.int64)
        b0 = np.zeros(0, np.int64)
        z0 = np.zeros(0)
        s0, ray0 = LP.hsolve(Aeq, beq, E0, b0, self._iters(Aeq))
        if s0 == "unknown":
            st["lp_unknown"] += 1
        if s0 == "infeasible":
            for K in SUB_RAY_K:
                yk = LP.sub_ray(Aeq, beq, K, self._iters)
                if yk is not None and LP.certify_either_sign(Aeq, beq, E0, b0, yk, z0, CERT_BITS):
                    st["lp_refuted"] += 1
                    return True
            if ray0 is not None and LP.certify_either_sign(Aeq, beq, E0, b0, ray0, z0, CERT_BITS):
                st["lp_refuted"] += 1
                return True
            st["lp_cert_fail"] += 1
        meq = Aeq.shape[0]
        s1, ray = LP.hsolve(Aeq, beq, Aub, bub, self._iters(sp.vstack([Aeq, Aub])))
        if s1 == "unknown":
            st["lp_unknown"] += 1
        if s1 != "infeasible":
            if s0 == "infeasible":
                st["lp_cert_fail"] += 1
            return False
        ok = ray is not None and LP.certify_either_sign(Aeq, beq, Aub, bub, ray[:meq], ray[meq:], CERT_BITS)
        if ok:
            st["lp_refuted"] += 1
        else:
            st["lp_cert_fail"] += 1
        return ok
