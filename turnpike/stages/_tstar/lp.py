"""Linear programs of the T* search and their exact certificates.

Every LP here is a feasibility question about a polytope

    P = { z : Aeq z = beq,  Aub z <= bub,  0 <= z <= 1 }

with integer data (the face system of a node, see `tree.py`).  The LPs are solved in floating point by HiGHS, and a
floating-point answer is NEVER used to refute anything by itself.  A node is refuted by an LP only through `certify`:
an integer Farkas vector Y = (Ye, Yu), Yu <= 0, with

    Y.b  >  max_{z in [0,1]^nv} (A^T Y).z ,

checked in exact integer arithmetic.  For every z in P we would have Y.b = Y.A z + Yu.(bub - Aub z) <= (A^T Y).z, so such
a Y proves P empty, hence (P contains every 0/1 completion of the node) that the node has no realization.

An LP that stops at its iteration / time limit, or whose dual ray does not pass the exact check, therefore only prunes
LESS: the search stays sound, it may just do more work.
"""
from __future__ import annotations

import math

import numpy as np
import scipy.sparse as sp


# ------------------------------------------------------------------------------------------------ exact certificate check

def _new_highs(mod=None):
    """A HiGHS instance with its one-time start-up banner suppressed (it is printed before options can be set)."""
    import os
    if mod is None:
        import highspy as mod
    try:
        fd = os.dup(1)
        with open(os.devnull, "w") as dn:
            os.dup2(dn.fileno(), 1)
            try:
                h = mod.Highs()
            finally:
                os.dup2(fd, 1)
                os.close(fd)
    except OSError:                       # no real stdout file descriptor (e.g. some notebooks)
        h = mod.Highs()
    h.setOptionValue("output_flag", False)
    return h

def _absmax(M) -> int:
    """Largest absolute entry of M (sparse or dense) as a Python integer (no int64 wrap-around)."""
    d = M.tocsr().data if sp.issparse(M) else np.asarray(M).ravel()
    return 0 if d.size == 0 else max(abs(int(d.max())), abs(int(d.min())))


def certify(Aeq, beq, Aub, bub, ye, yu, bits: int = 40) -> bool:
    """Exact check that the (floating) multipliers (ye, yu) round to an integer Farkas certificate for P.

    The vector is scaled so that its largest entry is about 2^bits, rounded to integers (yu clipped to <= 0), and the
    inequality Y.b > sum of the positive entries of A^T Y is decided in integers (int64 with an a-priori overflow bound,
    the final sums in Python integers).  Guards: non-finite input, tiny vectors and huge matrix entries return False.
    """
    if not (np.all(np.isfinite(ye)) and np.all(np.isfinite(yu))):
        return False
    s = max(float(np.max(np.abs(ye))) if ye.size else 0.0, float(np.max(np.abs(yu))) if yu.size else 0.0)
    if not (s > 0) or not math.isfinite(s):
        return False
    if max(_absmax(Aeq), _absmax(Aub)) * max(1, Aeq.shape[0] + Aub.shape[0]) >= 2 ** 61:
        return False
    colsum = np.asarray(abs(Aeq).sum(axis=0)).ravel() + (np.asarray(abs(Aub).sum(axis=0)).ravel() if Aub.shape[0] else 0)
    cmax = int(colsum.max()) if colsum.size else 1
    b = min(bits, 61 - max(1, cmax).bit_length())
    if b < 8:
        return False
    sc = float(2 ** b) / s
    if not math.isfinite(sc):
        return False
    Yef = np.rint(ye * sc)
    Yuf = np.minimum(np.rint(yu * sc), 0) if yu.size else np.zeros(0)
    lim = float(2 ** b)
    if not (np.all(np.isfinite(Yef)) and np.all(np.abs(Yef) <= lim) and np.all(np.isfinite(Yuf)) and np.all(np.abs(Yuf) <= lim)):
        return False
    Ye = Yef.astype(np.int64)
    Yu = Yuf.astype(np.int64) if yu.size else np.zeros(0, np.int64)
    rvec = Aeq.T.tocsr() @ Ye                     # |entries| <= 2^b * cmax < 2^61: exact in int64
    if Yu.size:
        rvec = rvec + Aub.T.tocsr() @ Yu
    Yb = int(np.dot(Ye.astype(object), np.asarray(beq).astype(object)))
    if Yu.size:
        Yb += int(np.dot(Yu.astype(object), np.asarray(bub).astype(object)))
    mx = int(np.sum(rvec[rvec > 0].astype(object))) if np.any(rvec > 0) else 0
    return Yb - mx > 0


def certify_either_sign(Aeq, beq, Aub, bub, ye, yu, bits: int) -> bool:
    """`certify` for (ye, yu) or for (-ye, -yu) (a dual ray is only defined up to sign)."""
    return certify(Aeq, beq, Aub, bub, ye, yu, bits) or certify(Aeq, beq, Aub, bub, -ye, -yu, bits)


# ------------------------------------------------------------------------------------------------ HiGHS wrappers
def _highs_model(A, row_lo, row_hi, col_lo, col_hi, cost):
    import highspy
    lp = highspy.HighsLp()
    lp.num_col_ = A.shape[1]
    lp.num_row_ = A.shape[0]
    lp.col_cost_ = cost
    lp.col_lower_ = col_lo
    lp.col_upper_ = col_hi
    lp.row_lower_ = row_lo
    lp.row_upper_ = row_hi
    lp.a_matrix_.format_ = highspy.MatrixFormat.kColwise
    lp.a_matrix_.start_ = A.indptr
    lp.a_matrix_.index_ = A.indices
    lp.a_matrix_.value_ = A.data
    return lp


def _stack(Aeq, beq, Aub, bub):
    mub = Aub.shape[0]
    A = sp.vstack([Aeq, Aub]).tocsc().astype(float) if mub else Aeq.tocsc().astype(float)
    lo = np.concatenate([np.asarray(beq, float), -np.inf * np.ones(mub)])
    hi = np.concatenate([np.asarray(beq, float), np.asarray(bub, float)])
    return A, lo, hi


def hsolve(Aeq, beq, Aub, bub, iters: int):
    """Feasibility of P by HiGHS dual simplex with an ITERATION limit (one thread; deterministic).

    Returns ('infeasible', dual ray or None) / ('feasible', None) / ('unknown', None).  Floating point: never cuts by itself.
    """
    import highspy
    nv = Aeq.shape[1]
    A, lo, hi = _stack(Aeq, beq, Aub, bub)
    h = _new_highs()
    h.setOptionValue("output_flag", False)
    h.setOptionValue("solver", "simplex")
    h.setOptionValue("threads", 1)
    h.setOptionValue("simplex_iteration_limit", int(iters))
    h.passModel(_highs_model(A, lo, hi, np.zeros(nv), np.ones(nv), np.zeros(nv)))
    h.run()
    s = h.getModelStatus()
    MS = highspy.HighsModelStatus
    if s == MS.kInfeasible:
        try:
            _, has, vals = h.getDualRay()
            ray = np.asarray(vals, float) if has else None
        except Exception:
            ray = None
        return "infeasible", ray
    if s == MS.kOptimal:
        return "feasible", None
    return "unknown", None


def iteration_limit(A, factor: int, floor: int) -> int:
    """The HiGHS iteration limit of a face LP: factor * (rows + cols) + floor (no clock in any decision)."""
    return factor * (A.shape[0] + A.shape[1]) + floor


def sub_ray(Aeq, beq, K: int, iters_of):
    """Dual ray of the sub-system formed by the rows of Aeq with at most K nonzeros (and the columns they meet).

    A sub-system of the equality stage: if it is infeasible, so is the face.  The ray is padded with zeros to all rows of
    Aeq; the caller uses it only through `certify` on the WHOLE equality system.  Returns the padded ray or None.
    """
    A = Aeq.tocsr()
    ri = np.nonzero(np.diff(A.indptr) <= K)[0]
    if len(ri) == 0 or len(ri) == A.shape[0]:
        return None
    As = A[ri]
    cols = np.unique(As.indices)
    if len(cols) == 0:
        return None
    As = As[:, cols].tocsr()
    s, ray = hsolve(As, beq[ri], sp.csr_matrix((0, len(cols)), dtype=np.int64), np.zeros(0, np.int64), iters_of(As))
    if s != "infeasible" or ray is None or len(ray) != len(ri):
        return None
    full = np.zeros(A.shape[0])
    full[ri] = ray
    return full


# ------------------------------------------------------------------------------------------------ the LP oracle of close_k
class FaceLP:
    """Exact LP refutation of a STATE of one child's face system (the level-0 test of close_k^LP, `closure.py`).

    A state (lo, hi, rep) fixes some columns (aliased columns follow their representative); the LP is the face system
    (face equality rows + all <= rows, box [0, 1]) with those bounds.  One persistent HiGHS model is re-solved warm (dual
    simplex) with the state applied as column bounds; a small pool of feasible LP points answers 'feasible' at once for
    any state whose bounds they satisfy.  'infeasible' is returned only after an exactly checked Farkas certificate of
    the substituted system (fixed columns moved to the right-hand side).

    dead(state) -> True (exactly refuted) / False (an LP point exists) / None (unknown: treated as NOT refuted).
    """

    WARM_TIME = 5.0
    IPM_TIME = 120.0
    FARKAS_TIME = 300.0
    POOL = 64

    def __init__(self, Aeq, beq, Aub, bub, nq_face: int):
        import highspy
        self.hs = highspy
        self.Af = sp.csc_matrix(Aeq[:nq_face, :], dtype=np.int64)
        self.bf = np.asarray(beq[:nq_face], np.int64)
        self.Au = sp.csc_matrix(Aub, dtype=np.int64)
        self.bu = np.asarray(bub, np.int64)
        A, lo, hi = _stack(self.Af, self.bf, self.Au, self.bu)
        nv = A.shape[1]
        self.nv = nv
        h = _new_highs()
        h.setOptionValue("output_flag", False)
        h.setOptionValue("solver", "simplex")
        h.setOptionValue("time_limit", self.IPM_TIME)
        h.passModel(_highs_model(A, lo, hi, np.zeros(nv), np.ones(nv), np.zeros(nv)))
        self.h = h
        self.curL = np.zeros(nv)
        self.curH = np.ones(nv)
        self.pool = []
        self.calls = 0
        self.unknown = 0

    def _fits(self, L, H):
        for x in self.pool:
            if np.all(x >= L - 1e-9) and np.all(x <= H + 1e-9):
                return True
        return False

    def _remember(self, x):
        self.pool.append(x)
        if len(self.pool) > self.POOL:
            self.pool.pop(0)

    def _substituted(self, st):
        lo, hi, rep = st
        Li = lo[rep]
        Hi = hi[rep]
        fixed = Li == Hi
        free = ~fixed
        f1 = np.nonzero(fixed & (Li == 1))[0]
        b1 = self.bf - np.asarray(self.Af[:, f1].sum(axis=1)).ravel().astype(np.int64)
        b2 = self.bu - np.asarray(self.Au[:, f1].sum(axis=1)).ravel().astype(np.int64)
        return self.Af[:, free].tocsr(), b1, self.Au[:, free].tocsr(), b2, free, Li

    def _solve(self, A, row_lo, row_hi, solver, tl, cost=None, col_lo=None, col_hi=None):
        hs = self.hs
        h = _new_highs(hs)
        h.setOptionValue("output_flag", False)
        h.setOptionValue("solver", solver)
        h.setOptionValue("time_limit", float(tl))
        if solver == "ipm":
            h.setOptionValue("run_crossover", "on" if cost is not None else "off")
        nv = A.shape[1]
        h.passModel(_highs_model(A, row_lo, row_hi,
                                 np.zeros(nv) if col_lo is None else col_lo,
                                 np.ones(nv) if col_hi is None else col_hi,
                                 np.zeros(nv) if cost is None else cost))
        h.run()
        return h

    def _certify_dead(self, A1, b1, A2, b2) -> bool:
        """Farkas multipliers from the LP  max b.Y - sum w  s.t.  A^T Y - w <= 0, w >= 0, -1 <= Ye <= 1, -1 <= Yu <= 0,
        then `certify` (exact).  True only with an exactly verified certificate."""
        me, mu, nv = A1.shape[0], A2.shape[0], A1.shape[1]
        At = sp.hstack([A1.T, A2.T, -sp.identity(nv)]).tocsc().astype(float)
        cost = -np.concatenate([b1.astype(float), b2.astype(float), -np.ones(nv)])
        clo = np.concatenate([-np.ones(me), -np.ones(mu), np.zeros(nv)])
        chi = np.concatenate([np.ones(me), np.zeros(mu), np.inf * np.ones(nv)])
        for solver in ("ipm", "simplex"):
            h = self._solve(At, -np.inf * np.ones(nv), np.zeros(nv), solver, self.FARKAS_TIME, cost=cost, col_lo=clo, col_hi=chi)
            if h.getModelStatus() != self.hs.HighsModelStatus.kOptimal:
                continue
            y = np.asarray(h.getSolution().col_value, float)
            for bits in (200, 60, 40):
                if certify(A1, b1, A2, b2, y[:me], y[me:me + mu], bits):
                    return True
        return False

    def dead(self, st):
        MS = self.hs.HighsModelStatus
        self.calls += 1
        lo, hi, rep = st
        L = lo[rep].astype(float)
        H = hi[rep].astype(float)
        if self._fits(L, H):
            return False
        ch = np.nonzero((L != self.curL) | (H != self.curH))[0]
        if len(ch):
            self.h.changeColsBounds(len(ch), ch.astype(np.int32), L[ch], H[ch])
            self.curL[ch] = L[ch]
            self.curH[ch] = H[ch]
        self.h.setOptionValue("time_limit", self.WARM_TIME)
        self.h.run()
        s = self.h.getModelStatus()
        if s == MS.kOptimal:
            self._remember(np.asarray(self.h.getSolution().col_value, float))
            return False
        A1, b1, A2, b2, free, Li = self._substituted(st)
        if A1.shape[1] == 0:
            return bool(np.any(b1 != 0) or np.any(b2 < 0))     # every column fixed: exact integer check
        if s != MS.kInfeasible:
            Ast, rlo, rhi = _stack(A1, b1, A2, b2)
            h = self._solve(Ast, rlo, rhi, "ipm", self.IPM_TIME)
            s2 = h.getModelStatus()
            if s2 == MS.kOptimal:
                x = Li.astype(float).copy()
                x[free] = np.asarray(h.getSolution().col_value, float)
                self._remember(x)
                return False
            if s2 != MS.kInfeasible:
                self.unknown += 1
                return None
        if self._certify_dead(A1, b1, A2, b2):
            return True
        self.unknown += 1
        return None
