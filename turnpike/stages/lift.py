"""Section 5: lift + factor (POLY when the lift height is poly(n); PSEUDO otherwise).

The lattice K(D).  Let S = supp D.  K(D) is the set of integer vectors v in Z^S with v(a) + v(b) = v(a + b)
whenever a, b, a + b are all in S (a = b allowed).  The identity v(d) = d lies in K(D).  Choose a reduced basis
k_1, ..., k_r of K(D) and put c(d) = (k_1(d), ..., k_r(d)) in Z^r.  If Y (min 0) realizes D, then
L(y) = c(y) (L(0) = 0) satisfies L(y') - L(y) = c(y' - y) for y < y' because y, y' - y, y' are all in S.  So every
realization lifts to Z^r and its 0/1 polynomial A_L divides the lifted polynomial
    P_L(u) = n + sum_d D(d) (u^{c(d)} + u^{-c(d)}),     P_L = A_L A_L^*,
and conversely every 0/1 factor decodes back through the table c(d) -> d.  c is injective (the identity is in K),
so enumerating the 0/1 factors of P_L (turnpike.polyfactor) gives ALL realizations; none proves NO.

Factoring.  Default (pseudo_poly): P_L is factored directly as a sparse polynomial in r variables (fmpz_mpoly;
fmpz_poly when r = 1).  Fastest on every family tried, but with no polynomial bound in r (dense size (2h+1)^r).
Strict mode (pseudo_poly off): the T3R route.  A random integer map pi: Z^r -> Z^3, injective on {+-c(d)} and never
0, keeps everything above exact (apply pi to L); the projected polynomial has 3 variables and degree poly(n, h).
work["mode"] records "direct" or "project3".

Computing K(D) fast.  Process S in increasing order.  A value s with no decomposition s = a + b (a, b in S) is
a free parameter; otherwise v(s) := v(a) + v(b) for its first decomposition.  This writes every v(s) as an integer
combination V[s] of the P free values (for a realizable D every non-gap distance decomposes, so P <= n - 1), and
K(D) is isomorphic to the integer kernel of the constraints V[a] + V[b] - V[c] = 0 over all relations.  That kernel
is found by iterated refinement: keep a rational basis N of the current solution space, evaluate all residuals at
once (vectorised: U = V N, residual = U[a] + U[b] - U[c]), and replace N by N * null(rows violated) until no
relation is violated.  Every refinement lowers dim N, the final check is exact, and saturation in Z^P plus LLL give
the reduced basis.  Height h = max |k_i(d)|.
"""
from __future__ import annotations

import random

import flint
import numpy as np

from ..distances import Instance
from ..options import Options
from ..polyfactor import enumerate_realizations, multivariate_P, univariate_P
from . import StageResult

NAME = "lift"
SECTION = 5
COST = "POLY"

INT64_SAFE = 2 ** 62
PROJECT_K = 3              # strict mode: T3R projection to this many coordinates (polynomial guarantee)


# ----------------------------------------------------------------------------------------------------------------
# 3-term relations
def _sweep_kernel():
    """numba two-pointer sweep (stages/_lift_kernel.py), or None when numba is unavailable."""
    try:
        from ._lift_kernel import three_term_nb
        return three_term_nb
    except Exception:                          # pragma: no cover
        return None


def three_term(S: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """All index triples (ia, ib, ic), ia <= ib, with S[ia] + S[ib] = S[ic] (in no particular order).  O(|S|^2):
    a compiled two-pointer sweep for large S, else vectorised binary searches."""
    M = len(S)
    W = int(S[-1])
    if W >= INT64_SAFE // 2:
        index = {int(s): i for i, s in enumerate(S)}
        trip = [(i, j, index[int(S[i]) + int(S[j])]) for i in range(M) for j in range(i, M)
                if int(S[i]) + int(S[j]) in index]
        arr = np.array(trip, dtype=np.int64).reshape(-1, 3)
    elif len(S) > 1500 and _sweep_kernel() is not None:
        arr = _sweep_kernel()(np.ascontiguousarray(S, dtype=np.int64))
    else:
        s = S.astype(np.int64)
        out = []
        half = int(np.searchsorted(s, W // 2, side="right"))       # a <= b and a + b <= W force a <= W/2
        chunk = max(1, 4_000_000 // max(M, 1))
        for i0 in range(0, half, chunk):
            i1 = min(half, i0 + chunk)
            j1 = int(np.searchsorted(s, W - s[i0], side="right"))  # columns with s[i0] + s[j] <= W
            sums = s[i0:i1, None] + s[None, i0:j1]
            pos = np.searchsorted(s, sums)
            np.minimum(pos, M - 1, out=pos)
            hit = (s[pos] == sums)
            hit &= np.arange(i0, j1)[None, :] >= np.arange(i0, i1)[:, None]
            ii, jj = np.nonzero(hit)
            out.append(np.stack([ii + i0, jj + i0, pos[ii, jj]], axis=1))
        arr = np.concatenate(out) if out else np.zeros((0, 3), dtype=np.int64)
    return arr[:, 0], arr[:, 1], arr[:, 2]


# ----------------------------------------------------------------------------------------------------------------
# exact lattice K(D)
def _parametrize(S, ia, ib, ic):
    """V (|S| x P): v(s) as an integer combination of the free values (values with no decomposition)."""
    M = len(S)
    first = np.full(M, -1, dtype=np.int64)
    if len(ic):
        first[ic[::-1]] = np.arange(len(ic) - 1, -1, -1)     # a decomposition of each value (the first listed)
    free = np.nonzero(first < 0)[0]
    P = len(free)
    big = int(S[-1]) // max(1, int(S[0])) >= INT64_SAFE // max(P, 1)     # entries <= W / min S (nonnegative)
    V = np.zeros((M, P), dtype=object if big else np.int64)
    V[free, np.arange(P)] = 1
    fa, fb = ia[first[first >= 0]], ib[first[first >= 0]]
    for s, a, b in zip(np.nonzero(first >= 0)[0], fa, fb):
        V[s] = V[a] + V[b]
    return V


def _matmul(A: np.ndarray, B: np.ndarray) -> np.ndarray:
    """Exact integer product (int64 when provably no overflow, else Python integers)."""
    if A.dtype != object and B.dtype != object and A.size and B.size:
        bound = int(np.abs(A).max()) * int(np.abs(B).max()) * A.shape[1]
        if bound < 2 ** 53:                    # float64 BLAS is exact below 2^53
            return np.rint(A.astype(np.float64) @ B.astype(np.float64)).astype(np.int64)
        if bound < INT64_SAFE:
            return A @ B
    return A.astype(object) @ B.astype(object)


def _to_np(m: flint.fmpz_mat) -> np.ndarray:
    rows = [[int(x) for x in r] for r in m.tolist()]
    big = any(abs(x) >= INT64_SAFE for r in rows for x in r)
    return np.array(rows, dtype=object if big else np.int64).reshape(m.nrows(), m.ncols())


def _lll_rows(B: np.ndarray) -> np.ndarray:
    return _to_np(flint.fmpz_mat([[int(x) for x in r] for r in B]).lll())


def _int_kernel(R: np.ndarray) -> np.ndarray:
    """Basis (columns, r x r') of {t in Z^r : R t = 0}, via the HNF of [R^T | I_r] (unimodular row operations:
    the rows whose R^T-part vanishes carry a basis of the kernel in the identity part).  The kernel is saturated."""
    m, r = R.shape
    H = flint.fmpz_mat([[int(R[j, i]) for j in range(m)] + [int(i == k) for k in range(r)] for i in range(r)]).hnf()
    rows = [[int(H[i, m + k]) for k in range(r)] for i in range(r) if all(H[i, j] == 0 for j in range(m))]
    if not rows:
        return np.zeros((r, 0), dtype=np.int64)
    return _lll_rows(np.array(rows, dtype=object)).T


def _kernel(V, ia, ib, ic, deadline):
    """Integer matrix N (P x r): a basis of {y in Z^P : V y satisfies every relation}, saturated in Z^P.

    Refinement: scan every relation with U = V N; collect violated residual rows (at most 8 r distinct ones);
    replace N by N * ker_Z(rows).  The intersection of a saturated lattice with a rational subspace is saturated,
    so N stays a saturated basis.  Stops after a full scan with no violation (exact)."""
    N = np.eye(V.shape[1], dtype=np.int64)
    if V.shape[1] == 1:                         # one free value: the identity spans K(D), so N = (1) (no scan)
        return N
    U = V
    T = len(ia)
    viol = None
    if T > 200_000:
        try:
            from ._lift_kernel import violations as viol
        except Exception:                      # pragma: no cover
            viol = None
    while N.shape[1] > 0:
        r = N.shape[1]
        found = []
        if viol is not None and U.dtype != object and int(np.abs(U).max()) < INT64_SAFE // 4:
            if deadline is not None and deadline.expired():
                return None
            idx = viol(np.ascontiguousarray(U), ia, ib, ic, 64 * r)
            if len(idx):
                found.append(np.unique(U[ia[idx]] + U[ib[idx]] - U[ic[idx]], axis=0))
        else:
            step = max(1, 4_000_000 // r)
            nfound = 0
            for lo in range(0, T, step):
                if deadline is not None and deadline.expired():
                    return None
                hi = min(T, lo + step)
                R = U[ia[lo:hi]] + U[ib[lo:hi]] - U[ic[lo:hi]]
                bad = R[(R != 0).any(axis=1)]
                if len(bad):
                    bad = np.unique(bad, axis=0) if bad.dtype != object else bad
                    found.append(bad)
                    nfound += len(bad)
                    if nfound >= 8 * r:
                        break
        if not found:
            break
        rows = np.concatenate(found)[: 8 * r]
        N = _matmul(N, _int_kernel(rows))
        if N.shape[1]:
            N = _lll_rows(N.T).T
        U = _matmul(V, N)
    return N


def lattice_basis(inst: Instance, deadline=None):
    """(S, C, work): S = supp D ascending, C (|S| x r) rows c(d) of an LLL-reduced basis of K(D); None on timeout."""
    S = np.array(inst.distinct, dtype=object if inst.W >= INT64_SAFE // 2 else np.int64)
    ia, ib, ic = three_term(S)
    V = _parametrize(S, ia, ib, ic)
    N = _kernel(V, ia, ib, ic, deadline)
    if N is None:
        return None
    K = _matmul(V, N).T
    K = _lll_rows(K)                                    # r x |S|, small entries
    C = np.ascontiguousarray(K.T)
    work = {"relations": int(len(ia)), "params": int(V.shape[1]), "rank": int(C.shape[1]),
            "height": int(np.abs(C).max()) if C.size else 0}
    return S, C, work


# ----------------------------------------------------------------------------------------------------------------
def project(C: np.ndarray, k: int, rng: random.Random) -> tuple[list, int]:
    """T3R projection: random integer pi: Z^r -> Z^k with entries in [-R, R] (R doubling from 1, a few tries per R)
    that is injective on {+-c(d)} and never 0.  Returns the projected keys pi(c(d)) and R."""
    Cobj = C.astype(object)
    R = 1
    while True:
        for _ in range(4):
            pi = np.array([[rng.randint(-R, R) for _ in range(k)] for _ in range(C.shape[1])], dtype=object)
            keys = [tuple(int(x) for x in v) for v in Cobj @ pi]
            allk = set(keys) | {tuple(-x for x in v) for v in keys}
            if len(allk) == 2 * len(keys) and (0,) * k not in allk:
                return keys, R
        R *= 2


# ----------------------------------------------------------------------------------------------------------------
def run(inst: Instance, opts: Options, ctx: dict) -> StageResult:
    n = inst.n
    if n <= 2:
        return StageResult("YES", solutions=[[0, inst.W] if n == 2 else [0]])
    deadline = ctx.get("deadline")
    lat = lattice_basis(inst, deadline)
    if lat is None:
        return StageResult("TIMEOUT", note="deadline while computing K(D)")
    S, C, work = lat
    r, h = work["rank"], work["height"]
    degree = 2 * h + 1                                  # per-variable degree of the lifted polynomial, plus one
    if not (h <= n ** opts.lift_max_height_exp or (opts.pseudo_poly and degree <= opts.lift_max_box)):
        return StageResult("UNDECIDED", work=work, note=f"skipped: height {h} > n^{opts.lift_max_height_exp}")
    if deadline is not None and deadline.expired():
        return StageResult("TIMEOUT", work=work)

    if opts.pseudo_poly or r <= PROJECT_K:
        keys = [tuple(int(x) for x in row) for row in C]
        work["mode"] = "direct" if opts.pseudo_poly else "project3"     # r <= 3: the projection is the identity
    else:
        keys, R = project(C, PROJECT_K, random.Random(opts.seed))
        r = PROJECT_K
        work.update({"mode": "project3", "proj_R": R, "proj_height": max(abs(x) for k in keys for x in k)})
    table = {}
    for key, d in zip(keys, S):
        table[key] = int(d)
        table[tuple(-x for x in key)] = -int(d)
    items = [(key, inst.counts[int(d)]) for key, d in zip(keys, S)]
    if r == 1:
        P, _ = univariate_P(n, [(key[0], c) for key, c in items])
        work["degree"] = P.degree()
        to_key = lambda e: (e,)
    else:
        P = multivariate_P(n, items, flint.fmpz_mpoly_ctx.get(tuple(f"u{i}" for i in range(r)), "lex"))
        work["terms"] = len(P)
        to_key = lambda e: e

    def decode(exps):
        """Lifted exponents -> points: every difference of two lifted points is +-c(d) for a distance d."""
        e0 = to_key(exps[0])
        ys = [0]
        for e in exps[1:]:
            d = table.get(tuple(a - b for a, b in zip(to_key(e), e0)))
            if d is None:
                return None
            ys.append(d)
        m = min(ys)
        return sorted(y - m for y in ys)

    res = enumerate_realizations(P, n, decode, inst, opts.all_solutions, deadline)
    work.update(res.work)
    if res.timed_out:
        return StageResult("TIMEOUT", work=work, note=res.work.get("why", "deadline"))
    if res.solutions:
        return StageResult("YES", solutions=res.solutions, work=work)
    return StageResult("NO", work=work, note="lifted polynomial has no 0/1 factor that decodes to a realization")
