"""Deterministic test families for the turnpike solver.

Every generator returns a list of `Case(name, points, distances, expected, note)`:
    points    : one realization (sorted, starts at 0) for YES cases, None for NO cases
    distances : the full multiset as a sorted list of n(n-1)/2 positive integers
    expected  : 'YES' or 'NO'
    note      : where the case comes from and, for NO cases, why it is certainly NO

NO cases carry an independent certificate, never the solver under test:
    parity  : for odd n the sum of all distances is even (sum_k x_k (2k - n + 1), every coefficient even), so a
              multiset with an odd sum has no realization;
    exhaust : the reference backtracking in bench/baseline.py finished its exhaustive search with no realization
              (that code is checked against brute-force enumeration of every point set with W <= 14).

    FAMILIES          name -> generator
    all_cases()       every case of every family
    load_corpus(fam)  the cases stored in tests/data/<fam>.json (written by `python tests/families.py --write`)
"""
from __future__ import annotations

import json
import math
import os
import random
import sys
from collections import Counter, namedtuple

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
DATA = os.path.join(HERE, "data")
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)

from turnpike.distances import canonical, distances_of  # noqa: E402

Case = namedtuple("Case", "name points distances expected note")


# ------------------------------------------------------------------------------------------------ helpers
def dlist(points) -> list:
    """Sorted distance multiset of `points` as a list."""
    c = distances_of(points)
    out = []
    for d in sorted(c):
        out.extend([d] * c[d])
    return out


def yes(name, points, note="") -> Case:
    pts = sorted(int(p) for p in points)
    assert len(set(pts)) == len(pts), name
    pts = [p - pts[0] for p in pts]
    return Case(name, pts, dlist(pts), "YES", note)


def product(A, B, M=None):
    """A + M*B with M > 2*max(A) so the copies cannot interact (a realization of A (x) B)."""
    A = sorted(A); B = sorted(B)
    if M is None:
        M = 2 * (A[-1] - A[0]) + 1
    return sorted({a + M * b for a in A for b in B})


def _baseline():
    bdir = os.path.join(ROOT, "bench")
    if bdir not in sys.path:
        sys.path.insert(0, bdir)
    import baseline
    return baseline


# ------------------------------------------------------------------------------------------------ YES families
def gen_random():
    out = []
    for n in (10, 30, 100, 300):
        for tag, W in (("sq", n * n), ("1e9", 10 ** 9)):
            rnd = random.Random(1000 * n + (1 if tag == "sq" else 2))
            pts = {0, W}
            while len(pts) < n:
                pts.add(rnd.randint(1, W - 1))
            out.append(yes(f"random_n{n}_W{tag}", pts, f"uniform random set, n={n}, W={W}, seed {1000 * n}"))
    return out


HOMOMETRIC_PAIRS = [
    ([0, 1, 4, 10, 12, 17], [0, 1, 8, 11, 13, 17]),
    ([0, 1, 2, 6, 8, 11], [0, 1, 6, 7, 9, 11]),        # found by exhaustive scan of W <= 14, n = 6
    ([0, 1, 2, 6, 10, 13], [0, 1, 4, 5, 11, 13]),
    ([0, 1, 3, 8, 10, 14], [0, 1, 7, 9, 11, 14]),
]


def homometric_products():
    """[(name, [realizations...])]: pairs and products A (x) B of pair members; every listed set has the same
    distance multiset (a product with a homometric factor stays homometric)."""
    out = []
    for i, (a, b) in enumerate(HOMOMETRIC_PAIRS):
        out.append((f"homometric_pair{i}", [a, b]))
    P0, P1 = HOMOMETRIC_PAIRS[0], HOMOMETRIC_PAIRS[1]
    M = 2 * 17 + 1
    out.append(("homometric_prod_00", [product(x, y, M) for x in P0 for y in P0]))
    out.append(("homometric_prod_01", [product(x, y, M) for x in P0 for y in P1]))
    out.append(("homometric_prod_10", [product(x, y, 23) for x in P1 for y in P0]))
    P2 = HOMOMETRIC_PAIRS[2]
    out.append(("homometric_prod_012", [product(product(x, y, M), z, M * 12)
                                        for x in P0 for y in P1 for z in P2]))
    return out


def gen_homometric():
    out = []
    for name, reals in homometric_products():
        d0 = dlist(reals[0])
        for r in reals[1:]:
            assert dlist(r) == d0, name
        k = len({tuple(canonical(r)) for r in reals})
        out.append(yes(name, reals[0], f"homometric: at least {k} realizations up to mirror"))
    return out


def zhang(nz, F, W):
    """Zhang's family (as quoted in Dakic 2000, Sec. 3.3): L = {0} u [1,nz] u [2nz+1,3nz] u F,
    R' = {0} u [1,nz] u [2nz+1,3nz] u G with F u G = [nz+1, 2nz] disjoint; X = L u (W - R')."""
    G = [t for t in range(nz + 1, 2 * nz + 1) if t not in F]
    base = [0] + list(range(1, nz + 1)) + list(range(2 * nz + 1, 3 * nz + 1))
    L = set(base) | set(F)
    R = set(base) | set(G)
    return sorted(L | {W - v for v in R})


def gen_zhang():
    out = []
    for nz in (3, 5, 7, 9, 11, 13):
        rnd = random.Random(77 + nz)
        masks = {
            "first": [nz + 1],
            "alt": [t for t in range(nz + 1, 2 * nz + 1) if (t - nz) % 2 == 1],
            "rand": sorted(rnd.sample(range(nz + 1, 2 * nz + 1), nz // 2)),
        }
        for mname, F in masks.items():
            X = zhang(nz, F, 1000 * nz)
            out.append(yes(f"zhang_nz{nz}_{mname}", X, f"Zhang family nz={nz}, F={F}, W={1000 * nz}"))
    return out


N20 = [0, 2, 4, 5, 6, 7, 8, 9, 12, 17, 18, 19, 22, 23, 26, 29, 31, 32, 41, 42]
N21 = [0, 1, 4, 5, 8, 9, 13, 17, 19, 20, 21, 23, 27, 29, 30, 31, 34, 39, 41, 44, 46]
HM_B = "11100011010110100001"


def hm(K, g, mu):
    """HM(K, g, mu): progression of step g on [0, K), 0/1 middle mu on [K, K+L), run on [W-K+1, W], W = 2K-1+L."""
    L = len(mu); W = 2 * K - 1 + L
    return sorted({0, W} | set(range(0, K, g)) | {K + j for j in range(L) if mu[j] == "1"}
                  | set(range(W - K + 1, W + 1)))


HM24B2 = hm(24, 3, HM_B * 2)


def gen_narrow():
    out = [yes("N20", N20, "hard narrow core N20"), yes("N21", N21, "hard narrow core N21"),
           yes("HM24B2", HM24B2, "HM(24,3,B^2), B=" + HM_B)]
    for K, r in ((12, 1), (18, 1), (24, 3), (30, 2)):
        out.append(yes(f"HM{K}B{r}", hm(K, 3, HM_B * r), f"HM({K},3,B^{r})"))
    return out


def hyperbolic(L):
    """Dot's hyperbolic (Fibonacci-height) chain, rebuilt from prose (lead T758); n = 5L + 9."""
    K = 2 * 3 ** L; Q = 10 * K + 1; g = (Q, Q * Q + 1); X = set()
    M = [[1, 0], [0, 1]]
    for _ in range(L + 1):
        gi = (g[0] * M[0][0] + g[1] * M[1][0], g[0] * M[0][1] + g[1] * M[1][1])
        X |= {gi[0] * a + gi[1] * b for a in range(3) for b in range(3)}
        M = [[M[0][0] + M[0][1], M[0][0] + 2 * M[0][1]], [M[1][0] + M[1][1], M[1][0] + 2 * M[1][1]]]
    return sorted(X)


def gen_rigid():
    out = [yes(f"doubling_L{L}", [0] + [2 ** j for j in range(L + 1)], f"doubling chain, W=2^{L}")
           for L in (8, 16, 32, 64)]
    for L in (2, 4, 8, 16):
        X = hyperbolic(L)
        assert len(X) == 5 * L + 9
        out.append(yes(f"hyperbolic_L{L}", X, f"Dot hyperbolic chain L={L}, n=5L+9"))
    return out


def symchain(core, L):
    core = sorted(core); core = [c - core[0] for c in core]; Wc = core[-1]
    t = 2 ** (L + 1); W = 2 * t + Wc
    return sorted({0, W} | {t + y for y in core} | {2 ** j for j in range(L + 1)} | {W - 2 ** j for j in range(L + 1)})


def gen_symchain():
    return [yes(f"symchain_HM24B2_L{L}", symchain(HM24B2, L),
                f"symmetric doubling chain around centred core HM24B2, L={L}") for L in (6, 10, 40)]


def mirrored_product(core, J, a=7):
    M = a * 2 ** J
    base = {M * x + e for x in core for e in (0, a)}
    W = max(base)
    return sorted(base | {a * 2 ** j for j in range(J + 1)} | {W - a * 2 ** j for j in range(J + 1)})


def gen_mirprod():
    return [yes(f"mirprod_N21_J{J}", mirrored_product(N21, J), f"7*2^{J}*N21 + {{0,7}} with mirrored chain 7*2^j")
            for J in (8, 10)]


def _primes_above(a, m):
    out = []; x = a
    while len(out) < m:
        x += 1
        if all(x % q for q in range(2, math.isqrt(x) + 1)):
            out.append(x)
    return out


def gen_skiena():
    out = []
    for n, seed in ((24, 1), (60, 2), (150, 3)):
        rnd = random.Random(seed)
        gaps = [100000 + rnd.randint(0, 17) for _ in range(n - 1)]
        X = [0]
        for v in gaps:
            X.append(X[-1] + v)
        out.append(yes(f"skiena_jitter_n{n}", X, f"gaps 100000+randint(0,17), seed {seed}"))
        gaps = _primes_above(10 ** 6, n - 1)
        random.Random(seed).shuffle(gaps)
        X = [0]
        for v in gaps:
            X.append(X[-1] + v)
        out.append(yes(f"skiena_primes_n{n}", X, f"gaps = shuffled consecutive primes above 1e6, seed {seed}"))
    return out


# ------------------------------------------------------------------------------------------------ NO families
def _certify_no(dists, n, cap=50_000):
    """Independent NO certificate or None."""
    if n % 2 == 1 and sum(dists) % 2 == 1:
        return "parity (odd n, odd distance sum)"
    st, _, nodes, _ = _baseline().backtrack(dists, node_cap=cap)
    if st == "NO":
        return f"exhaust (reference backtracking, {nodes} nodes)"
    if st == "YES":
        return "REALIZABLE"
    return None


def perturb(case: Case, rnd: random.Random, tries=6):
    """Move one copy of a distance d to d +- 1, keeping size; skip if the max repeats or a value becomes <= 0."""
    D = list(case.distances)
    W = D[-1]
    for _ in range(tries):
        i = rnd.randrange(len(D)); s = rnd.choice((-1, 1))
        d = D[i]; nd = d + s
        if nd <= 0:
            continue
        E = D[:i] + D[i + 1:] + [nd]
        E.sort()
        if E[-1] == E[-2]:
            continue
        cert = _certify_no(E, case_n(case))
        if cert and cert != "REALIZABLE":
            return E, f"{case.name}: one copy of {d} moved to {nd}; NO by {cert}"
    return None


def case_n(case):
    return (1 + math.isqrt(1 + 8 * len(case.distances))) // 2


def gen_no():
    rnd = random.Random(4242)
    out = []
    # odd n -> parity certificate whenever the move changes the parity (always: one value moves by 1);
    # even n -> exhaustive reference search with a node cap (dropped when the cap is hit)
    sources = (gen_random()[:4] + [c for c in gen_narrow() if c.name in ("N20", "N21", "HM24B2")]
               + gen_zhang()[:6] + [c for c in gen_zhang() if c.name in ("zhang_nz7_first", "zhang_nz9_first", "zhang_nz11_first", "zhang_nz13_first")] + gen_homometric()[:2] + gen_skiena()[:2] + gen_rigid()[4:6]
               + gen_symchain()[:1] + gen_mirprod()[:1])
    for c in sources:
        r = perturb(c, rnd)
        if r is not None:
            out.append(Case(f"no_perturb_{c.name}", None, r[0], "NO", r[1]))
    # random multisets of the right size (largest value unique)
    for n, W in ((7, 30), (10, 100), (11, 121), (12, 60), (31, 961), (101, 10 ** 4), (101, 10 ** 9)):
        for k in range(2):
            for _ in range(50):
                m = n * (n - 1) // 2
                E = sorted([rnd.randint(1, W - 1) for _ in range(m - 1)] + [W])
                cert = _certify_no(E, n)
                if cert and cert != "REALIZABLE":
                    out.append(Case(f"no_random_n{n}_W{W}_{k}", None, E, "NO", f"random multiset; NO by {cert}"))
                    break
    return out


FAMILIES = {
    "random": gen_random,
    "homometric": gen_homometric,
    "zhang": gen_zhang,
    "narrow": gen_narrow,
    "rigid": gen_rigid,
    "symchain": gen_symchain,
    "mirprod": gen_mirprod,
    "skiena": gen_skiena,
    "no": gen_no,
}


def all_cases():
    out = []
    for f in FAMILIES.values():
        out.extend(f())
    return out


# ------------------------------------------------------------------------------------------------ corpus I/O
def _to_json(c: Case) -> dict:
    n = (1 + math.isqrt(1 + 8 * len(c.distances))) // 2
    return {"name": c.name, "n": n, "W": c.distances[-1], "expected": c.expected, "note": c.note,
            "points": c.points, "distances": c.distances}


def write_corpus(families=None):
    os.makedirs(DATA, exist_ok=True)
    sizes = {}
    for fam in families or FAMILIES:
        cases = FAMILIES[fam]()
        path = os.path.join(DATA, f"{fam}.json")
        with open(path, "w") as fh:
            json.dump([_to_json(c) for c in cases], fh, separators=(",", ":"))
        sizes[fam] = (len(cases), os.path.getsize(path))
    return sizes


def load_corpus(family=None):
    """Cases from tests/data (all families, or one). Falls back to the generator when the file is missing."""
    fams = [family] if family else list(FAMILIES)
    out = []
    for fam in fams:
        path = os.path.join(DATA, f"{fam}.json")
        if os.path.exists(path):
            with open(path) as fh:
                for r in json.load(fh):
                    out.append(Case(r["name"], r["points"], r["distances"], r["expected"], r["note"]))
        else:
            out.extend(FAMILIES[fam]())
    return out


def family_of(name: str) -> str:
    for fam in FAMILIES:
        path = os.path.join(DATA, f"{fam}.json")
        if os.path.exists(path):
            with open(path) as fh:
                if any(r["name"] == name for r in json.load(fh)):
                    return fam
    raise KeyError(name)


if __name__ == "__main__":
    if "--write" in sys.argv:
        for fam, (k, b) in write_corpus().items():
            print(f"{fam:<12} {k:4d} cases  {b / 1e6:7.2f} MB")
    else:
        for c in all_cases():
            n = (1 + math.isqrt(1 + 8 * len(c.distances))) // 2
            print(f"{c.name:<32} n={n:<5} W={c.distances[-1]:<22} {c.expected}  {c.note}")
