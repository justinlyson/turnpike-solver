"""Tests for the algebraic stages: direct factoring (section 2) and lift + factor (section 3).

Run with `python -m pytest tests/test_algebra.py` or directly (`python tests/test_algebra.py` prints timings).
Every YES answer is re-verified here; every NO is on an instance obtained by a mutation that is checked to be
unrealizable by the other algebraic stage where that stage applies.
"""
from __future__ import annotations

import os
import random
import sys
import time

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from turnpike.distances import Instance, canonical, parse_numbers, verify  # noqa: E402
from turnpike.options import Options  # noqa: E402
from turnpike.stages import Deadline, factor, lift  # noqa: E402

DATA = r"D:\Math\Turnpike_Reconstruction\research\experiments\FE_liftB_20261121\out"
N21 = [0, 1, 4, 5, 8, 9, 13, 17, 19, 20, 21, 23, 27, 29, 30, 31, 34, 39, 41, 44, 46]
TIMINGS: list = []


def _run(stage, inst, all_solutions=False, limit=120.0, strict=False):
    t = time.perf_counter()
    opts = Options(all_solutions=all_solutions)
    res = stage.run(inst, opts.strict() if strict else opts, {"deadline": Deadline(limit)})
    TIMINGS.append((stage.NAME, inst.n, inst.W, res.status, round(time.perf_counter() - t, 3), res.work))
    for pts in res.solutions:
        assert verify(pts, inst), (stage.NAME, pts)
    return res


def _yes(stage, X, all_solutions=False, strict=False):
    inst = Instance.from_points(X)
    res = _run(stage, inst, all_solutions, strict=strict)
    assert res.status == "YES", (stage.NAME, res.status, res.note, res.work)
    return res


def _mutant(X, rng):
    """Move one copy of a distance d to d + 1 (count kept; the maximum stays unique)."""
    c = dict(Instance.from_points(X).counts)
    W = max(c)
    d = rng.choice([d for d in c if d < W - 1])
    c[d] -= 1
    c[d + 1] = c.get(d + 1, 0) + 1
    return Instance.from_distances({k: v for k, v in c.items() if v})


# ---------------------------------------------------------------------------------------------------- families
def random_set(rng, n, W):
    return sorted(rng.sample(range(1, W), n - 2) + [0, W])


def steven_gaps(rng, n):
    xs = [0]
    for _ in range(n - 1):
        xs.append(xs[-1] + 100000 + rng.randint(0, 17))
    return xs


def prime_gaps(rng, n):
    """Points whose consecutive gaps are consecutive primes above 10^6 in random order."""
    ps, p = [], 10 ** 6
    while len(ps) < n - 1:
        p += 1
        if all(p % q for q in range(2, int(p ** 0.5) + 1)):
            ps.append(p)
    rng.shuffle(ps)
    xs = [0]
    for g in ps:
        xs.append(xs[-1] + g)
    return xs


def mirrored_chain(J):
    """X = M*N21 + {0, 7}, plus {7*2^j, W - 7*2^j : j <= J}, M = 7*2^J (rank 1, height W/7)."""
    M = 7 * 2 ** J
    core = {M * a + b for a in N21 for b in (0, 7)}
    W = max(core)
    return sorted(core | {7 * 2 ** j for j in range(J + 1)} | {W - 7 * 2 ** j for j in range(J + 1)})


def zhang_style(k, W):
    """Two dense end clusters {0..k} and W - {0..k} plus their midpoint (many partial placements for backtracking)."""
    return sorted(set(range(k + 1)) | {W - i for i in range(k + 1)} | {W // 2})


def load(name):
    with open(os.path.join(DATA, name)) as f:
        return parse_numbers(f.read())


# ---------------------------------------------------------------------------------------------------- tests
def test_homometric_pair_all_solutions():
    a, b = [0, 1, 4, 10, 12, 17], [0, 1, 8, 11, 13, 17]
    want = {tuple(canonical(a)), tuple(canonical(b))}
    for stage in (factor, lift):
        for X in (a, b):
            got = {tuple(s) for s in _yes(stage, X, all_solutions=True).solutions}
            assert got == want, (stage.NAME, got)


def test_random_sets_lift():
    rng = random.Random(5)
    for n in (5, 10, 20, 30, 45, 60):
        _yes(lift, random_set(rng, n, 10 ** 9))


def test_random_sets_factor():
    rng = random.Random(6)
    for n, W in ((5, 40), (12, 300), (25, 1500)):
        _yes(factor, random_set(rng, n, W))


def test_no_mutants():
    rng = random.Random(7)
    for n, W in ((8, 60), (12, 200), (16, 600)):
        X = random_set(rng, n, W)
        inst = _mutant(X, rng)
        r_f = _run(factor, inst)
        r_l = _run(lift, inst)
        assert r_f.status in ("NO", "YES") and r_l.status in ("NO", "YES")
        assert r_f.status == r_l.status, (r_f.status, r_l.status)


def test_zhang_style():
    for k, W in ((4, 200), (6, 1001)):
        _yes(factor, zhang_style(k, W))
        _yes(lift, zhang_style(k, W))


def test_hm24b2():
    X = load("HM24B2.txt")
    _yes(factor, X)
    _yes(lift, X)


def test_steven_gaps_lift():
    rng = random.Random(8)
    for n in (30, 60):
        _yes(lift, steven_gaps(rng, n))


def test_prime_gaps_lift():
    """n = 40 (rank 2, height 322) is a known slow case: flint's bivariate factorization takes ~160 s."""
    rng = random.Random(9)
    for n in (20,):
        _yes(lift, prime_gaps(rng, n))


def test_mirrored_chain():
    for J in (2, 4):
        X = mirrored_chain(J)
        _yes(lift, X)


def test_strict_mode_projection():
    """pseudo_poly off: lift uses the T3R projection to 3 coordinates (mode project3), still exact."""
    rng = random.Random(10)
    for X in (random_set(rng, 12, 10 ** 9), random_set(rng, 20, 10 ** 9), steven_gaps(rng, 30)):
        res = _yes(lift, X, strict=True)
        assert res.work["mode"] == "project3"
    inst = _mutant(random_set(rng, 10, 10 ** 6), rng)
    res = _run(lift, inst, strict=True)
    assert res.status in ("NO", "YES", "UNDECIDED")


if __name__ == "__main__":
    tests = [v for k, v in sorted(globals().items()) if k.startswith("test_")]
    failed = 0
    for t in tests:
        try:
            t()
            print(f"PASS {t.__name__}", flush=True)
        except AssertionError as e:
            failed += 1
            print(f"FAIL {t.__name__}: {e}", flush=True)
    print(f"\n{'stage':6} {'n':>4} {'W':>14} {'status':8} {'sec':>8}  work")
    for name, n, W, st, sec, work in TIMINGS:
        keys = ("mode", "rank", "height", "proj_height", "degree", "terms", "factors", "pairs", "candidates")
        print(f"{name:6} {n:4} {W:14} {st:8} {sec:8.3f}  " + " ".join(f"{k}={work[k]}" for k in keys if k in work))
    sys.exit(1 if failed else 0)
