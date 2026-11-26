"""Correctness tests for turnpike.solve on the corpus in tests/data (see tests/families.py).

    pytest tests/test_solver.py                 # default: every case except the slow list
    TURNPIKE_SLOW=1 pytest tests/test_solver.py # everything

A wrong answer is a failure everywhere. UNDECIDED / TIMEOUT are allowed only in the strict-mode and time-limit tests.
"""
from __future__ import annotations

import os
import sys

import pytest

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
sys.path.insert(0, HERE)

import families  # noqa: E402
import turnpike  # noqa: E402
from turnpike.distances import Instance, canonical, distances_of, verify  # noqa: E402
from turnpike.options import Options  # noqa: E402

SLOW = os.environ.get("TURNPIKE_SLOW") == "1"
TIME = float(os.environ.get("TURNPIKE_TEST_TIMEOUT", "120"))
SLOW_CASES = {"random_n300_Wsq", "random_n300_W1e9", "homometric_prod_012", "symchain_HM24B2_L40"}
# solver took > 15 s on these in bench/results.md (2026-11-26, section 3 'factor'); drop from this set once fast
SLOW_CASES |= {"no_perturb_zhang_nz3_first", "no_perturb_zhang_nz3_rand"}

CASES = families.load_corpus()
YES = [c for c in CASES if c.expected == "YES"]
NO = [c for c in CASES if c.expected == "NO"]


def _param(cases):
    out = []
    for c in cases:
        marks = [pytest.mark.skip(reason="slow; set TURNPIKE_SLOW=1")] if (c.name in SLOW_CASES and not SLOW) else []
        out.append(pytest.param(c, id=c.name, marks=marks))
    return out


def _multiset(lst):
    m = {}
    for d in lst:
        m[d] = m.get(d, 0) + 1
    return m


def _check_yes_points(case, pts):
    inst = Instance.from_distances(case.distances)
    assert pts is not None, f"{case.name}: YES without points"
    assert verify(pts, inst), f"{case.name}: returned points do not verify"
    assert dict(distances_of(pts)) == _multiset(case.distances), f"{case.name}: wrong distance multiset"


# ------------------------------------------------------------------------------------------- default mode
@pytest.mark.parametrize("case", _param(YES))
def test_yes(case):
    rep = turnpike.solve(case.distances, Options(time_limit=TIME))
    assert rep.status == "YES", f"{case.name}: status {rep.status}, decided_by {rep.decided_by}\n{rep.summary()}"
    _check_yes_points(case, rep.points)
    assert list(rep.points) == canonical(rep.points), f"{case.name}: points not in canonical form"
    assert rep.n == families.case_n(case) and rep.W == case.distances[-1]


@pytest.mark.parametrize("case", _param(NO))
def test_no(case):
    rep = turnpike.solve(case.distances, Options(time_limit=TIME))
    assert rep.status == "NO", f"{case.name} ({case.note}): status {rep.status}\n{rep.summary()}"
    assert rep.points is None


def test_instance_input_accepted():
    c = YES[0]
    rep = turnpike.solve(Instance.from_distances(c.distances))
    assert rep.status == "YES"
    _check_yes_points(c, rep.points)


# ------------------------------------------------------------------------------------------- all solutions
@pytest.mark.parametrize("name,reals", [pytest.param(nm, r, id=nm, marks=(
    [pytest.mark.skip(reason="slow")] if nm in SLOW_CASES and not SLOW else []))
    for nm, r in families.homometric_products()])
def test_homometric_all_solutions(name, reals):
    D = families.dlist(reals[0])
    rep = turnpike.solve(D, Options(all_solutions=True, time_limit=TIME))
    assert rep.status == "YES", f"{name}: {rep.status}"
    inst = Instance.from_distances(D)
    got = {tuple(canonical(s)) for s in rep.solutions}
    for s in rep.solutions:
        assert verify(s, inst), f"{name}: a returned solution does not verify"
    assert len(got) == len(rep.solutions), f"{name}: duplicate solutions (up to mirror)"
    want = {tuple(canonical(r)) for r in reals}
    missing = want - got
    assert not missing, f"{name}: missing realizations {sorted(missing)[:2]}"


def test_homometric_pair_both():
    a, b = families.HOMOMETRIC_PAIRS[0]
    rep = turnpike.solve(families.dlist(a), Options(all_solutions=True))
    got = {tuple(canonical(s)) for s in rep.solutions}
    assert {tuple(canonical(a)), tuple(canonical(b))} <= got


# ------------------------------------------------------------------------------------------- strict mode
@pytest.mark.parametrize("case", _param(CASES))
def test_strict_never_wrong(case):
    rep = turnpike.solve(case.distances, Options(time_limit=TIME).strict())
    assert rep.status in (case.expected, "UNDECIDED", "TIMEOUT"), \
        f"{case.name}: strict mode WRONG ANSWER {rep.status} (expected {case.expected})"
    assert not rep.fallback_used, f"{case.name}: strict mode ran the exponential fallback"
    if rep.status == "YES":
        _check_yes_points(case, rep.points)


# ------------------------------------------------------------------------------------------- time limits
HARD = ["zhang_nz13_alt", "HM24B3", "HM30B2", "symchain_HM24B2_L40", "random_n300_W1e9", "homometric_prod_012",
        "no_random_n101_W1000000000_0", "no_perturb_mirprod_N21_J8"]


@pytest.mark.parametrize("limit", [0.001, 0.05, 0.5])
@pytest.mark.parametrize("name", HARD)
def test_time_limit_never_wrong(name, limit):
    case = next(c for c in CASES if c.name == name)
    rep = turnpike.solve(case.distances, Options(time_limit=limit))
    assert rep.status in (case.expected, "TIMEOUT"), f"{name} @ {limit}s: {rep.status} (expected {case.expected})"
    if rep.status == "YES":
        _check_yes_points(case, rep.points)
    # the limit is honoured (generous slack for setup and one indivisible step)
    assert rep.seconds <= limit + 10.0, f"{name}: ran {rep.seconds:.1f}s with a {limit}s limit"


# ------------------------------------------------------------------------------------------- corpus integrity
@pytest.mark.parametrize("fam", [f for f in families.FAMILIES if f != "no"] + ["no"])
def test_corpus_matches_generators(fam):
    if fam == "no" and not SLOW:
        pytest.skip("regenerating NO certificates takes ~40 s; set TURNPIKE_SLOW=1")
    stored = families.load_corpus(fam)
    fresh = families.FAMILIES[fam]()
    assert [(c.name, c.points, c.distances, c.expected) for c in stored] == \
           [(c.name, c.points, c.distances, c.expected) for c in fresh]
    for c in stored:
        if c.expected == "YES":
            assert families.dlist(c.points) == c.distances, c.name
