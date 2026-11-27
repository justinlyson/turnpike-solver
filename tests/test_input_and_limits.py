"""Regression tests for input handling, time limits and complete solution lists."""
from __future__ import annotations

import time
from decimal import Decimal
from fractions import Fraction

import numpy as np
import pytest

from turnpike import Instance, Options, solve
from turnpike.__main__ import main
from turnpike.distances import InvalidInput, distances_of

SMALL = "1,4,10,12,17,3,9,11,16,6,8,13,2,7,5"            # points 0 1 4 10 12 17


# ------------------------------------------------------------------------------------------- non-integer input
@pytest.mark.parametrize("bad", [
    {1.5: 1, 2: 1, 3.5: 1},
    [np.float32(1.5), np.float32(2), np.float32(3.5)],
    [Fraction(3, 2), 2, Fraction(7, 2)],
    [True, 2, 3],
    ["1", "2", "3"],
])
def test_non_integer_distances_rejected(bad):
    with pytest.raises(InvalidInput):
        Instance.from_distances(bad)


def test_non_integer_points_rejected():
    with pytest.raises(InvalidInput):
        Instance.from_points([0, 1.5, 3])


@pytest.mark.parametrize("ok", [[1.0, 2, 3], [np.int64(1), 2, 3], [Decimal("1"), 2, 3], {1: 1, 2: 1, 3.0: 1}])
def test_integer_valued_input_accepted(ok):
    r = solve(ok)
    assert r.status == "YES" and r.points == [0, 1, 3]


# ------------------------------------------------------------------------------------------- command line
def test_cli_bad_tokens_exit_2():
    assert main(["-d", "1,2,abc"]) == 2
    assert main(["-d", "1,2,3x"]) == 2
    assert main(["-d", "1,2.5,3"]) == 2


def test_cli_missing_file_exit_2(tmp_path):
    assert main([str(tmp_path / "nope.txt")]) == 2


@pytest.mark.parametrize("enc", ["utf-8", "utf-8-sig", "utf-16"])
def test_cli_file_encodings(tmp_path, enc):
    f = tmp_path / "d.txt"
    f.write_bytes(SMALL.encode(enc))
    assert main([str(f)]) == 0


def test_cli_large_decimal_token_exact():
    # 12345678901234567.0 must not be rounded through float (it would become ...568)
    assert main(["-d", "12345678901234567.0,1,12345678901234568"]) == 0


# ------------------------------------------------------------------------------------------- time limit inside T*
def test_time_limit_honoured_in_tstar():
    X = sorted(u + 40 * i for u in (0, 1, 4, 10, 12, 17) for i in range(13))
    t = time.perf_counter()
    r = solve(distances_of(X), Options(quick_search=False, pseudo_poly=False, lift=False, fallback=False, time_limit=3))
    assert r.status in ("TIMEOUT", "YES")
    assert time.perf_counter() - t < 15


# ------------------------------------------------------------------------------------------- complete lists
def test_fallback_all_more_than_4096_hits():
    # the compiled search reaches the first realization 4096 times before the second one
    X = sorted(u + 40 * i for u in (0, 1, 4, 10, 12, 17) for i in range(13))
    o = Options(all_solutions=True, pseudo_poly=False, lift=False, tstar=False)
    r = solve(distances_of(X), o)
    assert r.decided_by == "fallback" and r.complete is True
    assert len(r.solutions) == len(solve(distances_of(X), Options(all_solutions=True)).solutions) == 2


def test_complete_flag():
    D = distances_of([0, 1, 4, 10, 12, 17])
    assert solve(D).complete is None
    assert solve(D, Options(all_solutions=True)).complete is True
    r = solve(D, Options(all_solutions=True, quick_search=False, pseudo_poly=False, lift=False, fallback=False))
    if r.decided_by not in ("factor", "lift"):
        assert r.complete is False and "MAY BE INCOMPLETE" in r.summary()
