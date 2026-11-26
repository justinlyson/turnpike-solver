"""Self-contained tests of the T*_k stage (no research code needed).

Hoped-for results: YES lists return a verified realization; one-move mutants return NO (or UNDECIDED, never a wrong
YES); n21f36 = N21 + {36} is an (ST_1^LP+) violation (UNDECIDED, st_violation = 1) and is decided at k = 2;
all_solutions finds both members of a homometric pair.
Run:  python tests/test_tstar_port.py      (also collectable by pytest)
"""
from __future__ import annotations

import os
import random
import sys
from collections import Counter

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from turnpike.distances import Instance, canonical, distances_of, verify  # noqa: E402
from turnpike.options import Options  # noqa: E402
from turnpike.stages import Deadline, tstar  # noqa: E402

N21 = [0, 1, 4, 5, 8, 9, 13, 17, 19, 20, 21, 23, 27, 29, 30, 31, 34, 39, 41, 44, 46]
N20 = [0, 2, 4, 5, 6, 7, 8, 9, 12, 17, 18, 19, 22, 23, 26, 29, 31, 32, 41, 42]


def _run(inst, k=2, all_solutions=False):
    o = Options()
    o.tstar_k = k
    o.all_solutions = all_solutions
    return tstar.run(inst, o, {"deadline": Deadline(600)})


def _yes_lists():
    rnd = random.Random(7)
    out = [N21, N20, [0, 5, 13, 14, 22, 27, 39, 40, 51]]
    for _ in range(4):
        n = rnd.randint(10, 22)
        W = rnd.randint(2 * n, 6 * n)
        out.append(sorted({0, W} | set(rnd.sample(range(1, W), n - 2))))
    return out


def _mutant(X, rnd):
    D = distances_of(X)
    vals = sorted(D)
    while True:
        d = rnd.choice(vals[:-1])
        if d + 1 != vals[-1]:
            D[d] -= 1
            D[d + 1] += 1
            return Instance.from_distances(+D)


def test_yes():
    for X in _yes_lists():
        inst = Instance.from_points(X)
        r = _run(inst)
        assert r.status == "YES", (X, r.status, r.note)
        assert verify(r.solutions[0], inst)
        assert r.work["child_tests"] <= 1 + 2 * (inst.n - 1)


def test_no_mutants():
    rnd = random.Random(11)
    for X in _yes_lists():
        inst = _mutant(X, rnd)
        r = _run(inst)
        assert r.status in ("NO", "UNDECIDED"), (X, r.status)


def test_n21f36():
    inst = Instance.from_points(N21 + [36])
    r1 = _run(inst, k=1)
    assert r1.status == "UNDECIDED" and r1.work["st_violation"] == 1
    r2 = _run(inst, k=2)
    assert r2.status == "YES" and canonical(r2.solutions[0]) == canonical(N21 + [36])


def test_all_solutions():
    r = _run(Instance.from_points([0, 1, 4, 10, 12, 17]), all_solutions=True)
    assert r.status == "YES"
    assert sorted(canonical(X) for X in r.solutions) == sorted([canonical([0, 1, 4, 10, 12, 17]),
                                                                canonical([0, 1, 8, 11, 13, 17])])


if __name__ == "__main__":
    for f in (test_yes, test_no_mutants, test_n21f36, test_all_solutions):
        f()
        print(f.__name__, "ok", flush=True)
