"""Section 6: the T*_k search (propagation + exact face LP), polynomial under the conjecture (ST_k^LP+).

The search walks the T* class tree (candidate classes c = min(p, W - p) in increasing order, forced classes closed,
two children at each branching class).  Every child is tested by T*_k:

  * T*'s exact node test - counting, combinatorial oracle, degree peeling, exact 0/1 propagation on the SA_2 face
    system, and the face LP refuted only by an exactly checked integer Farkas certificate;
  * when two children pass that test: close_k^LP, k-level probing on the child's face system, with the exact LP as the
    refutation inside every leg (opts.tstar_k levels).

Theorem A-ST_k^LP: if (ST_k^LP+) holds - every child that does not extend to a realization, of a node that does, is
refuted by T*_k - then the DFS never backtracks and decides the list after at most 1 + 2(n - 1) child tests, each of
polynomial cost for fixed k.  The stage runs exactly that DFS with NO backtracking:

  YES        a realization was reached (re-verified by the orchestrator);
  NO         every child on the way was refuted by exact tests with no surviving sibling ever skipped: a proof;
  UNDECIDED  an entered child's subtree died while its node had another surviving child: the DFS would have to
             backtrack, i.e. (ST_k^LP+) failed on this list (work['st_violation'] = 1); the orchestrator falls back to
             the exhaustive search;
  TIMEOUT    the deadline expired.

Work counts: child_tests (children formed and tested), nodes (entered nodes), lp_calls (face LPs of T*'s test plus LP
decisions inside close_k), propagations (face propagations plus propagation calls inside close_k), close_tests,
lp_unknown (LP answers left undecided: they only prune less), st_violation.

opts.all_solutions: every surviving child is entered; all realizations reached are returned (up to mirror image).  A
surviving child whose subtree has no realization is still an (ST_k^LP+) violation.
"""
from __future__ import annotations

from collections import Counter

from . import StageResult, TimeUp
from ._tstar import Search, run_with_big_stack

NAME = "tstar"
SECTION = 6
COST = "POLY"


def _work(st, violation):
    return {
        "child_tests": int(st["child_tests"]),
        "nodes": int(st["nodes"]),
        "lp_calls": int(st["lp_calls"] + st["close_lp_calls"]),
        "propagations": int(st["propagations"] + st["close_propagations"]),
        "close_tests": int(st["close_tests"]),
        "close_refuted": int(st["close_refuted"]),
        "lp_unknown": int(st["lp_unknown"] + st["close_lp_unknown"]),
        "lp_cert_fail": int(st["lp_cert_fail"]),
        "st_violation": 1 if violation else 0,
    }


def run(inst, opts, ctx) -> StageResult:
    n = inst.n
    if n <= 1:
        return StageResult("YES", solutions=[[0]], work=_work(Counter(), False))
    deadline = ctx.get("deadline") if ctx else None
    k = int(getattr(opts, "tstar_k", 3))
    S = Search(inst.counts, n, k, deadline=deadline, all_solutions=bool(getattr(opts, "all_solutions", False)))
    try:
        status = run_with_big_stack(S.run)
    except TimeUp:
        return StageResult("TIMEOUT", work=_work(S.st, False), note="deadline reached in T*_%d" % k)
    if status == "VIOLATION":
        return StageResult("UNDECIDED", work=_work(S.st, True),
                           note="ST_k violation: a first-wrong child survived the test (k = %d)" % k)
    if status == "NO":
        return StageResult("NO", work=_work(S.st, False), note="T*_%d refuted every child" % k)
    sols, seen = [], set()
    for X in S.solutions:
        X = [int(p) - int(X[0]) for p in X]
        key = min(tuple(X), tuple(sorted(X[-1] - p for p in X)))
        if key not in seen:
            seen.add(key)
            sols.append(X)
    return StageResult("YES", solutions=sols, work=_work(S.st, False))
