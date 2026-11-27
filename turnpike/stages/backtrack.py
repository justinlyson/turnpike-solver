"""Classical backtracking (Skiena, Smith and Lemke, 1990), used twice:

  section 1  `quick`    : node budget max(64, n ** opts.quick_exp)  -> polynomial (POLY).
  section 7  `fallback` : no budget                                 -> exponential worst case (EXP).

Rule: the largest distance not yet explained must join one of the two end points, so the next point is y or W - y.
Place it, remove its distances to all placed points, recurse. With 0 and W placed first, the first choice y vs W - y
gives mirror images, so only y is tried at the root.

A run that finishes inside its budget is exact: it returns every realization (YES) or proves there is none (NO).
"""
from __future__ import annotations

import sys
import threading

import numpy as np

from . import StageResult


class _Stop(Exception):
    pass


# ---------------------------------------------------------------------------------------------- pure-Python search
def _search_py(counts: dict, n: int, W: int, node_limit, all_solutions: bool, deadline):
    cnt = dict(counts)
    cnt[W] -= 1
    keys = sorted((d for d in cnt), reverse=True)
    placed = [0, W]
    sols: list = []
    seen: set = set()                           # one point set is often reached along many paths: keep it once
    work = {"nodes": 0}
    limit = node_limit if node_limit is not None else -1

    def rec(ptr: int, root: bool) -> None:
        while ptr < len(keys) and cnt[keys[ptr]] == 0:
            ptr += 1
        if ptr == len(keys):
            key = tuple(sorted(placed))
            if key not in seen:
                seen.add(key)
                sols.append(list(key))
            if not all_solutions:
                raise _Stop
            return
        y = keys[ptr]
        cands = (y,) if (root or y == W - y) else (y, W - y)
        for p in cands:
            work["nodes"] += 1
            if work["nodes"] == limit:
                raise _Stop("limit")
            if deadline is not None and (work["nodes"] & 4095) == 0 and deadline.expired():
                raise _Stop("time")
            used = []
            ok = True
            for q in placed:
                d = p - q if p > q else q - p
                c = cnt.get(d, 0)
                if c == 0:
                    ok = False
                    break
                cnt[d] = c - 1
                used.append(d)
            if ok:
                placed.append(p)
                rec(ptr, False)
                placed.pop()
            for d in used:
                cnt[d] += 1

    box = {"reason": None, "error": None}

    def go():
        try:
            rec(0, True)
        except _Stop as e:
            box["reason"] = e.args[0] if e.args else None
        except BaseException as e:              # on the worker thread an error would otherwise vanish and read as NO
            box["error"] = e

    old = sys.getrecursionlimit()
    sys.setrecursionlimit(max(old, 4 * n + 1000))
    try:
        if n < 800:
            go()
        else:                                   # deep recursion: run on a thread with a large C stack
            old_ss = threading.stack_size()
            threading.stack_size(min(1 << 30, max(64 << 20, n * 64 << 10)))
            try:
                th = threading.Thread(target=go)
                th.start()
                th.join()
            finally:
                threading.stack_size(old_ss)
    finally:
        sys.setrecursionlimit(old)
    if box["error"] is not None:
        raise box["error"]
    return sols, work, box["reason"]


# ---------------------------------------------------------------------------------------------- compiled search
_KERNEL = None


def _kernel():
    """The numba kernel (stages/_bt_kernel.py), loaded on first use; None if numba is unavailable."""
    global _KERNEL
    if _KERNEL is None:
        try:
            from ._bt_kernel import _search_nb
            _KERNEL = _search_nb
        except Exception:                      # pragma: no cover
            _KERNEL = False
    return _KERNEL or None


def search(inst, node_limit=None, all_solutions=False, deadline=None, compiled=True):
    """Run the search. Returns (status, solutions, work) with status
    'YES' / 'NO' (exhaustive, exact), 'LIMIT' (node budget reached), 'TIMEOUT'."""
    n, W = inst.n, inst.W
    if n <= 2:
        return "YES", [[0, W]] if n == 2 else [[0]], {"nodes": 0}
    kern = _kernel() if (compiled and W < 2 ** 62) else None   # compiled runs cannot be interrupted by the clock
    if kern is not None:
        vals, cnt0 = inst.arrays()
        lim = -1 if node_limit is None else int(node_limit)
        st, nodes, sols = kern(vals, cnt0, n, W, lim, all_solutions, 4096 if all_solutions else 1)
        if st == 3:                             # more hits than the kernel's buffer holds: rerun uncapped
            sols, work, reason = _search_py(inst.counts, n, W, node_limit, all_solutions, deadline)
            work["compiled_overflow"] = 1
        else:
            sols = list({tuple(sorted(int(v) for v in row)): None for row in sols})
            sols = [list(t) for t in sols]
            work = {"nodes": int(nodes), "compiled": 1}
            if st == 2:
                return "LIMIT", sols, work
            return ("YES" if sols else "NO"), sols, work
    else:
        sols, work, reason = _search_py(inst.counts, n, W, node_limit, all_solutions, deadline)
    if reason == "limit":
        return "LIMIT", sols, work
    if reason == "time":
        return "TIMEOUT", sols, work
    return ("YES" if sols else "NO"), sols, work


# ---------------------------------------------------------------------------------------------- the two stages
class Quick:
    NAME = "quick"
    SECTION = 1
    COST = "POLY"

    @staticmethod
    def run(inst, opts, ctx) -> StageResult:
        budget = max(64, int(inst.n ** opts.quick_exp))
        # Phase 1: pure Python with a small budget (most inputs finish in about n nodes; loading the compiled kernel
        # costs ~1 s). Phase 2, only if needed: the compiled kernel with the full budget.
        # Below n = 150 the whole n^2 budget runs faster in Python than the kernel takes to load.
        small = budget if inst.n < 150 else min(budget, 20 * inst.n)
        st, sols, work = search(inst, small, opts.all_solutions, ctx.get("deadline"), compiled=inst.n >= 1000)
        if st == "LIMIT" and budget > small:
            st, sols, work = search(inst, budget, opts.all_solutions, ctx.get("deadline"), compiled=True)
        work["budget"] = budget
        if st == "LIMIT":
            return StageResult("UNDECIDED", work=work, note="node budget reached")
        if st == "TIMEOUT":
            return StageResult("TIMEOUT", work=work)
        if st == "YES" and opts.all_solutions:
            work["complete"] = 1          # exhaustive inside the budget: these are all realizations
        return StageResult(st, sols, work)


class Fallback:
    NAME = "fallback"
    SECTION = 7
    COST = "EXP"

    @staticmethod
    def run(inst, opts, ctx) -> StageResult:
        dl = ctx.get("deadline")
        st, sols, work = search(inst, opts.fallback_node_limit, opts.all_solutions, dl,
                                compiled=(dl is None or dl.t_end is None))
        if st == "LIMIT":
            return StageResult("UNDECIDED", work=work, note="fallback node limit reached")
        if st == "TIMEOUT":
            return StageResult("TIMEOUT", work=work)
        return StageResult(st, sols, work)
