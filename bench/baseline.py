"""Reference baseline: plain classical backtracking (Skiena, Smith and Lemke 1990), written independently of the package.

Rule: the largest distance not yet explained joins one of the two end points, so the next point is y or W - y.
With 0 and W placed, the first choice is fixed to y (the other choice gives the mirror image).
No pruning beyond the multiset check. A run that ends inside its caps is exact (YES with a realization, or NO).

    backtrack(distances, node_cap=None, time_cap=None, all_solutions=False) -> (status, points, nodes, solutions)
    status: 'YES' | 'NO' | 'CAP'   (CAP = node or time cap hit before a decision)
"""
from __future__ import annotations

import heapq
import math
import time
from collections import Counter


class _Cap(Exception):
    pass


def backtrack(distances, node_cap=None, time_cap=None, all_solutions=False):
    if isinstance(distances, dict):
        D = Counter({int(d): int(m) for d, m in distances.items() if m})
    else:
        D = Counter(int(d) for d in distances)
    total = sum(D.values())
    n = (1 + math.isqrt(1 + 8 * total)) // 2
    if n * (n - 1) // 2 != total or any(d <= 0 for d in D):
        return "NO", None, 0, []
    if n <= 1:
        return "YES", [0], 0, [[0]]
    W = max(D)
    if D[W] != 1:
        return "NO", None, 0, []
    D[W] -= 1
    heap = [-d for d in D if D[d] > 0]
    heapq.heapify(heap)
    pts = [0, W]
    sols = []
    nodes = 0
    t_end = None if time_cap is None else time.perf_counter() + time_cap

    def top():
        while heap and D[-heap[0]] <= 0:
            heapq.heappop(heap)
        return -heap[0] if heap else None

    def try_place(c):
        need = Counter()
        for p in pts:
            d = c - p if c > p else p - c
            if d == 0:
                return None
            need[d] += 1
            if need[d] > D[d]:
                return None
        return need

    def rec(root):
        nonlocal nodes
        nodes += 1
        if node_cap is not None and nodes > node_cap:
            raise _Cap
        if t_end is not None and (nodes & 1023) == 0 and time.perf_counter() > t_end:
            raise _Cap
        if len(pts) == n:
            sols.append(sorted(pts))
            return not all_solutions
        y = top()
        cands = [y] if root or y == W - y else [y, W - y]
        for c in cands:
            need = try_place(c)
            if need is None:
                continue
            for d, m in need.items():
                D[d] -= m
            pts.append(c)
            stop = rec(False)
            pts.pop()
            for d, m in need.items():
                if D[d] == 0:
                    heapq.heappush(heap, -d)
                D[d] += m
            if stop:
                return True
        return False

    try:
        rec(True)
    except _Cap:
        if sols and not all_solutions:
            return "YES", sols[0], nodes, sols
        return "CAP", (sols[0] if sols else None), nodes, sols
    if sols:
        return "YES", sols[0], nodes, sols
    return "NO", None, nodes, []
