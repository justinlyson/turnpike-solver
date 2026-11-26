"""Section 3: anchors + square root ("Vieta" step). Polynomial (POLY).

Idea. Points that come in mirror pairs {p, W - p} can be placed without knowing which side anything else is on.
Call them anchors A (always including 0 and W). Every other point q must satisfy |q - a| in D for all a in A, and because
A is mirror-symmetric, q and W - q pass this test together. So the test gives the folded set of pairs {q, W - q}; what
is left unknown is which side of each pair is used.

Write the unknown side choice as a polynomial b(z) = sum of z^(2q - W) over the chosen q (centred, doubled exponents).
Then b + b* = S (all folded pairs, known) and b b* = P (the autocorrelation of the non-anchor points, known from D once
the anchor distances are removed). So b and b* are the two roots of t^2 - S t + P: b = (S + sqrt(S^2 - 4P)) / 2.
The square root of a sparse polynomial is computed term by term from the top (each step fixes the next term), exactly
and in polynomial time. Pairs used on both sides form a known symmetric part F, handled with b b* = P - F^2 - F S.

Anchor finding here is the "symmetric peel": starting from {0, W}, repeatedly take the largest unexplained distance y
and add the pair {y, W - y} if all of its distances are available. This is a heuristic (the open Anchor Lemma is about
when it suffices); the stage is sound regardless, because every answer is verified and failure gives UNDECIDED.
"""
from __future__ import annotations

from collections import Counter

from . import StageResult

NAME = "vieta"
SECTION = 3
COST = "POLY"


def _mul(a: Counter, b: Counter) -> Counter:
    out: Counter = Counter()
    for i, x in a.items():
        for j, y in b.items():
            out[i + j] += x * y
    return Counter({k: v for k, v in out.items() if v})


def sparse_sqrt(q: Counter, max_terms: int):
    """Exact square root R of a sparse Laurent polynomial q with leading coefficient 1, computed from the top term down.
    Returns None if q is not a perfect square with at most max_terms terms."""
    if not q:
        return Counter()
    top_e = max(q)
    if q[top_e] != 1 or top_e % 2:
        return None
    top = top_e // 2
    r = Counter({top: 1})
    res = Counter(q)
    res[top_e] -= 1
    res = Counter({k: v for k, v in res.items() if v})
    while res:
        if len(r) > max_terms:
            return None
        t = max(res)
        k = t - top
        if k >= top:                          # would not be a lower term: not a square
            return None
        coef = res[t]
        if coef % 2:
            return None
        c = coef // 2
        res[2 * k] -= c * c
        for e, v in r.items():
            res[e + k] -= 2 * c * v
        res = Counter({e: v for e, v in res.items() if v})
        r[k] += c
        if r[k] == 0:
            del r[k]
    return r


def _peel_anchors(cnt: Counter, W: int):
    """Symmetric peel. Returns (anchors sorted, remaining counts)."""
    A = [0, W]
    cnt = Counter(cnt)
    cnt[W] -= 1
    keys = sorted((d for d in cnt if cnt[d] > 0), reverse=True)
    idx = 0
    while True:
        while idx < len(keys) and cnt[keys[idx]] <= 0:
            idx += 1
        if idx == len(keys):
            break
        y = keys[idx]
        p, q = y, W - y
        if p == q or p in A or q in A:
            break
        need = Counter()
        for a in A:
            need[abs(p - a)] += 1
            need[abs(q - a)] += 1
        need[abs(p - q)] += 1
        if any(cnt[d] < m for d, m in need.items()):
            break
        for d, m in need.items():
            cnt[d] -= m
        A += [p, q]
    return sorted(A), cnt


def run(inst, opts, ctx) -> StageResult:
    W, n = inst.W, inst.n
    D = Counter(inst.counts)
    A, _ = _peel_anchors(D, W)
    work = {"anchors": len(A)}
    if len(A) <= 2 or len(A) >= n:
        return StageResult("UNDECIDED", work=work, note="no useful anchor set")
    Aset = set(A)
    # 1. folded candidates
    cand = set()
    for d in D:
        for q in (d, W - d):
            if 0 < q < W and q not in Aset and all(D.get(abs(q - a), 0) > 0 for a in A):
                cand.add(q)
    pairs = sorted({min(q, W - q) for q in cand if W - q in cand})
    work["folded_pairs"] = len(pairs)
    # 2. remove anchor-anchor and anchor-core distances (side-independent because A is mirror-symmetric)
    dc = Counter(D)
    for j in range(len(A)):
        for i in range(j):
            dc[A[j] - A[i]] -= 1
    kp = {q: min(2, min(dc.get(abs(q - a), 0) for a in A)) for q in pairs}
    for q in pairs:
        for a in A:
            dc[abs(q - a)] -= kp[q]
    if min(dc.values(), default=0) < 0 or sum(kp.values()) + len(A) != n:
        return StageResult("UNDECIDED", work=work, note="anchor accounting failed")
    # 3. S (free pairs), F (pairs used on both sides), P (autocorrelation of non-anchor points), doubled exponents
    S, F = Counter(), Counter()
    for q in pairs:
        u = W - 2 * q
        if kp[q] == 1:
            S[u] += 1
            S[-u] += 1
        elif kp[q] == 2:
            F[u] += 1
            F[-u] += 1
    P = Counter({0: sum(kp.values())})
    for d, k in dc.items():
        if k > 0:
            P[2 * d] += k
            P[-2 * d] += k
    for e, v in _mul(F, F).items():
        P[e] -= v
    for e, v in _mul(F, S).items():
        P[e] -= v
    Q = _mul(S, S)
    for e, v in P.items():
        Q[e] -= 4 * v
    Q = Counter({k: v for k, v in Q.items() if v})
    R = sparse_sqrt(Q, 2 * len(pairs) + 2) if Q else Counter()
    if R is None:
        return StageResult("UNDECIDED", work=work, note="square root failed")
    two_b = Counter(S)
    for e, v in R.items():
        two_b[e] += v
    if any(v not in (0, 2) for v in two_b.values()):
        return StageResult("UNDECIDED", work=work, note="root is not 0/1")
    pts = sorted(set(A) | {(e + W) // 2 for e, v in two_b.items() if v} | {(e + W) // 2 for e in F})
    work["sqrt_terms"] = len(R)
    return StageResult("YES", [pts], work)
