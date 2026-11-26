"""Benchmark of the inputs that make the solver work hardest (larger and trickier than the test corpus).

    python bench/tricky.py [--timeout S] [--only SUBSTRING] [--no-mutants]

Each case runs in-process with a time limit; YES cases are followed by one or two no-solution variants (one distance
moved by 1). Prints the deciding section and the per-section seconds. Nothing here is used by the tests.
"""
from __future__ import annotations

import argparse
import os
import random
import sys
import time
from collections import Counter

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "tests"))

from families import HM_B, N20, N21, hm, mirrored_product, product, symchain, zhang  # noqa: E402
from turnpike import Instance, Options, solve  # noqa: E402
from turnpike.distances import distances_of  # noqa: E402


def cases():
    rnd = random.Random(5)
    out = []
    for nz in (15, 19, 25):
        out.append((f"zhang_nz{nz}_W{1000 * nz}", zhang(nz, [nz + 1], 1000 * nz)))
        out.append((f"zhang_nz{nz}_wide", zhang(nz, [nz + 1], 10 ** 9 + 7)))
    for K, r in ((36, 4), (42, 5), (48, 6)):
        out.append((f"HM{K}B{r}", hm(K, 3, HM_B * r)))
    out.append(("N21xN20_narrow", product(N21, N20)))
    out.append(("N21xN20_wide", product(N21, N20, 1_000_003)))
    out.append(("HM24B2xN20_wide", product(hm(24, 3, HM_B * 2), N20, 10 ** 7 + 19)))
    out.append(("mirprod_J16", mirrored_product(N21, 16)))
    out.append(("symchain_HM36_L40", symchain(hm(36, 3, HM_B * 4), 40)))
    for n, W in ((2000, 10 ** 9), (5000, 10 ** 12)):
        out.append((f"random_n{n}", sorted(rnd.sample(range(1, W), n - 2) + [0, W])))
    return out


def show(tag, r):
    secs = " ".join(f"{a.name}:{a.seconds:.2f}" for a in r.attempts[1:] if a.outcome != "SKIPPED")
    flag = "  ** FALLBACK **" if r.fallback_used else ""
    print(f"{tag:<34} n={r.n:<5} {r.status:<8} {str(r.decided_by):<10} {r.seconds:8.2f} s   {secs}{flag}", flush=True)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--timeout", type=float, default=120)
    ap.add_argument("--only", default="")
    ap.add_argument("--no-mutants", action="store_true")
    a = ap.parse_args()
    rnd = random.Random(1)
    for name, X in cases():
        if a.only and a.only not in name:
            continue
        inst = Instance.from_points(X)
        show(name, solve(inst, Options(time_limit=a.timeout)))
        if a.no_mutants:
            continue
        D = Counter(distances_of(X)) if len(X) < 1500 else None
        if D is None:
            continue
        W = max(D)
        ks = sorted(d for d in D if d != W)
        for t in range(2):
            d = rnd.choice(ks)
            c = Counter(D)
            c[d] -= 1
            c[d + 1] += 1
            if c[W] != 1 or d + 1 == W:
                continue
            show(f"  NO variant {t}", solve(dict(c), Options(time_limit=a.timeout)))


if __name__ == "__main__":
    main()
