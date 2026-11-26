"""Command line: python -m turnpike [FILE] [options]

Input: integers separated by spaces, commas or newlines, from FILE, from --distances, or from standard input.
By default they are the n(n-1)/2 distances; with --points they are points (the distances are computed first; useful
for testing).
"""
from __future__ import annotations

import argparse
import json
import sys

from .distances import InvalidInput, Instance, parse_numbers
from .options import Options
from .solver import solve


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(prog="turnpike", description="Reconstruct points on a line from pairwise distances.")
    ap.add_argument("file", nargs="?", help="input file (default: standard input)")
    ap.add_argument("--distances", "-d", help="the input numbers inline, e.g. '1,2,3'")
    ap.add_argument("--points", action="store_true", help="input numbers are points, not distances")
    ap.add_argument("--all", action="store_true", help="list every realization (up to mirror image)")
    ap.add_argument("--json", action="store_true", help="print the report as JSON")
    ap.add_argument("--strict", action="store_true",
                    help="polynomial sections only: no pseudo-polynomial factoring, no exhaustive fallback")
    ap.add_argument("--no-quick", action="store_true", help="skip budgeted backtracking")
    ap.add_argument("--no-pseudo", action="store_true", help="skip pseudo-polynomial tricks (direct factoring, large lifts)")
    ap.add_argument("--no-lift", action="store_true", help="skip lift + factor")
    ap.add_argument("--no-positivity", action="store_true", help="skip the positivity certificate (proves NO)")
    ap.add_argument("--no-vieta", action="store_true", help="skip anchors + square root")
    ap.add_argument("--no-tstar", action="store_true", help="skip propagation + exact LP (T*)")
    ap.add_argument("--no-fallback", action="store_true", help="skip the exhaustive fallback")
    ap.add_argument("-k", type=int, default=3, help="look-ahead depth of T* (default 3)")
    ap.add_argument("--time-limit", type=float, default=None, help="seconds; returns TIMEOUT when reached")
    ap.add_argument("--verbose", "-v", action="store_true")
    a = ap.parse_args(argv)

    if a.distances is not None:
        text = a.distances
    elif a.file:
        with open(a.file, encoding="utf-8") as f:
            text = f.read()
    else:
        text = sys.stdin.read()
    try:
        nums = parse_numbers(text)
        inst = Instance.from_points(nums) if a.points else Instance.from_distances(nums)
    except InvalidInput as e:
        print(f"invalid input: {e}", file=sys.stderr)
        return 2

    o = Options(all_solutions=a.all, tstar_k=a.k, time_limit=a.time_limit, verbose=a.verbose)
    if a.strict:
        o = o.strict()
    o.quick_search &= not a.no_quick
    o.pseudo_poly &= not a.no_pseudo
    o.lift &= not a.no_lift
    o.positivity &= not a.no_positivity
    o.vieta &= not a.no_vieta
    o.tstar &= not a.no_tstar
    o.fallback &= not a.no_fallback

    rep = solve(inst, o)
    print(rep.to_json() if a.json else rep.summary())
    return {"YES": 0, "NO": 1}.get(rep.status, 3)


if __name__ == "__main__":
    sys.exit(main())
