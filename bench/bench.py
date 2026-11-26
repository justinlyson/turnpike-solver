"""Benchmark: run the test corpus through turnpike.solve and through the plain backtracking baseline.

    python bench/bench.py [--family NAME] [--case SUBSTR] [--strict] [--timeout S] [--baseline-cap NODES]
                          [--no-baseline] [--inproc] [--out bench/results.md]

Each solve runs in its own child process (hard-killed by PID at timeout + 30 s) unless --inproc is given.
Columns: family, case, n, W, expected, status, check (ok / open / WRONG), decided_by, seconds, per-section seconds,
fallback_used, baseline status / nodes / seconds, speed-up (baseline seconds / solver seconds, when both decided).
"""
from __future__ import annotations

import argparse
import json
import os
import subprocess
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
for p in (ROOT, os.path.join(ROOT, "tests"), HERE):
    if p not in sys.path:
        sys.path.insert(0, p)

import families  # noqa: E402
from baseline import backtrack  # noqa: E402
from turnpike.distances import distances_of  # noqa: E402


def solve_one(case, strict: bool, timeout: float) -> dict:
    import turnpike
    from turnpike.options import Options
    opts = Options(time_limit=timeout)
    if strict:
        opts = opts.strict()
    t = time.perf_counter()
    try:
        rep = turnpike.solve(case.distances, opts)
    except Exception as e:  # report, do not crash the benchmark
        return {"status": "ERROR", "error": f"{type(e).__name__}: {e}", "seconds": time.perf_counter() - t}
    sec = time.perf_counter() - t
    d = rep.to_dict()
    ok_points = None
    if rep.points is not None:
        c = distances_of(rep.points)
        want = {}
        for x in case.distances:
            want[x] = want.get(x, 0) + 1
        ok_points = dict(c) == want
    sections = []
    for a in d.get("attempts", []):
        if a.get("outcome") != "SKIPPED":
            sections.append(f"{a.get('section')}.{a.get('name')}:{a.get('seconds', 0):.2f}")
    return {"status": rep.status, "decided_by": rep.decided_by, "seconds": sec, "sections": " ".join(sections),
            "fallback_used": rep.fallback_used, "guarantee_breach": rep.guarantee_breach, "points_ok": ok_points}


def check(case, r) -> str:
    st = r.get("status")
    if st in ("UNDECIDED", "TIMEOUT"):
        return "open"
    if st == "ERROR":
        return "ERROR"
    if st != case.expected:
        return "WRONG"
    if st == "YES" and not r.get("points_ok"):
        return "WRONG"
    return "ok"


def run_child(fam, name, strict, timeout):
    cmd = [sys.executable, os.path.abspath(__file__), "--one", fam, name, "--timeout", str(timeout)]
    if strict:
        cmd.append("--strict")
    try:
        p = subprocess.run(cmd, capture_output=True, text=True, timeout=timeout + 30)
    except subprocess.TimeoutExpired:   # subprocess.run kills the child by its own handle
        return {"status": "TIMEOUT", "decided_by": "hard-kill", "seconds": timeout + 30, "sections": "",
                "fallback_used": None}
    for line in reversed(p.stdout.splitlines()):
        if line.startswith("{"):
            return json.loads(line)
    return {"status": "ERROR", "error": (p.stderr or "")[-300:], "seconds": 0, "sections": ""}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--family")
    ap.add_argument("--case", help="substring filter on case names")
    ap.add_argument("--strict", action="store_true")
    ap.add_argument("--timeout", type=float, default=60.0)
    ap.add_argument("--baseline-cap", type=int, default=2_000_000)
    ap.add_argument("--no-baseline", action="store_true")
    ap.add_argument("--inproc", action="store_true")
    ap.add_argument("--out", default=os.path.join(HERE, "results.md"))
    ap.add_argument("--one", nargs=2, metavar=("FAMILY", "NAME"), help=argparse.SUPPRESS)
    a = ap.parse_args()

    if a.one:
        case = next(c for c in families.load_corpus(a.one[0]) if c.name == a.one[1])
        print(json.dumps(solve_one(case, a.strict, a.timeout)), flush=True)
        return

    fams = [a.family] if a.family else list(families.FAMILIES)
    rows = []
    hdr = ("| family | case | n | W | exp | status | check | decided_by | s | sections (s) | fallback "
           "| baseline | bl nodes | bl s | speed-up |")
    print(hdr)
    for fam in fams:
        for case in families.load_corpus(fam):
            if a.case and a.case not in case.name:
                continue
            n = families.case_n(case); W = case.distances[-1]
            r = solve_one(case, a.strict, a.timeout) if a.inproc else run_child(fam, case.name, a.strict, a.timeout)
            chk = check(case, r)
            if a.no_baseline:
                bl = ("-", 0, 0.0)
            else:
                t = time.perf_counter()
                st, pts, nodes, _ = backtrack(case.distances, node_cap=a.baseline_cap, time_cap=a.timeout)
                bl = (st, nodes, time.perf_counter() - t)
            sp = ""
            if chk == "ok" and bl[0] in ("YES", "NO") and r["seconds"] > 0:
                sp = f"{bl[2] / r['seconds']:.1f}x"
            Wtxt = str(W) if W < 10 ** 7 else f"{W:.2e}"
            row = (f"| {fam} | {case.name} | {n} | {Wtxt} | {case.expected} | {r.get('status')} | "
                   f"{'**' + chk + '**' if chk not in ('ok', 'open') else chk} | {r.get('decided_by')} | "
                   f"{r.get('seconds', 0):.3f} | {r.get('sections', '')} | {r.get('fallback_used')} | "
                   f"{bl[0]} | {bl[1]} | {bl[2]:.3f} | {sp} |")
            if r.get("error"):
                row += f" error: {r['error'][:120]!r}"
            print(row, flush=True)
            rows.append((row, chk, case))

    wrong = [c for _, k, c in rows if k == "WRONG"]
    with open(a.out, "w") as fh:
        fh.write(f"# Benchmark results\n\n{time.strftime('%Y-%m-%d %H:%M')}; strict={a.strict}; "
                 f"timeout={a.timeout}s; baseline node cap={a.baseline_cap}\n\n")
        fh.write(hdr + "\n|" + "---|" * 15 + "\n")
        for row, _, _ in rows:
            fh.write(row + "\n")
        fh.write(f"\n{len(rows)} cases; ok {sum(k == 'ok' for _, k, _ in rows)}, open "
                 f"{sum(k == 'open' for _, k, _ in rows)}, WRONG {len(wrong)}, ERROR "
                 f"{sum(k == 'ERROR' for _, k, _ in rows)}\n")
        for c in wrong:
            fh.write(f"\nWRONG ANSWER on {c.name}: distances = {c.distances}\n")
    print(f"\nwrote {a.out}; WRONG {len(wrong)}")


if __name__ == "__main__":
    main()
