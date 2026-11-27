"""The orchestrator: runs the sections in order until one decides, verifies every answer, and writes the report.

Order (each can be switched off in Options; cheap sections first):
  0  preprocess  divide all distances by their gcd; trivial n <= 2
  1  quick       classical backtracking with a polynomial node budget              POLY
  2  positivity  |A(e^it)|^2 >= 0: one certified negative angle proves NO           POLY
  3  vieta       mirror-pair anchors + exact sparse square root                     POLY
  4  factor      direct factoring of P(z) = A(z) A(1/z)                             PSEUDO (cost grows with W)
  5  lift        lift to the 3-term relation lattice, factor (T3U / T3R)            POLY when the lift height is poly(n)
  6  tstar       propagation + exact LP search, no backtracking allowed             POLY under the (ST_k) conjecture
  7  fallback    exhaustive backtracking                                            EXP (never needed on any tested input)

Soundness: a YES is returned only after the full distance multiset of the answer has been recomputed and compared with
the input; a NO only comes from a section that proves it (exhaustive search inside its budget, complete factor
enumeration, or an exact certificate).
"""
from __future__ import annotations

import importlib
import time
import traceback

from .distances import Instance, canonical, verify
from .options import Options
from .report import Report, StageAttempt
from .stages import Deadline, StageResult
from .stages import backtrack, vieta


# run order: (name, module or loader, cost class, option attribute that enables it)
_ORDER = [
    ("quick", lambda: backtrack.Quick, "POLY", "quick_search"),
    ("positivity", "positivity", "POLY", "positivity"),
    ("vieta", lambda: vieta, "POLY", "vieta"),
    ("factor", "factor", "PSEUDO", "pseudo_poly"),
    ("lift", "lift", "POLY", "lift"),
    ("tstar", "tstar", "POLY", "tstar"),
    ("fallback", lambda: backtrack.Fallback, "EXP", "fallback"),
]


def _stage(loader):
    """Import a stage module only when it is about to run (keeps start-up fast). None if absent from this build."""
    if callable(loader):
        return loader()
    try:
        return importlib.import_module(f".stages.{loader}", __package__)
    except ImportError:
        return None


def _is_complete(stage_name: str, res: StageResult) -> bool:
    """Does a YES from this stage list every realization (needed for all_solutions)?"""
    if stage_name in ("factor", "lift", "fallback"):
        return True
    return bool(res.work.get("complete"))


def solve(data, options: Options | None = None) -> Report:
    """Solve a turnpike instance.

    data: an Instance, a dict distance -> multiplicity, or an iterable of distances.
    Returns a Report (status YES / NO / UNDECIDED / TIMEOUT, points, which section decided, per-section work).
    """
    opts = options or Options()
    t0 = time.perf_counter()
    inst = data if isinstance(data, Instance) else Instance.from_distances(data)
    rep = Report(status="UNDECIDED", n=inst.n, W=inst.W)
    deadline = Deadline(opts.time_limit)

    # ---------------------------------------------------------------- 0. preprocess
    if inst.n <= 2:
        pts = [0] if inst.n == 1 else [0, inst.W]
        rep.status, rep.points, rep.solutions, rep.decided_by = "YES", pts, [pts], "preprocess"
        rep.seconds = time.perf_counter() - t0
        return rep
    g = inst.gcd()
    work_inst = inst.scaled_down(g) if g > 1 else inst
    rep.attempts.append(StageAttempt("preprocess", 0, "POLY", "DONE", 0.0, {"gcd": g}))

    found: dict = {}                                   # canonical tuple -> points (original scale)
    complete = False                                   # all_solutions: a section that lists every realization decided
    ctx = {"deadline": deadline}
    for sec, (name, loader, cost, flag) in enumerate(_ORDER, start=1):
        if not getattr(opts, flag):
            rep.attempts.append(StageAttempt(name, sec, cost, "SKIPPED", note=f"{flag} off"))
            continue
        if deadline.expired():
            rep.status = "TIMEOUT"
            break
        stage = _stage(loader)
        if stage is None:
            continue                                   # module not present in this build
        ts = time.perf_counter()
        try:
            res: StageResult = stage.run(work_inst, opts, ctx)
        except Exception as e:                         # a stage bug must never become a wrong answer
            rep.attempts.append(StageAttempt(name, sec, cost, "ERROR", time.perf_counter() - ts,
                                             note=f"{type(e).__name__}: {e}"))
            if opts.verbose:
                traceback.print_exc()
            continue
        att = StageAttempt(name, sec, cost, res.status, time.perf_counter() - ts, dict(res.work), res.note)
        rep.attempts.append(att)
        if name == "fallback":
            rep.fallback_used = True
        if res.work.get("st_violation"):
            rep.guarantee_breach = True

        if res.status == "YES":
            good = []
            for s in res.solutions:
                pts = [p * g for p in s]
                if verify(pts, inst):
                    good.append(canonical(pts))
            if not good:
                att.outcome = "ERROR"
                att.note = "returned points failed verification (ignored)"
                continue
            for s in good:
                found.setdefault(tuple(s), s)
            if not opts.all_solutions or _is_complete(name, res):
                rep.status, rep.decided_by = "YES", name
                complete = bool(opts.all_solutions)
                break
        elif res.status == "NO":
            if found:                                  # contradiction: a verified YES exists
                att.outcome = "ERROR"
                att.note = "claimed NO but a verified realization exists (ignored)"
                continue
            rep.status, rep.decided_by = "NO", name
            break
        elif res.status == "TIMEOUT":
            rep.status = "TIMEOUT"
            break

    if opts.all_solutions and found:
        rep.complete = complete
    if found and rep.status != "YES":
        # all_solutions requested but only partial sections answered: still a YES, flagged complete = False
        rep.status = "YES"
        rep.decided_by = rep.decided_by or next(a.name for a in rep.attempts if a.outcome == "YES")
    if found:
        rep.solutions = sorted(found.values())
        rep.points = rep.solutions[0]
    rep.seconds = time.perf_counter() - t0
    return rep
