"""Solver stages. Every stage module exposes

    NAME      short identifier used in reports
    SECTION   section number in the README
    COST      'POLY', 'PSEUDO' or 'EXP' (see turnpike.options)
    run(inst, opts, ctx) -> StageResult

Contract (checked by the orchestrator):
  * status 'YES': `solutions` holds at least one point list; every list is re-verified against the input.
  * status 'NO' : the stage has PROVED that no realization exists (exact arithmetic only).
  * status 'UNDECIDED': the stage could not decide within its rules; the next stage runs.
  * A stage never returns a wrong YES (verification would catch it) and never a NO without proof.
"""
from __future__ import annotations

import time
from dataclasses import dataclass, field


@dataclass
class StageResult:
    status: str                                  # 'YES' | 'NO' | 'UNDECIDED' | 'TIMEOUT'
    solutions: list = field(default_factory=list)
    work: dict = field(default_factory=dict)     # work counts (nodes, factors, LPs, ...)
    note: str = ""


class Deadline:
    """Wall-clock limit shared by all stages. `expired()` is cheap to call often."""

    def __init__(self, seconds: float | None):
        self.t_end = None if seconds is None else time.perf_counter() + seconds

    def expired(self) -> bool:
        return self.t_end is not None and time.perf_counter() > self.t_end

    def remaining(self) -> float | None:
        return None if self.t_end is None else max(0.0, self.t_end - time.perf_counter())


class TimeUp(Exception):
    pass
