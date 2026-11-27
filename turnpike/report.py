"""The solver's report: the answer, which section decided it, and what every section did."""
from __future__ import annotations

import json
from dataclasses import asdict, dataclass, field

ISSUE_URL = "https://github.com/justinlyson/turnpike-solver/issues"


@dataclass
class StageAttempt:
    name: str
    section: int
    cost: str                # POLY | PSEUDO | EXP
    outcome: str             # YES | NO | UNDECIDED | SKIPPED | TIMEOUT | ERROR
    seconds: float = 0.0
    work: dict = field(default_factory=dict)
    note: str = ""


@dataclass
class Report:
    status: str                                   # YES | NO | UNDECIDED | TIMEOUT
    points: list | None = None                    # one realization (canonical: starts at 0, lexicographically <= mirror)
    solutions: list = field(default_factory=list) # all realizations found (up to mirror), when requested
    complete: bool | None = None                  # all_solutions only: True = `solutions` is every realization;
                                                  # False = only partial sections answered, more may exist
    decided_by: str | None = None
    attempts: list = field(default_factory=list)
    fallback_used: bool = False                   # the exponential section ran
    guarantee_breach: bool = False                # T* met a first-wrong child it could not refute (an (ST_k) violation)
    seconds: float = 0.0
    n: int = 0
    W: int = 0

    @property
    def needs_report(self) -> bool:
        return self.fallback_used or self.guarantee_breach

    def to_dict(self) -> dict:
        d = asdict(self)
        d["needs_report"] = self.needs_report
        return d

    def to_json(self) -> str:
        return json.dumps(self.to_dict())

    def summary(self) -> str:
        lines = [f"status      : {self.status}   (n = {self.n}, W = {self.W}, {self.seconds:.3f} s)"]
        if self.points is not None:
            pts = self.points if len(self.points) <= 40 else self.points[:20] + ["..."] + self.points[-5:]
            lines.append(f"points      : {pts}")
        if self.complete is False:
            lines.append(f"realizations: {len(self.solutions)} found (up to mirror image); the list MAY BE INCOMPLETE:"
                         " no section that lists every realization finished")
        elif len(self.solutions) > 1:
            lines.append(f"realizations: {len(self.solutions)} (up to mirror image)")
        lines.append(f"decided by  : {self.decided_by}")
        lines.append("sections    :")
        for a in self.attempts:
            w = ", ".join(f"{k}={v}" for k, v in a.work.items())
            lines.append(f"  [{a.section}] {a.name:<12} {a.cost:<6} {a.outcome:<9} {a.seconds:8.3f} s  {w}"
                         + (f"  ({a.note})" if a.note else ""))
        if self.needs_report:
            lines.append("")
            lines.append("*** This input needed more than the polynomial-time sections. ***")
            if self.guarantee_breach:
                lines.append("    T* met a wrong branch it could not refute within its k levels.")
            if self.fallback_used:
                lines.append("    The exhaustive (exponential worst-case) fallback was used.")
            lines.append("    We have never seen this happen. Please report the input at")
            lines.append(f"    {ISSUE_URL}")
        return "\n".join(lines)
