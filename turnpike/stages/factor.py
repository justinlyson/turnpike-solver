"""Section 4: direct factoring of the autocorrelation polynomial (PSEUDO: cost poly(W / gcd)).

With g = gcd(D) every realization lies in gZ, so divide by g and factor
    P(z) = z^W' (n + sum_d D(d) (z^{d/g} + z^{-d/g})),   W' = W / g,
over Z.  `turnpike.polyfactor` enumerates the 0/1 factors A with P = A A^*; their exponent sets are the
realizations (times g).  The enumeration is exhaustive, so an empty result proves NO.
"""
from __future__ import annotations

from ..distances import Instance
from ..options import Options
from ..polyfactor import enumerate_realizations, univariate_P
from . import Deadline, StageResult

NAME = "factor"
SECTION = 4
COST = "PSEUDO"

SHARE = 0.5        # fraction of the remaining time this stage may use, so later stages still get to run


def run(inst: Instance, opts: Options, ctx: dict) -> StageResult:
    if inst.n <= 2:
        return StageResult("YES", solutions=[[0, inst.W] if inst.n == 2 else [0]])
    g = inst.gcd()
    width = inst.W // g
    if not opts.pseudo_poly:
        return StageResult("UNDECIDED", note="skipped: pseudo-polynomial stages are off")
    if width > opts.max_factor_width:
        return StageResult("UNDECIDED", note=f"skipped: W/gcd = {width} > max_factor_width")
    deadline = ctx.get("deadline")
    small = inst.scaled_down(g) if g > 1 else inst
    P, _ = univariate_P(inst.n, small.counts.items())

    def decode(exps):
        return [e - exps[0] for e in exps]

    remaining = None if deadline is None else deadline.remaining()
    share = None if remaining is None else Deadline(remaining * SHARE)
    res = enumerate_realizations(P, inst.n, decode, small, opts.all_solutions, share)
    work = {"degree": P.degree(), **res.work}
    if res.timed_out:
        if deadline.expired():
            return StageResult("TIMEOUT", work=work, note=res.work.get("why", "deadline"))
        return StageResult("UNDECIDED", work=work, note=f"stopped at {SHARE:.0%} of the remaining time: "
                                                        + res.work.get("why", "enumeration"))
    sols = [[g * x for x in pts] for pts in res.solutions]
    if sols:
        return StageResult("YES", solutions=sols, work=work)
    return StageResult("NO", work=work, note="factorization: no 0/1 factor A with P = A A*")
