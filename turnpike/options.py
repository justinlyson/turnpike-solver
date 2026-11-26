"""Solver options.

Every stage can be switched off. Stages are marked by their worst-case cost class:

  POLY        polynomial in n and log W (unconditionally, or under the stated conjecture for T*)
  PSEUDO      polynomial in the numbers themselves (W or the lift height), not in log W.
              Very fast in practice for moderate widths; switch off with pseudo_poly=False to run the
              strictly polynomial pipeline only.
  EXP         exponential worst case (the exhaustive fallback). Never needed on any input we have tested;
              if it is ever used the report says so and asks for the input to be reported.
"""
from __future__ import annotations

from dataclasses import dataclass


@dataclass
class Options:
    # ---------------------------------------------------------------- what to return
    all_solutions: bool = False          # enumerate every realization (up to mirror image), not just one

    # ---------------------------------------------------------------- quick search (POLY: node budget n^quick_exp)
    quick_search: bool = True
    quick_exp: float = 2.0               # node budget = max(64, n ** quick_exp)

    # ---------------------------------------------------------------- positivity certificate (POLY, proves NO only)
    positivity: bool = True

    # ---------------------------------------------------------------- direct factoring (PSEUDO: cost poly(W))
    pseudo_poly: bool = True             # master switch for every PSEUDO trick (direct factoring, large lift boxes)
    max_factor_width: int = 3_000        # direct factoring only when W / gcd <= this (flint: ~4 s at degree 2000)

    # ---------------------------------------------------------------- lift + factor (POLY when the height is poly(n))
    lift: bool = True
    lift_max_height_exp: float = 3.0     # strict mode: lift height must be <= n ** this
    lift_max_box: int = 20_000_000       # PSEUDO mode: allow any lift whose projected degree is <= this

    # ---------------------------------------------------------------- anchors + square root (POLY)
    vieta: bool = True

    # ---------------------------------------------------------------- propagation + exact LP search, T*_k (POLY under (ST_k))
    tstar: bool = True
    tstar_k: int = 3

    # ---------------------------------------------------------------- exhaustive fallback (EXP)
    fallback: bool = True
    fallback_node_limit: int | None = None

    # ---------------------------------------------------------------- misc
    time_limit: float | None = None      # seconds; when hit the status is TIMEOUT (never a wrong answer)
    seed: int = 1
    verbose: bool = False

    def strict(self) -> "Options":
        """Copy with every non-polynomial speed trick disabled (PSEUDO off, fallback off)."""
        o = Options(**self.__dict__)
        o.pseudo_poly = False
        o.fallback = False
        return o
