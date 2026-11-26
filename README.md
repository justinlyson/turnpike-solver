# turnpike-solver

Reconstruct integer points on a line from the multiset of their pairwise distances (the **turnpike problem**,
also called the partial digest problem).

Give it the n(n−1)/2 distances, in any order. It returns a set of points with exactly those distances, or a proof
that none exists. It also reports which section of the algorithm decided the answer.

```bash
pip install -e .
python -m turnpike -d "1,4,10,12,17,3,9,11,16,6,8,13,2,7,5" --all
```

```
status      : YES   (n = 6, W = 17, 0.004 s)
points      : [0, 1, 4, 10, 12, 17]
realizations: 2 (up to mirror image)
decided by  : quick
```

```python
from turnpike import solve, Options
rep = solve([1, 4, 10, 12, 17, 3, 9, 11, 16, 6, 8, 13, 2, 7, 5])
rep.status, rep.points, rep.decided_by
```

A short paper describing the approach, what is proved and the measurements: [paper/turnpike_solver.pdf](paper/turnpike_solver.pdf).

## Guarantees

- **Every YES is verified.** The full distance multiset of the answer is recomputed and compared with the input before
  it is returned.
- **Every NO is proved.** It comes from an exhaustive search inside its budget, a complete enumeration of polynomial
  factorizations, or an exact certificate. Floating-point results are never trusted for a NO.
- **The worst case is not proved polynomial.** Sections 1–6 run in polynomial time (section 4 is pseudo-polynomial and
  can be switched off). Sections 1–5 are partial: each decides some inputs and hands the rest on. Section 6 decides every input in polynomial time *if* one open conjecture
  holds (see below). Section 7 is an exhaustive fallback that is exponential in the worst case. **On every input we have
  tested, including every hard family we know of and every one we built to break the method, section 7 has never been
  needed.** If it ever is, the solver says so and asks you to report the input.

## The sections

The solver runs the sections in this order and stops at the first one that decides.

| # | Section | What it does | Cost class | Flag to disable |
|---|---------|--------------|------------|-----------------|
| 0 | preprocess | Divide all distances by their gcd; handle n ≤ 2. | polynomial | — |
| 1 | quick | Classical backtracking (Skiena–Smith–Lemke): the largest unexplained distance must touch an end point. Node budget n². | polynomial (budgeted) | `--no-quick` |
| 2 | positivity | If X is a realization, f(θ) = \|A(e^{iθ})\|² = n + 2 Σ D(d) cos(dθ) is never negative. One angle with f(θ) < 0, certified with rigorous ball arithmetic, proves NO. In our tests it catches nearly every no-solution input in milliseconds. | polynomial | `--no-positivity` |
| 3 | vieta | Finds mirror pairs {p, W−p} ("anchors"), folds the rest into pairs, and recovers which side each point is on as a root of t² − S t + P, using an exact sparse square root. | polynomial | `--no-vieta` |
| 4 | factor | Factors P(z) = A(z)·A(1/z) directly and enumerates the 0/1 factor pairs. Complete. | **pseudo-polynomial** (grows with the width W) | `--no-pseudo` |
| 5 | lift | Lifts the distances to the lattice of their additive relations a + b = c, then factors the lifted polynomial. Complete. By default it factors directly in all r lattice coordinates (fastest in practice). In `--strict` mode it first projects to 3 random coordinates, the version with a proof of polynomial cost. | polynomial when the lift height is polynomial (projected version) | `--no-lift` |
| 6 | tstar | Depth-first search where every branch is tested by propagation plus an exact linear program, looking k levels ahead. Never backtracks: if a wrong branch survives the test, it stops and reports. | polynomial if conjecture (ST_k) holds | `--no-tstar`, `-k` |
| 7 | fallback | Exhaustive backtracking. | **exponential** worst case | `--no-fallback` |

`--strict` turns off every non-polynomial trick: section 4, direct factoring in many lift coordinates, large lifts,
and section 7. A strict run can return UNDECIDED, but never a wrong answer.

### Why these sections

- **Narrow inputs** (small W): section 4 factors directly. This includes the classic hard families for backtracking,
  such as Zhang's.
- **Inputs with a short lift**: section 5. Many wide inputs with additive structure are secretly low-dimensional. For
  example, gaps of "100000 plus small noise" decompose into (number of gaps spanned, total noise), and the lift sees that.
- **Rigid wide inputs** (long chains of relations such as 1, 2, 4, 8, …): their lift is tall, but the chains act as
  anchors. Section 1 (one-sided chains) or section 3 (mirror-symmetric chains) handles them.
- **No-solution inputs**: section 2 almost always proves NO at once. Sections 4–6 handle the rest.
- **Everything else**: section 6.

### The open conjecture

Section 6 is a search tree in which each branch is checked by a test (propagation plus an exact LP over k levels).
Conjecture (ST_k) says that for a fixed k, every wrong branch is rejected by the test. If it holds, the search makes
at most 1 + 2(n−1) branch tests, each polynomial, so the whole algorithm is polynomial. Every other step of that
argument has a written proof that has been independently reviewed. (ST_k) itself is the one open statement. With
k = 3 we know of no counterexample.

Implementation note: section 6's LPs are solved with HiGHS, and every answer it uses for a refutation is confirmed
exactly. Some of its LP calls have time limits. A limit can only turn an answer into "unknown" (less pruning), never
into a wrong NO.

## Options

```
python -m turnpike [FILE] [-d "..."] [--points] [--all] [--json] [--strict]
                   [--no-quick] [--no-positivity] [--no-vieta] [--no-pseudo] [--no-lift] [--no-tstar] [--no-fallback]
                   [-k K] [--time-limit SECONDS]
```

- `--points` reads the numbers as points and computes their distances first (handy for testing).
- `--all` lists every realization, up to mirror image. Different realizations with the same distances
  (homometric sets) exist.
- `--time-limit` returns TIMEOUT when the time is up, never a wrong answer.
- Exit codes: 0 = YES, 1 = NO, 3 = UNDECIDED or TIMEOUT, 2 = invalid input.

## Test results

- **Corpus** (`tests/data`, 97 cases): random sets up to n = 300 and W = 10⁹; homometric pairs and their products;
  Zhang's family, which is exponential for plain backtracking; hard narrow cores (HM(K,3,B^r)); rigid chains (doubling
  up to 2⁶⁴, hyperbolic); mirror-symmetric chain attacks up to W ≈ 4·10¹²; mirrored-chain products; 37 no-solution
  inputs. Every case is decided correctly, in 1.8 s total. In `--strict` mode every case is also decided, using
  polynomial sections only (quick 70, lift 17, positivity 10). **The fallback was never used.**
- **Exhaustive cross-check:** on 12,000 random small inputs (solvable and perturbed), in default and strict mode, the
  status and the complete solution set match exhaustive search exactly; likewise on 4,000 more (including narrow and
  wide products) in four modes: default, strict, all solutions, and algebra only.
- **Speed against plain backtracking:** Zhang nz = 13 takes 0.15 s against 58 s (708,664 nodes). Zhang no-solution
  variants take 0.1–0.5 s, where plain backtracking does not finish in 60 s.
- **Hard and large cases** (`bench/tricky.py`): products of hard cores (n = 420 and n = 1040; on the n = 420 ones plain
  backtracking is still unfinished after 10⁷–10⁸ nodes) take 1.7–2.3 s; Zhang nz = 25 at W = 10⁹ takes 0.9 s;
  mirror-symmetric chains (n = 172) take 1.3 s; random n = 2000 takes 0.4 s and n = 5000 takes 3.2 s. Each
  no-solution variant (one distance moved by 1) is decided in 0.1–1.1 s.

```bash
python -m pytest tests          # ~240 tests, under a minute (TURNPIKE_SLOW=1 adds the slow cases)
python bench/bench.py           # benchmark table, written to bench/results.md
```

## If the fallback runs

The report ends with:

```
*** This input needed more than the polynomial-time sections. ***
```

We have never seen this happen. Please open an issue with the input (the distances, or the points if you know them).
Such an input would be a counterexample to the conjecture behind section 6, or a case where the earlier sections need
improving. Either way we want to see it.

## Layout

```
turnpike/
  distances.py      instance parsing, validation, verification
  options.py        all switches, with their cost classes
  solver.py         the orchestrator (section order, verification, report)
  report.py         the report and its text / JSON output
  polyfactor.py     shared factoring core: P = A A*, enumeration of 0/1 factor pairs
  stages/
    backtrack.py    sections 1 and 7 (compiled kernel _bt_kernel.py via numba when available)
    positivity.py   section 2
    vieta.py        section 3
    factor.py       section 4
    lift.py         section 5 (compiled relation sweep _lift_kernel.py via numba when available)
    tstar.py        section 6 (helpers in _tstar/)
tests/              test corpus (hard families, NO inputs, homometric pairs) and tests
bench/              benchmark against plain backtracking (bench.py) and the hardest/largest cases (tricky.py)
paper/              short paper (LaTeX source and PDF)
```

## Requirements

Python ≥ 3.10, numpy, scipy, python-flint (fast exact factoring), highspy (LP). numba is optional: it gives a faster
backtracking kernel.

## Background

- S. Skiena, W. Smith, P. Lemke, *Reconstructing sets from interpoint distances*, SoCG 1990: the backtracking
  algorithm and the factoring approach.
- The turnpike problem is not known to be solvable in polynomial time in the worst case. This package is a research
  implementation of an approach under development.

## License and citation

Apache License 2.0 (see `LICENSE`). Anyone who redistributes this code or builds on it must keep the copyright and
`NOTICE` file, which credit the author. If you use it in research, please cite it using `CITATION.cff`
(GitHub shows a "Cite this repository" button).
