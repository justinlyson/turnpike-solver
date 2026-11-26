"""Section 2: positivity certificate: a fast, exact proof of NO.

If X realizes D, then with A(z) = sum_{x in X} z^x,

    f(theta) = |A(e^{i theta})|^2 = n + 2 * sum_d D(d) cos(d theta)  >=  0   for every theta.

So one angle theta with f(theta) < 0 proves that no realization exists. (A random perturbation of a valid list almost
always has such an angle: on every no-solution mutant we generated, one was found.)

Search: candidate angles are theta = 2 pi p / q with integers p, q, so d*theta mod 2 pi = 2 pi ((d p) mod q) / q is
reduced EXACTLY in integer arithmetic, for any size of d. For a grid size N, folding the distances mod N and taking one
FFT gives f at ALL N angles 2 pi k / N at once (cost O(|D| + N log N)). Narrow lists: N >= 8W covers every angle that
matters. Wide lists: a few grids of a few hundred thousand angles (sizes 2^18, 3^11, 5^8: different angle sets), then a
small number of random angles. Certification: the sum is re-evaluated with ball arithmetic (python-flint arb),
which gives a rigorous enclosure; NO is returned only if the whole enclosure is below 0.

Cost: O(|D| log W) per angle; polynomial in n and log W (POLY). It can only prove NO; otherwise UNDECIDED.
"""
from __future__ import annotations

import numpy as np

from . import StageResult

NAME = "positivity"
SECTION = 2
COST = "POLY"

_FFT_MAX = 1 << 21          # FFT size limit for the exhaustive grid
_WIDE_GRIDS = (1 << 18, 3 ** 11, 5 ** 8)   # aliased FFT grids for wide lists
_Q_WIDE = 2147483647        # prime 2^31 - 1: products (d mod q) * p fit in int64
_WORK = 3 * 10 ** 6         # random-angle budget, in angle x distance evaluations


def certify(counts: dict, n: int, p: int, q: int, prec: int = 128) -> bool:
    """Rigorous check that f(2 pi p / q) < 0 (ball arithmetic)."""
    import flint
    old = flint.ctx.prec
    flint.ctx.prec = prec
    try:
        s = flint.arb(n)
        for d, m in counts.items():
            r = (d * p) % q
            s += 2 * m * flint.arb.cos_pi(flint.arb(flint.fmpq(2 * r, q)))
        return bool(s < 0)
    finally:
        flint.ctx.prec = old


def _grid_candidates(counts: dict, n: int, N: int, k: int = 5):
    """f at all angles 2 pi j / N (exact folding of d mod N, one FFT); up to k most negative (j, N, value)."""
    a = np.zeros(N)
    ds = np.fromiter((d % N for d in counts), dtype=np.int64, count=len(counts))
    np.add.at(a, ds, np.fromiter(counts.values(), dtype=np.float64, count=len(counts)))
    f = n + 2.0 * np.real(np.fft.fft(a))
    idx = np.argpartition(f, min(k, N - 1))[:k]
    idx = idx[np.argsort(f[idx])]
    return [(int(i), N, float(f[i])) for i in idx if f[i] < 0]


def _random_candidates(counts: dict, n: int, seed: int, k: int = 5):
    rng = np.random.default_rng(seed)
    q = _Q_WIDE
    ds = np.array([d % q for d in counts], dtype=np.int64)
    ms = np.array(list(counts.values()), dtype=np.float64)
    samples = int(min(1 << 15, max(256, _WORK // max(1, len(ds)))))
    best = []
    for start in range(0, samples, 256):
        ps = rng.integers(1, q, size=256, dtype=np.int64)
        r = (ds[None, :] * ps[:, None]) % q                      # exact: both factors < 2^31
        f = n + 2.0 * (np.cos(2 * np.pi * r / q) @ ms)
        for j in np.argsort(f)[:k]:
            if f[j] < 0:
                best.append((int(ps[j]), q, float(f[j])))
        if best:
            break
    best.sort(key=lambda t: t[2])
    return best[:k]


def run(inst, opts, ctx) -> StageResult:
    counts, n, W = inst.counts, inst.n, inst.W
    tried = []
    if 8 * W <= _FFT_MAX:
        grids = [1 << int(np.ceil(np.log2(max(64, 8 * W))))]
    else:
        grids = list(_WIDE_GRIDS)
    for N in grids:
        tried.append(N)
        for p, q, val in _grid_candidates(counts, n, N):
            if certify(counts, n, p, q):
                return StageResult("NO", work={"grids": tried, "theta": f"2pi*{p}/{q}", "f": round(val, 6)},
                                   note="f(theta) = |A|^2 < 0 at a certified angle")
    if 8 * W > _FFT_MAX:
        for p, q, val in _random_candidates(counts, n, opts.seed):
            if certify(counts, n, p, q):
                return StageResult("NO", work={"grids": tried, "theta": f"2pi*{p}/{q}", "f": round(val, 6)},
                                   note="f(theta) = |A|^2 < 0 at a certified angle")
    return StageResult("UNDECIDED", work={"grids": tried}, note="no certified negative angle found")
