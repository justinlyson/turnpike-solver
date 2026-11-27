"""Distance multisets: parsing, validation, normalisation and verification.

A turnpike instance is the multiset D of the n(n-1)/2 pairwise distances of n distinct integer points.
Everything in the solver works on an `Instance`; every answer is checked with `verify` before it is returned.
"""
from __future__ import annotations

import math
import re
from collections import Counter
from decimal import Decimal, InvalidOperation
from typing import Iterable, Sequence

import numpy as np


class InvalidInput(ValueError):
    """The input cannot be a turnpike distance multiset (wrong count, non-integer, non-positive, ...)."""


def _as_int(x, what: str = "value") -> int:
    """x as a Python int, or InvalidInput if x is not exactly an integer (1.5, np.float32(2.5), Fraction(3, 2), True,
    '7' are all rejected; 2.0, np.int64(2), Decimal('2') are accepted)."""
    if isinstance(x, (bool, np.bool_, str, bytes)):
        raise InvalidInput(f"{what} {x!r} is not an integer")
    try:
        v = int(x)
        exact = bool(v == x)
    except (TypeError, ValueError, OverflowError, ArithmeticError):
        exact = False
    if not exact:
        raise InvalidInput(f"{what} {x!r} is not an integer")
    return v


def _all_ints(xs) -> bool:
    return all(type(x) is int or isinstance(x, np.integer) for x in xs)


def _int_list(xs, what: str) -> list:
    xs = list(xs)
    return [int(x) for x in xs] if _all_ints(xs) else [_as_int(x, what) for x in xs]


def distances_of(points: Sequence[int]) -> Counter:
    """Multiset of pairwise distances of `points` (exact integers; numpy fast path when values fit in int64)."""
    pts = sorted(_int_list(points, "point"))
    n = len(pts)
    if n < 2:
        return Counter()
    if pts[-1] - pts[0] < 2 ** 62 and n > 64:
        a = np.asarray(pts, dtype=np.int64)
        iu = np.triu_indices(n, 1)
        diff = (a[None, :] - a[:, None])[iu]
        vals, cnts = np.unique(diff, return_counts=True)
        return Counter(dict(zip(map(int, vals), map(int, cnts))))
    c: Counter = Counter()
    for j in range(1, n):
        pj = pts[j]
        for i in range(j):
            c[pj - pts[i]] += 1
    return c


def canonical(points: Iterable[int]) -> list[int]:
    """Shift to start at 0 and pick the lexicographically smaller of the set and its mirror image."""
    pts = sorted(int(p) for p in points)
    if not pts:
        return []
    pts = [p - pts[0] for p in pts]
    w = pts[-1]
    mir = sorted(w - p for p in pts)
    return min(pts, mir)


class Instance:
    """A validated distance multiset.

    counts : distance -> multiplicity (all distances positive integers)
    n      : number of points, from |D| = n(n-1)/2
    W      : the largest distance (the span; the two end points are 0 and W)
    """

    def __init__(self, counts: dict | None = None, n: int = 0, W: int = 0, arrays=None):
        self._counts = counts            # distance -> multiplicity (built lazily from arrays when needed)
        self._vc = arrays                # (distinct ascending, multiplicities) as int64 arrays, or None
        self.n = n
        self.W = W
        self._sorted = None
        self._arr = None

    def __repr__(self) -> str:
        return f"Instance(n={self.n}, W={self.W})"

    # ------------------------------------------------------------------ construction
    @classmethod
    def from_distances(cls, dists: Iterable[int] | dict) -> "Instance":
        if isinstance(dists, dict):
            c = Counter()
            for d, m in dists.items():                 # 2 and 2.0 are the same distance: add, don't overwrite
                c[_as_int(d, "distance")] += _as_int(m, "multiplicity")
            c = Counter({d: m for d, m in c.items() if m != 0})
        else:
            dl = list(dists)
            if not _all_ints(dl):
                dl = [_as_int(d, "distance") for d in dl]
            if len(dl) > 5000 and all(isinstance(d, (int, np.integer)) for d in dl) and max(dl) < 2 ** 62 and min(dl) > -(2 ** 62):
                v, k = np.unique(np.asarray(dl, dtype=np.int64), return_counts=True)
                return cls._from_arrays(v, k.astype(np.int64))
            c = Counter(int(d) for d in dl)
        total = sum(c.values())
        if any(m < 0 for m in c.values()):
            raise InvalidInput("negative multiplicity")
        if any(d <= 0 for d in c):
            raise InvalidInput("distances must be positive (points must be distinct)")
        n = (1 + math.isqrt(1 + 8 * total)) // 2
        if n * (n - 1) // 2 != total:
            raise InvalidInput(f"{total} distances is not n(n-1)/2 for any n")
        W = max(c) if c else 0
        if c and c[W] != 1:
            raise InvalidInput(f"the largest distance {W} occurs {c[W]} times; it must occur once")
        return cls(counts=dict(c), n=n, W=W)

    @classmethod
    def _from_arrays(cls, v: np.ndarray, k: np.ndarray) -> "Instance":
        """Validated construction from sorted distinct values v and multiplicities k (numpy fast path)."""
        total = int(k.sum())
        if v.size and v[0] <= 0:
            raise InvalidInput("distances must be positive (points must be distinct)")
        n = (1 + math.isqrt(1 + 8 * total)) // 2
        if n * (n - 1) // 2 != total:
            raise InvalidInput(f"{total} distances is not n(n-1)/2 for any n")
        if v.size and k[-1] != 1:
            raise InvalidInput(f"the largest distance {int(v[-1])} occurs {int(k[-1])} times; it must occur once")
        return cls(n=n, W=int(v[-1]) if v.size else 0, arrays=(v, k))

    @classmethod
    def from_points(cls, points: Iterable[int]) -> "Instance":
        pts = sorted(_int_list(points, "point"))
        if len(set(pts)) != len(pts):
            raise InvalidInput("points must be distinct")
        n = len(pts)
        if n > 64 and pts[-1] - pts[0] < 2 ** 62:
            a = np.asarray(pts, dtype=np.int64)
            iu = np.triu_indices(n, 1)
            v, k = np.unique((a[None, :] - a[:, None])[iu], return_counts=True)
            return cls._from_arrays(v, k.astype(np.int64))
        return cls.from_distances(distances_of(pts))

    # ------------------------------------------------------------------ views
    @property
    def counts(self) -> dict:
        """distance -> multiplicity."""
        if self._counts is None:
            v, k = self._vc
            self._counts = dict(zip(v.tolist(), k.tolist()))
        return self._counts

    @property
    def total(self) -> int:
        return self.n * (self.n - 1) // 2

    @property
    def distinct(self) -> list[int]:
        """Distinct distances, ascending."""
        if self._sorted is None:
            self._sorted = self._vc[0].tolist() if self._vc is not None else sorted(self.counts)
        return self._sorted

    def arrays(self):
        """(distinct distances ascending, multiplicities) as int64 arrays (W < 2^62); cached."""
        if self._vc is None:
            m = len(self.counts)
            v = np.fromiter(self.counts.keys(), dtype=np.int64, count=m)
            c = np.fromiter(self.counts.values(), dtype=np.int64, count=m)
            o = np.argsort(v, kind="stable")
            self._vc = (v[o], c[o])
        return self._vc

    def sorted_array(self) -> np.ndarray:
        """All distances with multiplicity, ascending, as int64 (W < 2^62)."""
        if self._arr is None:
            v, c = self.arrays()
            self._arr = np.repeat(v, c)
        return self._arr

    def gcd(self) -> int:
        if self._vc is not None:
            return int(np.gcd.reduce(self._vc[0])) or 1
        g = 0
        for d in self.counts:
            g = math.gcd(g, d)
        return g or 1

    def scaled_down(self, g: int) -> "Instance":
        if self._vc is not None and self._counts is None:
            return Instance(n=self.n, W=self.W // g, arrays=(self._vc[0] // g, self._vc[1]))
        return Instance(counts={d // g: m for d, m in self.counts.items()}, n=self.n, W=self.W // g)

    def as_list(self) -> list[int]:
        out = []
        for d in self.distinct:
            out.extend([d] * self.counts[d])
        return out


def verify(points: Sequence[int], inst: Instance) -> bool:
    """True iff `points` are n distinct integers whose distance multiset is exactly inst.counts."""
    if len(points) != inst.n or len(set(points)) != inst.n:
        return False
    pts = sorted(int(p) for p in points)
    if pts[-1] - pts[0] < 2 ** 62 and inst.n > 64:          # fast path: compare sorted int64 arrays
        a = np.asarray(pts, dtype=np.int64)
        iu = np.triu_indices(inst.n, 1)
        got = np.sort((a[None, :] - a[:, None])[iu])
        return np.array_equal(got, inst.sorted_array())
    return dict(distances_of(pts)) == inst.counts


_SPLIT = re.compile(r"[\s,;\[\]\(\)\{\}]+")


def parse_numbers(text: str) -> list[int]:
    """Read integers separated by whitespace, commas, semicolons or brackets."""
    out = []
    for tok in _SPLIT.split(text.replace(chr(0xFEFF), "").strip()):
        if not tok:
            continue
        try:
            v = int(tok)
        except ValueError:
            try:                                    # '12.0', '1e3': exact decimal arithmetic, never float
                f = Decimal(tok)
            except InvalidOperation:
                raise InvalidInput(f"'{tok}' is not a number") from None
            if not f.is_finite() or f != f.to_integral_value():
                raise InvalidInput(f"'{tok}' is not an integer") from None
            v = int(f)
        out.append(v)
    return out
