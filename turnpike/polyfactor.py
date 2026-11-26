"""Realizations from the factorization of the autocorrelation polynomial.

For a point set X let A(z) = sum_{x in X} z^x.  Then A(z) A(1/z) = n + sum_d D(d) (z^d + z^-d) =: P(z), so
after the shift by z^W the polynomial P is palindromic of degree 2W and every realization of D is a 0/1 factor
A of P with P = A * A^*, where A^*(z) = z^{deg A} A(1/z) is the reciprocal.  The same holds verbatim for Laurent
polynomials in several variables (the lifted polynomial of `turnpike.stages.lift`), with A^* reversing every
variable.

Enumeration.  Factor P into irreducibles over Z.  Reciprocation permutes the irreducible factors:
  * a self-reciprocal factor g (g^* = +-g) occurs with even multiplicity 2m and A contains exactly g^m;
  * a reciprocal pair (f, f^*), f^* != +-f, occurs as f^e (f^*)^e and A contains f^j (f^*)^(e-j), 0 <= j <= e.
So the candidates are indexed by the tuples (j_1, ..., j_s).  Replacing every j_i by e_i - j_i turns A into
+-A^*, the mirror image of the same point set, so only tuples with j <= e - j (lexicographically) are tried.
A candidate is a realization candidate iff all its coefficients lie in {0, 1} up to a global sign.  Because
||A||^2 = P[middle] = n and A(1)^2 = P(1) = n^2, any factor with all coefficients in {-1, 0, 1} automatically has a
single sign and n terms, so the test reduces to "max |coefficient| <= 1".

Soundness.  Every realization is reached by some tuple, so a complete enumeration with no verified candidate
PROVES that D has no realization (exact integer arithmetic throughout).  Candidates are decoded to point sets by a
caller-supplied map and re-verified against the instance before being accepted.

Speed: towers.  Shortcut for large univariate P when a single realization is wanted.  A factor of the form F(z^s)
divides every polyphase component P_r(w) of P(z) = sum_r z^r P_r(z^s), so H(w) = gcd_r P_r(w) collects all of them
for that s, using gcds of degree deg P / s (a cheap gcd mod a prime screens each s first; the largest s with H != 1
is taken).  Product sets A(z) = A1(z) A2(z^s) give such factors.  Write P = R(z) H(z^s).  When deg R < 2s, a product
A1(z) G(z^s) with deg A1 < s has no colliding terms, so it is 0/1 iff A1 and G are; the 0/1 solutions of
R = A1 A1^* and H = G G^* are found separately (recursively, small degrees) and combined.  This route can only
produce a (verified) YES; when it finds nothing, the full factorization runs as usual.

Univariate leaves are first screened by their lowest and highest `LOW` coefficients (truncated products,
`fmpz_poly.mul_low`), which must already be 0/1; the full product is formed only for leaves passing the screen.
"""
from __future__ import annotations

import multiprocessing as mp
from dataclasses import dataclass, field
from typing import Callable, Iterable, Optional

import flint

from .distances import Instance, canonical, verify

LOW = 48                 # number of low/high coefficients used to screen univariate candidates
INLINE_DEGREE = 4000     # univariate P up to this degree is factored in-process (fast, not interruptible)
INLINE_TERMS = 20000     # multivariate P up to this many terms likewise
COARSE_DEGREE = 1500     # univariate P above this degree tries the F(z^s) block shortcut first
_SCREEN_P = 2147483647   # prime for the gcd screen


@dataclass
class FactorOutcome:
    """Result of an enumeration.

    complete  : the enumeration finished: every tuple was examined, or it stopped at the first solution because
                all solutions were not requested.  complete with `solutions` empty PROVES that D has no realization.
    solutions : distinct realizations, canonical form, as point lists of the instance passed in
    timed_out : the deadline expired during the enumeration
    """

    complete: bool = False
    solutions: list = field(default_factory=list)
    work: dict = field(default_factory=dict)
    timed_out: bool = False


# ----------------------------------------------------------------------------------------------------------------
# building P
def univariate_P(n: int, items: Iterable[tuple[int, int]]) -> tuple[flint.fmpz_poly, int]:
    """P(z) * z^m = n z^m + sum mult (z^{m+e} + z^{m-e}) for (e, mult) in items (e != 0); returns (poly, m)."""
    items = [(abs(e), c) for e, c in items]
    m = max(e for e, _ in items)
    co = [0] * (2 * m + 1)
    co[m] = n
    for e, c in items:
        co[m + e] += c
        co[m - e] += c
    return flint.fmpz_poly(co), m


def multivariate_P(n: int, items: Iterable[tuple[tuple, int]], ctx: flint.fmpz_mpoly_ctx) -> flint.fmpz_mpoly:
    """Multivariate analogue: n y^off + sum mult (y^{off+v} + y^{off-v}), off_i = max |v_i| (no variable divides P)."""
    items = list(items)
    k = ctx.nvars()
    off = [max(abs(v[i]) for v, _ in items) for i in range(k)]
    terms: dict = {tuple(off): n}
    for v, c in items:
        p = tuple(o + x for o, x in zip(off, v))
        q = tuple(o - x for o, x in zip(off, v))
        terms[p] = terms.get(p, 0) + c
        terms[q] = terms.get(q, 0) + c
    return ctx.from_dict(terms)


# ----------------------------------------------------------------------------------------------------------------
# reciprocals
def _rev_uni(f: flint.fmpz_poly) -> flint.fmpz_poly:
    return flint.fmpz_poly(f.coeffs()[::-1])


def _rev_multi(f: flint.fmpz_mpoly) -> flint.fmpz_mpoly:
    deg = f.degrees()
    return f.context().from_dict({tuple(d - e for d, e in zip(deg, m)): c for m, c in f.to_dict().items()})


def _is_one_bit(f) -> bool:
    """All coefficients in {-1, 0, 1}."""
    if isinstance(f, flint.fmpz_poly):
        return f.height_bits() <= 1
    return all(abs(int(c)) == 1 for c in f.coeffs())


def _exponents(f) -> list:
    if isinstance(f, flint.fmpz_poly):
        return [i for i, c in enumerate(f.coeffs()) if c != 0]
    return sorted(f.monoms())


# ----------------------------------------------------------------------------------------------------------------
# factoring under a deadline
def _encode(P):
    if isinstance(P, flint.fmpz_poly):
        return ("u", [int(c) for c in P.coeffs()], None)
    return ("m", {m: int(c) for m, c in P.to_dict().items()}, P.context().names())


def _decode_poly(kind, data, names):
    if kind == "u":
        return flint.fmpz_poly(data)
    return flint.fmpz_mpoly_ctx.get(tuple(names), "lex").from_dict(data)


def _factor_child(conn, kind, data, names):
    """Child-process entry: factor and send back (content, [(encoded factor, multiplicity)])."""
    c, facs = _decode_poly(kind, data, names).factor()
    conn.send((int(c), [(_encode(f)[1], int(e)) for f, e in facs]))
    conn.close()


def factor(P, deadline=None):
    """(content, [(f, e)]) as flint's factor(), or None if the deadline expires first.

    A flint factorization cannot be interrupted, so large ones run in a child process that is killed (only
    that child, by its own handle) when the remaining time runs out."""
    small = P.degree() <= INLINE_DEGREE if isinstance(P, flint.fmpz_poly) else len(P) <= INLINE_TERMS
    remaining = None if deadline is None else deadline.remaining()
    if small or remaining is None:
        c, facs = P.factor()
        return int(c), [(f, int(e)) for f, e in facs]
    kind, data, names = _encode(P)
    recv, send = mp.get_context("spawn").Pipe(duplex=False)
    proc = mp.get_context("spawn").Process(target=_factor_child, args=(send, kind, data, names), daemon=True)
    proc.start()
    send.close()
    try:
        if not recv.poll(remaining):
            return None
        c, facs = recv.recv()
    except EOFError:
        raise RuntimeError("factoring child process died") from None
    finally:
        if proc.is_alive():
            proc.kill()
        proc.join()
        recv.close()
    return c, [(_decode_poly(kind, f, names), e) for f, e in facs]


# ----------------------------------------------------------------------------------------------------------------
def _spread(g: flint.fmpz_poly, s: int) -> flint.fmpz_poly:
    """g(z^s)."""
    co = g.coeffs()
    out = [0] * (s * (len(co) - 1) + 1)
    out[::s] = co
    return flint.fmpz_poly(out)


def tower(P: flint.fmpz_poly, deadline=None):
    """(H, s, R) with P = R(z) H(z^s), s the largest with H = gcd_r P_r nonconstant, provided deg R < 2s; else None."""
    co = [int(c) for c in P.coeffs()]
    d = len(co) - 1
    for s in range(d // 2, 1, -1):
        if deadline is not None and (s & 255) == 0 and deadline.expired():
            return None
        comps = []
        for r in range(s):
            if any(co[r::s]):
                comps.append(r)
                if len(comps) == 2:
                    break
        if len(comps) < 2:
            continue
        g = flint.nmod_poly([x % _SCREEN_P for x in co[comps[0]::s]], _SCREEN_P)
        if g.gcd(flint.nmod_poly([x % _SCREEN_P for x in co[comps[1]::s]], _SCREEN_P)).degree() < 1:
            continue
        h = None
        for r in range(s):                               # exact gcd over Z of all nonzero components
            c = co[r::s]
            if any(c):
                h = flint.fmpz_poly(c) if h is None else h.gcd(flint.fmpz_poly(c))
                if h.degree() < 1:
                    break
        if h is None or h.degree() < 1:
            continue
        if h.coeffs()[-1] < 0:
            h = -h
        R, rem = divmod(P, _spread(h, s))
        if rem != 0 or R.degree() >= 2 * s:
            return None
        return h, s, R
    return None


def zero_one_factors(Q: flint.fmpz_poly, deadline=None, depth: int = 0):
    """0/1 polynomials A with A A^* = Q, as ascending exponent lists, one per mirror pair (complete when found by
    full factorization; possibly partial through towers).  None on deadline."""
    if Q.degree() <= 0:
        return [[0]]
    if Q.degree() > COARSE_DEGREE and depth < 6:
        t = tower(Q, deadline)
        if t is not None:
            h, s, R = t
            Gs = zero_one_factors(h, deadline, depth + 1)
            Rs = zero_one_factors(R, deadline, depth + 1)
            if Gs is None or Rs is None:
                return None
            out, seen = [], set()
            for r in Rs:
                m = max(r)
                for rr in {tuple(r), tuple(sorted(m - x for x in r))}:
                    for g in Gs:
                        A = sorted(x + s * y for x in rr for y in g)
                        c = tuple(canonical(A))
                        if c not in seen:
                            seen.add(c)
                            out.append(A)
            return out
    o = _enumerate(Q, factor(Q, deadline), None, None, None, True, deadline)
    if o.timed_out:
        return None
    return [list(c) for c in o.solutions]


# ----------------------------------------------------------------------------------------------------------------
def enumerate_realizations(P, n: int, decode: Callable[[list], Optional[list]], inst: Instance,
                           all_solutions: bool, deadline=None) -> FactorOutcome:
    """All 0/1 factors A of P with P = A A^* (up to mirror image), decoded and verified against `inst`.

    P       : fmpz_poly (palindromic, nonzero constant term) or fmpz_mpoly (no variable divides P)
    decode  : exponent list of a 0/1 candidate (ints or tuples, ascending) -> point list, or None if impossible
    """
    uni = isinstance(P, flint.fmpz_poly)
    if uni and not all_solutions and P.degree() > COARSE_DEGREE:
        t = tower(P, deadline)
        if t is not None:
            h, s, R = t
            Gs = zero_one_factors(h, deadline, 1)
            Rs = zero_one_factors(R, deadline, 1) if Gs else None
            for r in Rs or []:
                m = max(r)
                for rr in (r, sorted(m - x for x in r)):
                    for g in Gs:
                        pts = decode(sorted(x + s * y for x in rr for y in g))
                        if pts is not None and verify(pts, inst):       # verified: done (a tower NO proves nothing)
                            return FactorOutcome(complete=True, solutions=[canonical(pts)],
                                                 work={"route": "tower", "tower_s": s, "tower_deg": h.degree()})
    if deadline is not None and deadline.expired():
        return FactorOutcome(timed_out=True, work={"why": "deadline"})
    fac = factor(P, deadline)
    return _enumerate(P, fac, n, decode, inst, all_solutions, deadline)


def _enumerate(P, fac, n, decode, inst, all_solutions, deadline) -> FactorOutcome:
    out = FactorOutcome()
    uni = isinstance(P, flint.fmpz_poly)
    rev = _rev_uni if uni else _rev_multi
    if fac is None:
        out.timed_out = True
        out.work["why"] = "deadline during factorization"
        return out
    content, facs = fac
    out.work["factors"] = sum(e for _, e in facs)
    out.work["distinct_factors"] = len(facs)
    if abs(content) != 1:                         # P = A A^* with A 0/1 forces P primitive (Gauss)
        out.complete = True
        out.work["why"] = "content"
        return out

    # group into self-reciprocal factors and reciprocal pairs
    one = flint.fmpz_poly([1]) if uni else P.context().constant(1)
    base = one
    pairs: list[tuple] = []                            # (f, f^*, e)
    used = [False] * len(facs)
    for i, (f, e) in enumerate(facs):
        if used[i]:
            continue
        used[i] = True
        fr = rev(f)
        if fr == f or fr == -f:
            if e % 2:
                out.complete = True
                out.work["why"] = "odd self-reciprocal multiplicity"
                return out
            base *= f ** (e // 2)
            continue
        j = next((j for j in range(i + 1, len(facs)) if not used[j] and (facs[j][0] == fr or facs[j][0] == -fr)), None)
        if j is None or facs[j][1] != e:               # cannot happen for a palindromic P; treat as not A A^*
            out.complete = True
            out.work["why"] = "unpaired factor"
            return out
        used[j] = True
        pairs.append((f, fr, e))
    out.work["pairs"] = len(pairs)
    total = 1
    for _, _, e in pairs:
        total *= e + 1
    out.work["tuples"] = (total + 1) // 2

    # options per pair: O_j = f^j (f^*)^(e-j), with exact reversal rev(O_j) = (f^*)^j f^(e-j)
    opts = []
    for f, fr, e in pairs:
        opts.append([f ** j * fr ** (e - j) for j in range(e + 1)])
    base_rev = rev(base)
    if uni:
        t = min(LOW, P.degree() // 2 + 1)
        lo = lambda g: g.truncate(t)
        lopts = [[lo(o) for o in row] for row in opts]
        lbase, lbase_rev = lo(base), lo(base_rev)

    seen: set = set()
    tried = screened = 0
    s = len(pairs)
    choice = [0] * s

    def accept(A) -> bool:
        """Leaf with 0/1 coefficients up to sign: decode, verify, record.  Returns True to stop.
        inst None: no decoding or instance check, record the exponent set itself."""
        if inst is None:
            pts = _exponents(A)
        else:
            pts = decode(_exponents(A))
            if pts is None or not verify(pts, inst):
                return False
        c = tuple(canonical(pts))
        if c not in seen:
            seen.add(c)
            out.solutions.append(list(c))
        return not all_solutions

    def leaf_uni(la, lr) -> bool:
        nonlocal screened
        for g in (la, lr):
            co = g.coeffs()
            sgn = 1 if co[0] > 0 else -1
            if any(c != 0 and c != sgn for c in co):
                return False
        screened += 1
        A = base
        for i in range(s):
            A *= opts[i][choice[i]]
        return _is_one_bit(A) and accept(A)

    def dfs(i, tied, prod, la, lr) -> bool:
        """Depth-first over pairs; `tied`: prefix equals its mirror so far (enforces j <= e - j lexicographically)."""
        nonlocal tried
        if deadline is not None and deadline.expired():
            raise _Expired
        if i == s:
            tried += 1
            if uni:
                return leaf_uni(la, lr)
            return _is_one_bit(prod) and accept(prod)
        e = pairs[i][2]
        for j in range(e + 1):
            if tied and j > e - j:
                break
            choice[i] = j
            if uni:
                nla = la.mul_low(lopts[i][j], t)
                nlr = lr.mul_low(lopts[i][e - j], t)
                stop = dfs(i + 1, tied and j == e - j, None, nla, nlr)
            else:
                stop = dfs(i + 1, tied and j == e - j, prod * opts[i][j], None, None)
            if stop:
                return True
        return False

    try:
        if uni:
            dfs(0, True, None, lbase, lbase_rev)
        else:
            dfs(0, True, base, None, None)
        out.complete = True
    except _Expired:
        out.timed_out = True
    out.work["candidates"] = tried
    if uni:
        out.work["screened"] = screened
    return out


class _Expired(Exception):
    pass
