"""Compiled (numba) helper for stages/lift.py: all 3-term relations S[i] + S[j] = S[k] (i <= j) by a two-pointer
sweep, O(|S|^2) simple steps and no temporaries.  Imported lazily; lift.py has a numpy fallback with the same output
set (order may differ)."""
import numba as _nb
import numpy as np


@_nb.njit(cache=True)
def _sweep(s, W, out, fill):
    M = s.shape[0]
    cnt = 0
    for i in range(M):
        if 2 * s[i] > W:
            break
        k = i
        for j in range(i, M):
            t = s[i] + s[j]
            if t > W:
                break
            while s[k] < t:
                k += 1
            if s[k] == t:
                if fill:
                    out[cnt, 0] = i
                    out[cnt, 1] = j
                    out[cnt, 2] = k
                cnt += 1
    return cnt


def three_term_nb(s):
    W = s[-1]
    dummy = np.zeros((1, 3), np.int64)
    c = _sweep(s, W, dummy, False)
    out = np.empty((c, 3), np.int64)
    _sweep(s, W, out, True)
    return out


@_nb.njit(cache=True)
def violations(U, ia, ib, ic, maxbad):
    """Indices t (at most maxbad, ascending) with U[ia[t]] + U[ib[t]] != U[ic[t]] (rows of the int64 matrix U)."""
    T = ia.shape[0]
    r = U.shape[1]
    out = np.empty(maxbad, np.int64)
    c = 0
    for t in range(T):
        a = ia[t]
        b = ib[t]
        k = ic[t]
        for q in range(r):
            if U[a, q] + U[b, q] != U[k, q]:
                out[c] = t
                c += 1
                break
        if c == maxbad:
            break
    return out[:c]
