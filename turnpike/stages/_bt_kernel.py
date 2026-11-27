"""Compiled (numba) version of the backtracking search in backtrack.py; same algorithm, iterative DFS.
Imported lazily: loading numba costs about half a second.

Speed tricks (same answers and node counts as the Python version):
  * distance -> index lookups through a bucket index: the m distinct distances are spread over m buckets by value
    (bucket = floor(d * m / (W + 1))), so a lookup scans about one entry, with 4 bytes of index per distance
    (buckets holding more than 8 entries, from clustered distances, are binary searched);
  * the distances consumed along the current path live in one flat stack (size n(n-1)/2), not an n x n array.
"""
import numba as _nb
import numpy as np

@_nb.njit(cache=True)
def _bucket_index(vals, W):
    """start[b] = first index j with bucket(vals[j]) >= b, for bucket(d) = int(d * scale), scale = m / (W + 1)."""
    m = vals.shape[0]
    scale = m / (W + 1.0)
    start = np.zeros(m + 2, np.int32)
    for j in range(m):
        start[np.int64(vals[j] * scale) + 1] += 1
    for b in range(m + 1):
        start[b + 1] += start[b]
    return start, scale


@_nb.njit(cache=True)
def _search_nb(vals, cnt0, n, W, node_limit, want_all, max_sols):
    """Iterative DFS. vals: distinct distances ascending (int64); cnt0: multiplicities.
    Returns (status, nodes, sols[n_found, n]): status 0 = exhausted, 1 = stopped at first solution, 2 = node limit,
    3 = solution buffer full (want_all only; the list is incomplete and must not be used)."""
    m = vals.shape[0]
    start, scale = _bucket_index(vals, W)

    cnt = cnt0.copy()
    cnt[m - 1] -= 1                              # W is the largest distance
    placed = np.zeros(n, np.int64)
    placed[1] = W
    npl = 2
    total = n * (n - 1) // 2
    ptr_st = np.zeros(n + 1, np.int64)           # per depth: scan pointer into vals (descending)
    ch_st = np.zeros(n + 1, np.int64)            # per depth: next choice (0: y, 1: W - y)
    base = np.zeros(n + 2, np.int64)             # per depth: start of its used distances in `used`
    used = np.zeros(total + 1, np.int64)         # flat stack of consumed distance indices along the path
    sols = np.zeros((max_sols, n), np.int64)
    nsol = 0
    nodes = 0
    depth = 0
    ptr_st[0] = m - 1
    top = 0
    while depth >= 0:
        p_ = ptr_st[depth]
        while p_ >= 0 and cnt[p_] == 0:
            p_ -= 1
        ptr_st[depth] = p_
        if p_ < 0:                               # all distances explained: a realization
            if nsol == max_sols:                 # buffer full: stop, the caller reruns without a cap
                return 3, nodes, sols[:nsol]
            for i in range(n):
                sols[nsol, i] = placed[i]
            nsol += 1
            if not want_all:
                return 1, nodes, sols[:nsol]
            depth -= 1
            if depth >= 0:                       # undo the placement made at this depth
                npl -= 1
                for t in range(base[depth], top):
                    cnt[used[t]] += 1
                top = base[depth]
            continue
        y = vals[p_]
        nch = 1 if (depth == 0 or y == W - y) else 2
        c = ch_st[depth]
        if c >= nch:                             # both choices tried: backtrack
            ch_st[depth] = 0
            depth -= 1
            if depth >= 0:
                npl -= 1
                for t in range(base[depth], top):
                    cnt[used[t]] += 1
                top = base[depth]
            continue
        ch_st[depth] = c + 1
        pt = y if c == 0 else W - y
        nodes += 1
        if nodes == node_limit:
            return 2, nodes, sols[:nsol]
        base[depth] = top
        ok = True
        for i in range(npl):
            q = placed[i]
            d = pt - q if pt > q else q - pt
            j = -1
            bk = np.int64(d * scale)
            lo = start[bk]
            hi = start[bk + 1]
            while hi - lo > 8:                   # clustered distances (wide products): binary search the bucket
                mid = (lo + hi) >> 1
                if vals[mid] <= d:
                    lo = mid
                else:
                    hi = mid
            for jj in range(lo, hi):
                v = vals[jj]
                if v == d:
                    j = jj
                    break
                if v > d:
                    break
            if j < 0 or cnt[j] == 0:
                ok = False
                break
            cnt[j] -= 1
            used[top] = j
            top += 1
        if not ok:
            for t in range(base[depth], top):
                cnt[used[t]] += 1
            top = base[depth]
            continue
        placed[npl] = pt
        npl += 1
        depth += 1
        ptr_st[depth] = p_
        ch_st[depth] = 0
    return 0, nodes, sols[:nsol]
