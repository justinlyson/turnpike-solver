"""DFS on the T*_k tree without backtracking (Theorem A-ST_k^LP, PROOF_STPLUS 8.12).

At an entered node with branching class c, both children ("c" and "W - c"; only "c" at the first branching class, by
mirror symmetry) are formed with their closure and tested:

  1. T*'s node test (tree.ClassTree.test: exact rules, propagation, certified face LP);
  2. if at least TWO children pass step 1 and k >= 1: close_k^LP on each of them (closure.CloseK).

Step 2 is skipped when at most one child passes step 1.  This does not change the outcome: at a node that extends to a
realization some child does too, and a sound test never refutes it, so a single passing child is that child; at a node
that does not extend, entering the child only postpones the refutation to its own subtree (finite, every leaf exact).

The children that pass are the node's SURVIVORS; the search enters the first one.  The conjecture (ST_k^LP+) says every
first-wrong child (a child that does not extend, of a node that does) is refuted, so on a node that extends the first
survivor extends too and the search reaches a realization after at most 1 + 2(n - 1) child tests.  If an entered
survivor's subtree dies while the node had another survivor, the search would have to backtrack: that is reported as
a violation of (ST_k^LP+) (or of the list having no realization while a wrong child survived) and the search stops.

Outcomes:  YES - a realization (verified by the caller);  NO - the root, or every survivor chain, was refuted with no
survivor ever skipped (every refutation is exact, so this is a proof);  VIOLATION - see above.

`backtrack=True` turns the search into the ordinary complete DFS of the research code (exponential worst case).  It is
used only by the equivalence tests; the stage always runs with backtrack=False.
"""
from __future__ import annotations

import sys
import threading
from collections import Counter

from .closure import CloseK
from .tree import ClassTree


class Violation(Exception):
    """An entered survivor died while its node had another survivor (the strict DFS would have to backtrack)."""


class Search:
    def __init__(self, counts, n, k, deadline=None, all_solutions=False, backtrack=False):
        self.k = int(k)
        self.all = bool(all_solutions)
        self.backtrack = bool(backtrack)
        self.deadline = deadline
        self.st = Counter()
        self.T = ClassTree(Counter(counts), n, deadline, self.st)
        self.T.keep_face = self.k >= 1
        self.solutions = []
        self.first_done = False

    # ------------------------------------------------------------------ one node
    def _children(self, ci, sides):
        """Test the children of the node at branching class ci.  Returns (leaf realization or None, survivors)."""
        T, st = self.T, self.st
        passed = []
        for side in sides:
            st["child_tests"] += 1
            s, ci2, trail = T.apply(ci, side)
            if s == "complete":
                X = T.leaf()
                if X is not None:
                    if not self.all:
                        T.undo(trail)
                        return X, []
                    self.solutions.append(X)
            elif s == "branch":
                if T.test(ci2):
                    passed.append((side, T.last_face))
            T.undo(trail)
        if self.k >= 1 and len(passed) >= 2:
            keep = []
            for side, face in passed:
                if face is None:                       # passed before the face stage (nothing left to decide)
                    keep.append((side, face))
                    continue
                st["close_tests"] += 1
                ck = CloseK(T, face, self.deadline)
                dead = ck.refutes(self.k)
                st["close_lp_calls"] += ck.lp.calls
                st["close_lp_unknown"] += ck.lp.unknown
                st["close_propagations"] += ck.node.calls
                if dead:
                    st["close_refuted"] += 1
                else:
                    keep.append((side, face))
            passed = keep
        return None, [side for side, _ in passed]

    def _node(self, ci):
        """Search below an entered node.  Returns True if a realization was found below it (in all-solutions mode: if
        at least one was), False if the subtree died."""
        T, st = self.T, self.st
        st["nodes"] += 1
        if self.first_done:
            sides = (0, 1)
        else:
            sides = (0,)
            self.first_done = True
        n_before = len(self.solutions)
        X, surv = self._children(ci, sides)
        if X is not None:
            self.solutions.append(X)
            return True
        leaves = len(self.solutions) - n_before          # complete children (all-solutions mode only)
        found_any = leaves > 0
        for side in surv:
            s, ci2, trail = T.apply(ci, side)
            try:
                ok = self._node(ci2)
            finally:
                T.undo(trail)
            if ok:
                found_any = True
                if not self.all:
                    return True
            elif not self.backtrack and len(surv) + leaves >= 2:
                raise Violation()
        return found_any

    # ------------------------------------------------------------------ the whole search
    def run(self):
        """Returns 'YES' (self.solutions filled), 'NO' (proved) or 'VIOLATION'."""
        T, st = self.T, self.st
        if T.D[T.W] != 1:
            return "NO"
        trail = []
        for p in (0, T.W):
            a = T.place(p)
            if a is None:
                return "NO"
            trail.append((p, a))
        s, ci = T.closure(0, trail)
        if s == "dead":
            return "NO"
        if s == "complete":
            X = T.leaf()
            if X is None:
                return "NO"
            self.solutions.append(X)
            return "YES"
        if not T.test(ci):
            return "NO"
        try:
            ok = self._node(ci)
        except Violation:
            return "VIOLATION"
        return "YES" if ok else "NO"


def run_with_big_stack(fn, *a):
    """The search recurses once per branching class (up to about n levels): run it on a thread with a 256 MB stack."""
    box = {}
    if sys.getrecursionlimit() < 10 ** 6:
        sys.setrecursionlimit(10 ** 6)
    old = threading.stack_size(256 * 1024 * 1024 - 4096)
    try:
        def target():
            try:
                box["r"] = fn(*a)
            except BaseException as e:          # re-raised in the caller's thread
                box["e"] = e
        th = threading.Thread(target=target)
        th.start()
        th.join()
    finally:
        threading.stack_size(old)
    if "e" in box:
        raise box["e"]
    return box["r"]
