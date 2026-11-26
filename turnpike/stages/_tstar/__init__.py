"""Internals of the T*_k stage (section 5): propagation + exact face LP search.

Modules
  lp         HiGHS wrappers and the exact Farkas certificate check (`certify`); the warm LP oracle of close_k (`FaceLP`).
  propagate  compiled exact 0/1 bound propagation (plain, and with aliasing of pair columns).
  tree       the T* class tree: candidate classes, forced closure, T*'s exact node test on the SA_2 face system.
  closure    close_k^LP: k-level probing over the face system + Q rows, with the certified LP as the base refutation.
  search     the DFS on the T*_k tree without backtracking (strict), reporting (ST_k^LP+) violations.

Soundness rule kept throughout: a node or a leg is refuted only by an exact rule, an exact propagation, or a Farkas
certificate verified in integer arithmetic.  A floating-point LP answer can only make the search prune LESS.
"""
from .search import Search, Violation, run_with_big_stack

__all__ = ["Search", "Violation", "run_with_big_stack"]
