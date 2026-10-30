# Development Log — Bugs Found and Fixed

The errors found while building and validating the solver, what gave them
away, and which results they invalidate. Kept because several of them only
appear in one scenario, and because results computed before a fix must not
be reused.

---

## 1. Core solver (found on the lid-driven cavity)

| Bug | Symptom | Fix | Where |
|---|---|---|---|
| Wrong sign of the PPE right-hand side | divergence / NaN | $b = -(\rho/\Delta t)\,(\nabla\cdot\mathbf{u}^*)\,V$, consistent with the positive-definite matrix | `PressurePoissonSystem::computeRHS` |
| Matrix sign convention | CG fails (matrix not positive-definite) | diagonal $+\sum a_{nb}$, off-diagonal $-a_{nb}$ | `LaplaceOperator::assembleToTriplets` |
| Reference pressure pinned on the diagonal only | non-symmetric matrix, CG breaks down | zero row **and** column 0, diagonal 1, `b[0] = 0` every step | `LaplaceOperator::pinReference`, [main.cpp:378](../../src/main.cpp#L378) |
| Ghost value $u_{ghost} = u_{wall}$ at walls | wrong wall condition (diffusion stencil off by a factor 4/3 at the first node) | image method $u_{ghost} = 2u_{wall} - u_P$ with the wall distance doubled | `MomentumSolver::uRHS/vRHS`, `DiffusionOperator` |
| Zone corner nodes inside the domain with zero cells | zero-width zones, $S/\Delta = 0/0$ = NaN in the mesh | for single-zone cavities put every `node_*` at $(L, H)$ | `input/*/data.txt` |
| $v$ on the north wall set from the lid value | spurious flux through the lid | normal velocity on every wall is 0; the lid value only enters the tangential ghost of $u$ | `applyBoundaryConditions` |

Performance choice from the same period: the CG pressure solve is
warm-started from the previous step's pressure.

---

## 2. Immersed obstacle (found validating the square cylinder, October 2026)

These only show up with an immersed obstacle; the cavity outputs are
byte-identical before and after the fixes, and `operator_verify` reproduces
the same convergence data. **Cylinder results computed before these fixes
are invalid.**

| Bug | Symptom | Fix | Where |
|---|---|---|---|
| Tangential ghost at an inlet used the inflow speed | wrong shear next to the inlet | `ghostValue(INLET)` returns 0 | [MomentumSolver.cpp:35](../../src/MomentumSolver.cpp#L35) |
| Viscous force components swapped | non-zero mean lift on a symmetric body ($C_L = 9.3$ at step 1), drag missing its friction part | $\mu\,\partial u/\partial y$ on top/bottom faces → drag; $\mu\,\partial v/\partial x$ on front/rear faces → lift | [ObstacleForces.cpp](../../src/postprocessing/ObstacleForces.cpp) |
| Obstacle-surface velocities entered the pressure equation | $\lvert\nabla\cdot\mathbf{u}\rvert \approx U/D$ in the cells next to the body | BCs re-applied right after `advance()`, before the PPE; `advance()` skips surface and solid nodes | [main.cpp:372](../../src/main.cpp#L372), [MomentumSolver.cpp:359](../../src/MomentumSolver.cpp#L359) |
| Inverted obstacle wall distances | image distance $D + h$ instead of $h$ → nearly free-slip cylinder; friction drag decreasing (even negative) under refinement; spurious leading-edge separation at $Re = 40$ | a solid neighbour to the west means the wall is the obstacle's **east** face (`obstacle_xmax`), etc. | `DiffusionOperator`, `MomentumSolver::uRHS/vRHS` |
| Diverged run reported as converged | NaN fields written as a "steady" result: `std::max(du, NaN)` drops the NaN | explicit `std::isfinite` check; the case stops with an error and writes nothing | [main.cpp:403-420](../../src/main.cpp#L403-L420) |
| Velocities inside the solid left at non-physical values | polluted exported fields (never used by the stencils) | held at 0 in `applyBoundaryConditions` | [MomentumSolver.cpp:139-145](../../src/MomentumSolver.cpp#L139-L145) |

After these fixes the steady $Re = 40$ case matches Breuer et al.
($C_D$ 1.716 vs 1.708, $L_r/D$ 2.143 vs 2.151) —
[../scenarios/square_cylinder.md §8](../scenarios/square_cylinder.md#8-validation-results-square_cylinder_breuer).

---

## 3. Lessons for future changes

- **Check the explicit time-step limit** before refining a grid or raising
  `dt`: $\Delta t \lesssim h_{\min}^2/(8\nu)$ for AB2 in 2D, with $\nu$
  replaced by $\max(Pr, 1)$ in the heated cavity
  ([../theory/equations.md §6](../theory/equations.md#6-stability-of-the-explicit-scheme)).
- **Symmetric problems are a free test**: a non-zero mean lift on the
  cylinder, or an asymmetric steady wake, means a bug, not physics.
- **Grid refinement must improve results monotonically**; a quantity that
  drifts away from the reference under refinement (here, the friction drag)
  points to a boundary-condition error.
- **Re-run the cavity cases after touching shared code** and compare the
  outputs byte by byte: they exercise none of the obstacle branches.
