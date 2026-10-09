# The Solver Pipeline — Class by Class

Theory notes on how the discrete equations ([equations.md](equations.md)) are
implemented on the staggered grid ([mesh.md](mesh.md)) as a sequence of C++
classes, each owning one piece of the per-timestep pipeline.

---

## 0. Overview: the per-timestep pipeline

```
                      ┌───────────────────────────────────────────────────────┐
 u^n, v^n            │ 1. MomentumSolver::advance            (predictor)      │
 ────────────────────→│    u* = u^n + dt·[ν∇²u − (u·∇)u]    (AB2)            │
                      │                                                       │
 θ^n (DHC only)       │ 2. addBuoyancy                        (DHC only)      │
 ────────────────────→│    v* += dt·Ra·Pr·θ                                   │
                      │                                                       │
                      │ 3. PressurePoissonSystem + CG         (projection)    │
                      │    b = (ρ/dt)·∇·u*                                    │
                      │    solve  ∇²p = b   via Conjugate Gradient            │
                      │    u^{n+1} = u* − (dt/ρ)·∇p                           │
                      │                                                       │
                      │ 4. EnergySolver::advance              (DHC only)      │
                      │    θ^{n+1} = θ^n + dt·[α∇²θ − (u^{n+1}·∇)θ]  (AB2)   │
                      └───────────────────────────────────────────────────────┘
                                          │
                                          └─→ u^{n+1}, v^{n+1}, p, θ^{n+1}
```

- **Order matters.** Buoyancy is injected into $v^*$ *after* the momentum
  predictor and *before* the pressure solve, so the projection step makes the
  **final** field divergence-free **including** the buoyant forcing. The
  energy equation is advanced *last*, using the already-corrected,
  divergence-free velocity — keeping scalar transport consistent with a
  solenoidal flow field.
- This ordering is also why the DHC extension (step 2 and step 4) never had
  to modify `MomentumSolver` or `PressurePoissonSystem` — see §5.
- Driving loop: [src/main.cpp:225-271](../../src/main.cpp#L225-L271).

---

## 1. Discrete operators — the shared building blocks

Four small, independently-testable classes implement the discrete spatial
operators, each used by one or more solvers above:

| Operator | Maps | Discretization | Class |
|---|---|---|---|
| Laplacian $\nabla^2$ | P → P | FVM, 2nd-order (2nd-order only on **uniform** meshes — see caveat below) | [`LaplaceOperator`](../../src/operators/LaplaceOperator.hpp) |
| Gradient $\nabla$ | P → (U,V) | FD, 2nd-order centered | [`GradientOperator`](../../src/operators/GradientOperator.hpp) |
| Divergence $\nabla\cdot$ | (U,V) → P | FVM, 2nd-order (outward-flux integration) | [`DivergenceOperator`](../../src/operators/DivergenceOperator.hpp) |
| Diffusion $\nu\nabla^2$ | (U,V) → (U,V) | FD, 2nd-order + ghost-cell walls | [`DiffusionOperator`](../../src/operators/DiffusionOperator.hpp) |

- `GradientOperator` and `DivergenceOperator` are purely geometric — no
  `dt`/`ρ` factors. `PressurePoissonSystem` applies those factors itself
  when assembling the RHS and projecting (§3).
- `DiffusionOperator` is reused directly by `MomentumSolver` (§2); its ghost
  values are supplied once at construction from the boundary data.
- **Verification:** [`operator_verify`](../../src/verification/verify_operators.cpp)
  checks every operator against a Method of Manufactured Solutions (MMS):
  a known analytic field is differentiated by hand, the discrete operator is
  applied to its sampled values, and the error is measured at increasing mesh
  resolution. `plotters/plot_convergence.py` turns
  `output/convergence.csv` into `output/convergence_rates.png`, expecting a
  slope of ≈2.0 (2nd order) in log-log error-vs-$h$ plots.
- **Known limitation (stretched meshes):** on a **non-uniform** (stretched)
  mesh, `LaplaceOperator`'s stencil `a_E = Se/(x[i+1]-x[i])` degrades to
  **first-order** accuracy (measured rate ≈1.0 instead of ≈2.0 — see
  `docs/OPERATORS_REPORT.md`, Suite 10). The fix (using `2/(h_e+h_w)` instead
  of `1/dx`) is identified but not yet applied to the production stencil —
  worth knowing before trusting convergence-order claims on highly stretched
  production grids.

---

## 2. `MomentumSolver`: the predictor step

[src/MomentumSolver.hpp](../../src/MomentumSolver.hpp) /
[.cpp](../../src/MomentumSolver.cpp) — implements
$u^* = u^n + \Delta t\,[\nu\nabla^2 u - (u\cdot\nabla)u]$ (and the `v`
analogue).

- **Construction** ([MomentumSolver.cpp:51](../../src/MomentumSolver.cpp#L51)):
  reads the wall boundary data once and derives the four ghost values used by
  every wall-adjacent stencil for the rest of the run
  (`u_ghost_north/south`, `v_ghost_west/east`) via `ghostValue`/`extractGhost`
  ([MomentumSolver.cpp:28-50](../../src/MomentumSolver.cpp#L28-L50)). A
  `DiffusionOperator` is constructed with these same ghosts and owned as a
  member, so the diffusive part of the stencil is never duplicated.
- **Diffusion:** delegated entirely to `DiffusionOperator` (§1) — 2nd-order
  finite differences on the non-uniform grid, ghost-cell image method at
  walls.
- **Convection:** selectable at construction,
  `ConvectionScheme::CDS` (central, 2nd order, can oscillate at high cell
  Peclet number) or `ConvectionScheme::UPWIND` (1st order, more diffusive but
  unconditionally stable with respect to convection).
- **`uRHS`/`vRHS`** ([MomentumSolver.cpp:119](../../src/MomentumSolver.cpp#L119),
  [:195](../../src/MomentumSolver.cpp#L195)): compute diffusion − convection
  at one interior node; only interior nodes are updated by `advance` —
  boundary nodes retain their prescribed values.
- **Time advance** ([MomentumSolver.cpp:273](../../src/MomentumSolver.cpp#L273)):
  Adams–Bashforth 2, falling back to Euler on the very first call (no
  previous-step RHS buffer yet) — see [equations.md §3](equations.md#3-time-integration).
- **`applyBoundaryConditions`** must be called once before the first
  `advance()` and again after every projection step, so boundary velocity
  nodes always carry their prescribed values (walls, moving lid) rather than
  stale predictor values.

---

## 3. `PressurePoissonSystem`: the elliptic constraint

[src/PressurePoissonSystem.hpp](../../src/PressurePoissonSystem.hpp) /
[.cpp](../../src/PressurePoissonSystem.cpp) — assembles and solves
$\nabla^2 p = (\rho/\Delta t)\nabla\cdot u^*$, then projects.

Per-step usage (also see §0 diagram):
```cpp
poisson.setVelocityFields(u_star, v_star);
auto b = poisson.computeRHS(dt, rho);
b[0] = 0.0;                                   // keep pinned row consistent
p_sol = CG_solver.solve(A, b, p_sol);         // warm-started, see §4
poisson.projectVelocity(p_sol, dt, rho, u_star, v_star);
```

- **`assembleMatrix`** ([PressurePoissonSystem.cpp:105](../../src/PressurePoissonSystem.cpp#L105)):
  builds the discrete Laplacian as a list of `Triplet{row, col, value}`
  (Eigen-compatible layout, see the comment in
  [PressurePoissonSystem.hpp](../../src/PressurePoissonSystem.hpp)), delegating
  the actual FVM stencil to `LaplaceOperator` (§1). Built **once** per case
  (the matrix doesn't change — only the RHS does, every step).
- **Singular system / pinning**
  ([PressurePoissonSystem.cpp:119](../../src/PressurePoissonSystem.cpp#L119)):
  with homogeneous Neumann conditions on every wall (a closed cavity has no
  pressure outlet), the discrete Laplacian is singular — constant pressure
  solves $\nabla^2 p = 0$ trivially, so the matrix has a 1D null space.
  `pinReferencePressure` regularizes this by forcing row **and column** 0 to
  the identity (`p[0] = 0`), not just the diagonal — zeroing only the
  diagonal would leave the matrix inconsistent with a zeroed RHS and was a
  real bug fixed during development (see `PROJECT_SUMMARY.md`).
- **`computeRHS`** ([PressurePoissonSystem.cpp:147](../../src/PressurePoissonSystem.cpp#L147)):
  integrates the predictor's divergence over each cell by outward-normal
  flux summation (FVM), reusing the same staggered index mapping as `mesh.md §2`:
  $$
  b_P = \frac{\rho}{\Delta t}\Big[u^*_E S_e - u^*_W S_w + v^*_N S_n - v^*_S S_s\Big]
  $$
- **`projectVelocity`** ([PressurePoissonSystem.cpp:190](../../src/PressurePoissonSystem.cpp#L190)):
  applies the correction
  $u_{next} = u^* - (\Delta t/\rho)(p_E - p_W)/\text{dist}_{EW}$ (and the `v`
  analogue) at interior velocity nodes only; boundary nodes are left
  untouched, preserving their prescribed values.

---

## 4. The linear solver: warm-started Conjugate Gradient

- `CG` (in [src/solver/sparse/](../../src/solver/sparse/solver_class_sparse.hpp))
  solves the symmetric positive-definite system `A·p = b` built above.
  SPD is guaranteed by construction (`LaplaceOperator`'s FVM stencil
  produces a symmetric, diagonally-dominant-after-pinning matrix), so CG
  converges without needing a general-purpose (non-SPD) solver.
- **Warm start:** `CG_solver.solve(A, b, p_sol)` is called with the
  **previous time step's pressure solution** as the initial guess, every
  step. Since the pressure field typically changes only slightly from one
  step to the next (especially near steady state), this drastically reduces
  the number of CG iterations needed per step compared to starting from
  zero — this is a deliberate performance choice.
- `CSRMatrix::fromTriplets` converts the assembled triplets into
  compressed-sparse-row storage once per case, used by every `multiply(A,v)`
  inside CG for the remainder of that run.
- Other solvers (`PCG`, `GS`, `SOR`) exist in the same file for
  experimentation/comparison but the production loop in `main.cpp` uses
  plain `CG`.

---

## 5. `EnergySolver` + `Buoyancy.hpp`: bolting on a new transported scalar

This pair is the template for **extending** the solver with any new
advected/diffused scalar without touching the momentum/pressure core — worth
reading even if you don't care about heat transfer specifically.

- **[`EnergySolver`](../../src/EnergySolver.hpp)** mirrors `MomentumSolver`'s
  structure almost line for line, but for a scalar living on the **pressure**
  grid (cell-centered, like $p$) rather than on a velocity face:
  - Construction ([EnergySolver.cpp:17](../../src/EnergySolver.cpp#L17)):
    reads per-wall temperature BC type/value from `pMesh.bound_*_T`.
  - Convective flux ([`tRHS`](../../src/EnergySolver.cpp#L48)): because
    $\theta$ sits on P and the staggered velocities already sit exactly on P's
    faces, the face-normal velocity needed for the flux requires **no
    interpolation** — `u[uIdx(i+1,j)]` is already the east-face velocity of
    cell $(i,j)$. Upwind (default, recommended at high Rayleigh number) or
    CDS, same as momentum.
  - Diffusive term: 2nd-order FD with the same ghost-cell image method as
    `DiffusionOperator` (§1), but with isothermal/adiabatic ghost formulas
    instead of velocity ones (see [equations.md §4](equations.md#4-boundary-conditions--the-ghost-cell-image-method)).
  - Time advance ([`advance`](../../src/EnergySolver.cpp#L118)): identical AB2
    scheme and Euler-start fallback as `MomentumSolver::advance`.
- **[`Buoyancy.hpp`](../../src/Buoyancy.hpp)** (`addBuoyancy`, header-only):
  a free function, not a class — called explicitly from the time loop between
  the momentum predictor and the pressure solve (§0). It:
  1. Linearly interpolates $\theta$ (on P) to each interior V-node (between
     the two adjacent P cells) since buoyancy acts on the $v$-velocity, which
     lives on a different grid than $\theta$.
  2. Adds `dt·Ra·Pr·θ_v` to that V-node — **only interior nodes** (`l=2..M`),
     leaving wall (no-penetration) V-nodes untouched.
- **Why this design works without touching the core:** `MomentumSolver` and
  `PressurePoissonSystem` have no knowledge of buoyancy or temperature at
  all — buoyancy is injected as a plain increment on an already-public
  velocity vector, and the pressure projection downstream automatically
  re-enforces incompressibility on the buoyancy-modified field. To add a
  different transported scalar (concentration, salinity, …), the same
  three-piece pattern applies: (1) a cell-centered AB2 scalar solver mirroring
  `EnergySolver`, (2) a small header-only coupling function mirroring
  `Buoyancy.hpp` that injects its source term into the predictor before the
  pressure solve, (3) a one-line call in the time loop — no change to
  `MomentumSolver` or `PressurePoissonSystem`.

---

## 6. Stability: this is a fully explicit scheme

- Both the momentum predictor and the energy update are **explicit** in time
  (no implicit diffusion solve) — only the pressure Poisson equation is
  solved implicitly (as an elliptic constraint, not a time-marched implicit
  diffusion term).
- This caps the usable time step by the explicit diffusion stability limit:
  $$
  \Delta t \lesssim \frac{\Delta x^2}{2\,\max(\alpha,\nu)}
  $$
- For the lid-driven cavity ($\nu \sim \text{Re}^{-1}$, moderate $Re$),
  typical production time steps (`dt ~ 1e-3`) are comfortably stable.
- For the Differentially Heated Cavity, `nu = Pr` (often $O(1)$) and
  $\alpha=1$ make diffusion much stiffer relative to the grid spacing, so
  production runs need `dt ~ 1e-5` — using the lid-cavity-scale `dt`
  diverges to NaN (the CG solve fails to converge first). See
  [../scenarios/differentially_heated_cavity.md](../scenarios/differentially_heated_cavity.md)
  for concrete numbers and validated settings.

---

## Key files quick-reference

| Piece | Header | Implementation |
|---|---|---|
| Mesh generation | [src/mesh/Mesh_config.h](../../src/mesh/Mesh_config.h) | [src/mesh/MeshGenerator.cpp](../../src/mesh/MeshGenerator.cpp) |
| Operators | `src/operators/*.hpp` | `src/operators/*.cpp` |
| Momentum predictor | [src/MomentumSolver.hpp](../../src/MomentumSolver.hpp) | [src/MomentumSolver.cpp](../../src/MomentumSolver.cpp) |
| Pressure / projection | [src/PressurePoissonSystem.hpp](../../src/PressurePoissonSystem.hpp) | [src/PressurePoissonSystem.cpp](../../src/PressurePoissonSystem.cpp) |
| Linear solver | [src/solver/sparse/solver_class_sparse.hpp](../../src/solver/sparse/solver_class_sparse.hpp) | [src/solver/sparse/solver_sparse.cpp](../../src/solver/sparse/solver_sparse.cpp) |
| Energy equation | [src/EnergySolver.hpp](../../src/EnergySolver.hpp) | [src/EnergySolver.cpp](../../src/EnergySolver.cpp) |
| Buoyancy coupling | [src/Buoyancy.hpp](../../src/Buoyancy.hpp) | header-only |
| Main driver / time loop | — | [src/main.cpp](../../src/main.cpp) |
| Operator verification (MMS) | — | [src/verification/verify_operators.cpp](../../src/verification/verify_operators.cpp) |
