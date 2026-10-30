# The Solver Pipeline — Class by Class

Theory notes on how the discrete equations ([equations.md](equations.md)) are
implemented on the staggered grid ([mesh.md](mesh.md)) as a sequence of C++
classes, each owning one piece of the per-timestep pipeline.

---

## 0. Overview: the per-timestep pipeline

```
 u^n, v^n  (θ^n for the DHC)
     │
     ▼
 1. MomentumSolver::advance                  predictor, interior nodes only
       u* = u^n + dt·[3/2 R^n − 1/2 R^{n−1}],   R = ν∇²u − (u·∇)u   (AB2)
     │
 2. MomentumSolver::applyBoundaryConditions  walls, inlet, outlet, obstacle on u*
     │
 3. addBuoyancy                              DHC only:  v* += dt·Ra·Pr·θ^n
     │
 4. PressurePoissonSystem + CG / PCG         projection
       b = −(ρ/dt)·(∇·u*)·V   →   solve A p = b   (A = −∇², SPD)
       u^{n+1} = u* − (dt/ρ)·∇p              interior nodes only
     │
 5. MomentumSolver::applyBoundaryConditions  reset boundary nodes, refresh outlet
     │
 6. EnergySolver::advance                    DHC only:  θ^{n+1} = θ^n + dt·[∇²θ − ∇·(u^{n+1}θ)]  (AB2)
     │
 7. ObstacleForces::computeForces            cylinder only:  C_D, C_L → forces.csv
     │
     ▼
 u^{n+1}, v^{n+1}, p^{n+1}, θ^{n+1}          every 100 steps: ‖Δu‖∞/dt vs ss_tol, NaN check
```

- **Step 2 matters.** The pressure equation must see the prescribed
  boundary fluxes, not the predictor's values on boundary/obstacle nodes —
  see [boundary_conditions.md §8](boundary_conditions.md#8-when-the-conditions-are-applied-in-a-time-step).
- **Order matters.** Buoyancy is injected into $v^*$ *after* the momentum
  predictor and *before* the pressure solve, so the projection step makes the
  **final** field divergence-free **including** the buoyant forcing. The
  energy equation is advanced *last*, using the already-corrected,
  divergence-free velocity — keeping scalar transport consistent with a
  solenoidal flow field.
- This ordering is also why the DHC extension (steps 3 and 6) never had
  to modify `MomentumSolver` or `PressurePoissonSystem` — see §5.
- Driving loop: [src/main.cpp:360-427](../../src/main.cpp#L360-L427). Where every
  boundary condition is enforced: [boundary_conditions.md](boundary_conditions.md).

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
  `output/verification/data/convergence.csv` into
  `output/verification/figures/convergence_rates.{pdf,png}`, expecting a
  slope of ≈2.0 (2nd order) in log-log error-vs-$h$ plots.
- **Known limitation (stretched meshes):** on a **non-uniform** (stretched)
  mesh, `LaplaceOperator`'s stencil `a_E = Se/(x[i+1]-x[i])` degrades to
  **first-order** accuracy (measured rate ≈1.0 instead of ≈2.0 — see
  [operators_verification.md](../reference/operators_verification.md), Suite 10). The fix (using `2/(h_e+h_w)` instead
  of `1/dx`) is identified but not yet applied to the production stencil —
  worth knowing before trusting convergence-order claims on highly stretched
  production grids.

---

## 2. `MomentumSolver`: the predictor step

[src/MomentumSolver.hpp](../../src/MomentumSolver.hpp) /
[.cpp](../../src/MomentumSolver.cpp) — implements
$u^* = u^n + \Delta t\,[\nu\nabla^2 u - (u\cdot\nabla)u]$ (and the `v`
analogue).

- **Construction** ([MomentumSolver.cpp:58](../../src/MomentumSolver.cpp#L58)):
  reads the wall boundary data once and derives the four ghost values used by
  every wall-adjacent stencil for the rest of the run
  (`u_ghost_north/south`, `v_ghost_west/east`) via `ghostValue`/`extractGhost`
  ([MomentumSolver.cpp:26-53](../../src/MomentumSolver.cpp#L26-L53)). A
  `DiffusionOperator` is constructed with these same ghosts and owned as a
  member, so the diffusive part of the stencil is never duplicated.
- **Diffusion:** delegated entirely to `DiffusionOperator` (§1) — 2nd-order
  finite differences on the non-uniform grid, ghost-cell image method at
  walls.
- **Convection:** selectable at construction,
  `ConvectionScheme::CDS` (central, 2nd order, can oscillate at high cell
  Peclet number) or `ConvectionScheme::UPWIND` (1st order, more diffusive but
  unconditionally stable with respect to convection).
- **`uRHS`/`vRHS`** ([MomentumSolver.cpp:156](../../src/MomentumSolver.cpp#L156),
  [:257](../../src/MomentumSolver.cpp#L257)): compute diffusion − convection
  at one interior node; only interior nodes are updated by `advance` —
  boundary nodes retain their prescribed values.
- **Time advance** ([MomentumSolver.cpp:359](../../src/MomentumSolver.cpp#L359)):
  Adams–Bashforth 2, falling back to Euler on the very first call (no
  previous-step RHS buffer yet) — see [equations.md §4](equations.md#4-time-integration).
  Nodes on or inside an immersed obstacle are skipped.
- **`applyBoundaryConditions`** ([MomentumSolver.cpp:84](../../src/MomentumSolver.cpp#L84))
  is called once before the first `advance()`, then twice per step: right
  after `advance()` (so the pressure equation sees the prescribed boundary
  fluxes) and after the projection. What it sets on each side is listed in
  [boundary_conditions.md §3](boundary_conditions.md#3-velocity-normal-component-nodes-on-the-boundary).

---

## 3. `PressurePoissonSystem`: the elliptic constraint

[src/PressurePoissonSystem.hpp](../../src/PressurePoissonSystem.hpp) /
[.cpp](../../src/PressurePoissonSystem.cpp) — assembles and solves
$\nabla^2 p = (\rho/\Delta t)\nabla\cdot u^*$, then projects.

Per-step usage (also see §0 diagram):
```cpp
poisson.setVelocityFields(u_star, v_star);
auto b = poisson.computeRHS(dt, rho);
if (!is_cylinder) b[0] = 0.0;                 // keep pinned row consistent
p_sol = is_cylinder ? PCG_solver.solve(A, b, p_sol)   // warm-started, see §4
                    : CG_solver.solve(A, b, p_sol);
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
  real bug fixed during development (see
  [development_log.md](../reference/development_log.md)). The cylinder
  scenarios do not pin: their Dirichlet outlet row already fixes the constant
  ([boundary_conditions.md §5](boundary_conditions.md#5-pressure)).
- **`computeRHS`** ([PressurePoissonSystem.cpp:147](../../src/PressurePoissonSystem.cpp#L147)):
  integrates the predictor's divergence over each cell by outward-normal
  flux summation (FVM), reusing the same staggered index mapping as `mesh.md §2`:
  $$
  b_P = -\frac{\rho}{\Delta t}\Big[u^*_E S_e - u^*_W S_w + v^*_N S_n - v^*_S S_s\Big]
  $$
  The minus sign matches the positive-definite matrix ($A\,p \approx -V\nabla^2 p$),
  so $A\,p = b$ is the discrete form of $\nabla^2 p = (\rho/\Delta t)\,\nabla\cdot\mathbf{u}^*$.
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
- **Which solver runs:** plain `CG` for the two cavities (serial); `PCG`
  (Jacobi-preconditioned CG, OpenMP-parallel, threads set by
  `OMP_NUM_THREADS`) for the cylinder scenarios, whose stretched meshes make
  the matrix much worse conditioned. Selected in
  [main.cpp:385-386](../../src/main.cpp#L385-L386).
- **Stopping criterion** (both): $\lVert r\rVert_2 < 10^{-6}\lVert b\rVert_2$,
  at most $\max(200, 5\sqrt{n})$ iterations
  ([solver_sparse.cpp:195](../../src/solver/sparse/solver_sparse.cpp#L195),
  [:293](../../src/solver/sparse/solver_sparse.cpp#L293)). If the cap is
  reached a `CG WARNING`/`PCG WARNING` is printed and the time step continues
  with the current iterate; `tools/run_case.sh` filters the PCG warnings out
  of the stored log.
- `GS`, `SOR` and `Jacobi` exist in the same file for experimentation and
  are not used by the time loop.

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
    instead of velocity ones (see [boundary_conditions.md §6](boundary_conditions.md#6-temperature-dhc-only)).
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

- Both the momentum predictor and the energy update are **explicit** in time;
  only the pressure equation is solved implicitly (as a constraint).
- The time step is therefore limited by the explicit diffusion and
  convection limits given in [equations.md §6](equations.md#6-stability-of-the-explicit-scheme),
  $\Delta t \lesssim h_{\min}^2/(8\,\max(\nu,\alpha))$ and
  $\Delta t \lesssim h_{\min}/|\mathbf{u}|_{\max}$.
- Lid-driven cavity and cylinder: $\nu = 1/Re$ is small, `dt ~ 1e-3` is
  stable on the shipped grids (except $Re = 100$ on the 120×120 cavity grid,
  which needs `dt = 3e-4`).
- Differentially heated cavity: $\nu = Pr$ and $\alpha = 1$ are $O(1)$, so
  `dt ~ 1e-5` or smaller is needed — the lid-cavity `dt` diverges to NaN
  (the CG solve fails first). Validated settings:
  [../scenarios/differentially_heated_cavity.md](../scenarios/differentially_heated_cavity.md).

---

## 7. The immersed solid obstacle (Square Cylinder scenario)

How an interior solid obstacle is bolted onto the same core solver used by
the two closed-cavity scenarios, without rewriting the mesh generator or the
per-timestep pipeline in §0.

- **Marking solid cells**: `MeshConfig::mark_solid_block(pMesh, x0, x1, y0, y1)`
  ([src/mesh/MeshGenerator.cpp](../../src/mesh/MeshGenerator.cpp)) runs once,
  right after `build_connectivity(pMesh)` and before the pressure matrix is
  assembled. It snaps the requested bounds to the nearest existing mesh face
  (so the obstacle boundary coincides exactly with P-cell faces), flags every
  P-cell inside as `is_solid`, and re-derives `bc_east/west/north/south` (and
  `n_east/...= -1`) at every fluid/solid interface — **the exact same
  mechanism `build_connectivity` already uses for the outer domain edge**,
  just triggered by an interior condition instead of a grid-extreme index.
  Never called for Lid-Driven Cavity or DHC.
- **Pressure**: `LaplaceOperator` already looped per-cell checking `bc_*`
  flags to decide whether to add a neighbor's contribution — a solid cell
  only needed a 2-line early return (trivial identity row, decoupling its
  pressure value), and fluid cells adjacent to the obstacle automatically get
  a zero-flux face there with no other change. `DivergenceOperator` gets the
  same kind of early return (solid cells contribute no divergence).
- **Momentum (the part that needed real new logic)**: unlike the pressure
  operator, `DiffusionOperator` and `MomentumSolver::uRHS/vRHS` never checked
  `bc_*` — they used raw index-range comparisons against the *global* mesh
  extent to decide between reading a real neighbor or using the wall
  ghost-image formula. Extending this to an *interior* obstacle required a
  shared helper,
  [`src/mesh/ObstacleMask.hpp`](../../src/mesh/ObstacleMask.hpp)
  (`uNodeBlocked`/`vNodeBlocked`/`uNodeOnObstacleSurface`/`vNodeOnObstacleSurface`),
  queried at every neighbor lookup: if a U/V neighbor falls strictly inside
  the solid block, the stencil now falls into the *same* ghost-image branch
  already used for outer walls, just referencing the obstacle's snapped face
  coordinate (`pMesh.obstacle_xmin/xmax/ymin/ymax`) instead of
  `uMesh_.yvc[0]`/`[M]`. For LDC/DHC, `is_solid` is never true anywhere, so
  every one of these new branches is dead code and the stencils are
  byte-identical to before.
- **Obstacle-surface velocity nodes**: `MomentumSolver::applyBoundaryConditions`
  gained a second pass, zeroing any U/V node that sits exactly on the
  fluid/solid interface (no-penetration) — same treatment as the outer-wall
  boundary nodes, just driven by `uNodeOnObstacleSurface`/`vNodeOnObstacleSurface`
  instead of a fixed index.
- **Inlet/outlet**: see
  [boundary_conditions.md](boundary_conditions.md) for how the
  `INLET`/`INLET_PARABOLIC`/`OUTLET` boundary types are handled — including a
  bug fixed along the way: `ghostValue(INLET, faceVal)` previously returned
  the prescribed *normal* inflow speed as a *tangential* ghost value, which
  would have corrupted the cross-stencil at an inlet.
- **Linear solver choice**: the near-obstacle mesh's high cell aspect ratio
  raises the pressure Laplacian's condition number enough that plain CG
  oscillates without converging in the iteration budget that works fine for
  LDC/DHC's milder meshes — this scenario uses Jacobi-preconditioned CG
  (`PCG`) instead, selected in `main.cpp` by scenario. LDC/DHC are
  unaffected (still plain `CG`).
- **Forces**: [`ObstacleForces`](../../src/postprocessing/ObstacleForces.hpp)
  surface-integrates pressure and viscous shear over the obstacle's 4 flat
  faces using the fluid-side cell values (never the solid cell's pinned,
  physically meaningless pressure) into `Cd`/`Cl`.

---

## 8. Possible extension: turbulence modelling (RANS / LES)

The solver currently resolves laminar flow only. A natural extension is an
eddy-viscosity turbulence model, either RANS (Spalart–Allmaras or $k$–$\omega$
SST) or LES (Smagorinsky or dynamic Smagorinsky). Both replace the molecular
viscosity with an effective one, $\nu_{\mathrm{eff}} = \nu + \nu_t$. The
architecture already supports most of this:
- The extra transported scalars of a RANS model ($\tilde\nu$, or $k$ and
  $\omega$) are cell-centred like $\theta$, so they could reuse the
  `EnergySolver` pattern of §5. That means upwind convection, ghost-cell walls
  and AB2, added to the time loop without touching the pressure projection.
- An LES model needs no extra transport equation: $\nu_t = (C_s\Delta)^2|\bar S|$
  is computed from the resolved strain rate at each step.

The main code change is in `DiffusionOperator`, which today assumes a
constant scalar `nu_`. It would need a spatially varying $\nu_t$ interpolated
to the staggered faces. It would also need the full stress divergence
$\nabla\cdot[\nu_{\mathrm{eff}}(\nabla\mathbf{u}+\nabla\mathbf{u}^{\mathsf T})]$,
because the transpose term no longer vanishes when the viscosity varies in
space.

Three numerical constraints follow from the current design (§6):
- The fully explicit scheme ties the time step to
  $\Delta x^2/\nu_{\mathrm{eff}}$, and $\nu_t$ is often orders of magnitude
  larger than $\nu$. This makes implicit diffusion close to a prerequisite.
- The stiff source and destruction terms of $k$–$\omega$ need a point-implicit
  treatment to keep $k,\omega>0$.
- Wall-resolved models require $y^+\lesssim 1$, which the existing `sinh`
  stretching can provide, but the Laplacian is only first order on such grids
  (see [mesh.md](mesh.md)). Spalart–Allmaras and SST also need the wall
  distance field.

Finally, LES is intrinsically three-dimensional: 2D turbulence has an inverse
energy cascade, so a 2D LES is not physically meaningful. The realistic
targets for this 2D code are therefore steady RANS and unsteady RANS (URANS).
A good benchmark is the square cylinder at $Re = 22\,000$ (Lyn et al. 1995),
which extends the existing Square Cylinder scenario. A genuine LES would
require a spanwise direction, i.e. a 3D extension of the mesh and operators.

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
| Obstacle-mask helpers | [src/mesh/ObstacleMask.hpp](../../src/mesh/ObstacleMask.hpp) | header-only |
| Drag/lift/Strouhal | [src/postprocessing/ObstacleForces.hpp](../../src/postprocessing/ObstacleForces.hpp) | [src/postprocessing/ObstacleForces.cpp](../../src/postprocessing/ObstacleForces.cpp) |
| Main driver / time loop | — | [src/main.cpp](../../src/main.cpp) |
| Operator verification (MMS) | — | [src/verification/verify_operators.cpp](../../src/verification/verify_operators.cpp) |
