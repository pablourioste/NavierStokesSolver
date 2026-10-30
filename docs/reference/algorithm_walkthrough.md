# NavierStokesSolver — Algorithm Walkthrough

A 2D incompressible-flow solver using **Chorin's projection method** on a
staggered MAC grid. This document follows the core algorithm through the code
for the **lid-driven cavity**, function by function. The other scenarios (the
Differentially Heated Cavity and the Square Cylinder) reuse the same core; see
[theory/solver.md](../theory/solver.md) and [scenarios/](../scenarios/) for what they
add.

---

## Table of Contents

1. [Overview](#1-overview)
2. [File layout](#2-file-layout)
3. [The staggered MAC grid](#3-the-staggered-mac-grid)
4. [Main data structures](#4-main-data-structures)
5. [Global algorithm: Chorin's projection method](#5-global-algorithm-chorins-projection-method)
6. [Mesh generation — MeshGenerator](#6-mesh-generation--meshgenerator)
7. [Discrete operators — src/operators/](#7-discrete-operators--srcoperators)
8. [Predictor step — MomentumSolver](#8-predictor-step--momentumsolver)
9. [Corrector step — PressurePoissonSystem](#9-corrector-step--pressurepoissonsystem)
10. [Linear solvers — solver_sparse](#10-linear-solvers--solver_sparse)
11. [Main loop — main.cpp](#11-main-loop--maincpp)
12. [Call graph](#12-call-graph)
13. [Operator verification](#13-operator-verification)
14. [Post-processing and validation](#14-post-processing-and-validation)
15. [Numerical parameters and stability](#15-numerical-parameters-and-stability)

---

## 1. Overview

The solver integrates the 2D incompressible Navier–Stokes equations:

```
∂u/∂t + (u·∇)u = -∇p/ρ + ν∇²u
∇·u = 0
```

**Test case:** square cavity with a moving lid (*lid-driven cavity*).
- Domain: L = H = 1.0
- Top (north) lid: prescribed velocity u = U_lid
- Other walls: no-slip (u = v = 0)
- Validation: reference data of Ghia et al. (1982), Re ∈ {100, 400, 1000, …, 10000}

**Key numerical choices:**
- Spatial discretisation: 2nd-order finite differences / finite volumes on a staggered MAC grid
- Time integration: explicit Adams–Bashforth 2 (AB2), with an explicit Euler first step
- Pressure–velocity coupling: Chorin projection (fractional step)
- Pressure linear system: warm-started Conjugate Gradient (CG)

---

## 2. File layout

Only the files used by the core algorithm are listed; the full tree is in the
[README](../../README.md#repository-layout) and in [code_map.md](code_map.md).

```
NavierStokesSolver/
├── CMakeLists.txt                        # CMake build (C++17), three executables
├── run.sh                                # Configure + build + run helper
│
├── input/
│   ├── scenario.txt                      # Selects the active scenario (read first)
│   └── lid_driven_cavity/
│       ├── data.txt                      # Domain geometry and zone corner nodes
│       ├── stretch.txt                   # Cell counts and stretching factors (beta)
│       ├── sim_params.txt                # U_lid, dt, rho, max_steps, ss_tol, solve_energy
│       ├── boundaries.txt                # Boundary conditions (WALL, WALL_FIXED_VALUE, …)
│       ├── reynolds_cases.txt            # Reynolds numbers to sweep
│       └── ghia_reference.txt            # Tabulated data of Ghia et al. (1982)
│
├── output/lid_driven_cavity/             # Simulation results (see README "Output layout")
│   ├── mesh/                             # mesh_P/U/V.vtk (ParaView)
│   └── Re_<N>/
│       ├── fields/                       # pressure.vtk, u_velocity.vtk, v_velocity.vtk
│       └── data/                         # centerline_u.csv, centerline_v.csv
│
├── src/
│   ├── main.cpp                          # Driver: mesh → time loop → VTK/CSV output
│   ├── mainUtils.hpp                     # Interactive mesh inspector (debugging)
│   ├── MomentumSolver.hpp/.cpp           # Predictor step (AB2, diffusion + convection)
│   ├── PressurePoissonSystem.hpp/.cpp    # Corrector step (PPE, projection)
│   │
│   ├── mesh/
│   │   ├── Mesh_config.h                 # GridData, InternalCell, BoundaryFace, MeshConfig
│   │   └── MeshGenerator.cpp             # Generation and connectivity of the P/U/V meshes
│   │
│   ├── operators/                        # Standalone discrete operators
│   │   ├── LaplaceOperator.hpp/.cpp      # FVM Laplacian ∇²p  (P → P)
│   │   ├── GradientOperator.hpp/.cpp     # Discrete gradient ∇p  (P → U,V)
│   │   ├── DivergenceOperator.hpp/.cpp   # Discrete divergence ∇·u  (U,V → P)
│   │   └── DiffusionOperator.hpp/.cpp    # Viscous diffusion ν∇²u  (U,V → U,V)
│   │
│   ├── solver/sparse/
│   │   ├── solver_class_sparse.hpp       # CSR format, Vector alias, solver declarations
│   │   ├── solver_sparse.cpp             # Jacobi, GS, SOR, CG, PCG
│   │   └── main_solver_sparse.cpp        # Standalone solver test program (not built)
│   │
│   └── verification/
│       └── verify_operators.cpp          # 10-suite MMS test (executable: operator_verify)
│
├── plotters/
│   ├── plot_results.py                   # Lid-driven cavity figures vs Ghia et al.
│   └── plot_convergence.py               # Operator convergence-rate figure (MMS)
│
└── docs/
    ├── reference/algorithm_walkthrough.md     # This document
    ├── reference/operators_verification.md    # Detailed operator and test report
    └── legacy/Lid_driven_structured.cpp       # Legacy single-file reference implementation
```

**Executables built in `build/`:**
- `mesh_gen` — main solver (source: `src/main.cpp` + the rest of `src/`)
- `reynolds_sweep` — the same driver under a separate name (`src/reynolds_sweep.cpp` includes `main.cpp`)
- `operator_verify` — standalone operator verification (source: `src/verification/verify_operators.cpp`)

---

## 3. The staggered MAC grid

On a MAC (*Marker-and-Cell*) grid each physical variable lives at a different
place in the cell:

```
      v(i,j+1)
        |
u(i,j)--P(i,j)--u(i+1,j)
        |
      v(i,j)
```

| Variable | Position | Indices | Flat index |
|----------|----------|---------|------------|
| Pressure `p` | Cell centre | i=1..N, j=1..M | `pIdx(i,j) = (j-1)·N + (i-1)` |
| Velocity `u` | Centre of the vertical face | k=1..N+1, j=1..M | `uIdx(k,j) = (j-1)·(N+1) + (k-1)` |
| Velocity `v` | Centre of the horizontal face | i=1..N, l=1..M+1 | `vIdx(i,l) = (l-1)·N + (i-1)` |

This arrangement couples the velocity divergence and the pressure gradient
naturally, **removing the checkerboard oscillations** that appear on
collocated grids.

Boundary nodes hold prescribed values and are **not updated** by the time
advance:
- U at k=1 (West wall) and k=N+1 (East wall)
- V at l=1 (South wall) and l=M+1 (North wall)

---

## 4. Main data structures

Defined in [src/mesh/Mesh_config.h](../../src/mesh/Mesh_config.h).

### `InternalCell` — internal cell of the fluid domain
```cpp
struct InternalCell {
    int global_id;              // Global linear index (0-based)
    double x, y;                // Centroid coordinates
    double Vol;                 // Cell volume
    int n_east, n_west,         // Neighbour cell IDs (-1 at a boundary)
        n_north, n_south;
    bool bc_east, bc_west,      // Does the face touch a boundary?
         bc_north, bc_south;
    bool is_solid = false;      // Immersed obstacle cell (Square Cylinder only)
};
```

### `BoundaryFace` — boundary face
```cpp
struct BoundaryFace {
    double x, y;                // Face centre
    double Area;                // Face area
    int neighbor_id;            // Adjacent internal cell
    BCType type;                // WALL_ADIABATIC, WALL_FIXED_VALUE, INLET, OUTLET, SYMMETRY, WALL_ISOTHERMAL
    double value;               // Prescribed value (velocity, temperature, …)
    std::string boco;           // Name of the boundary-condition group
};
```

### `GridData` — a complete mesh (P, U or V)
```cpp
struct GridData {
    std::vector<double> xvc, yvc;                    // Face positions
    std::vector<double> x, y;                        // Cell centres (ghosts included)
    std::vector<std::vector<double>> Se, Sw, Sn, Ss; // E/W/N/S face areas
    std::vector<std::vector<double>> V;              // Cell volumes
    std::vector<int> id_to_i, id_to_j;               // Global index → (i,j)
    int num_active_nodes;                            // Number of internal cells
    int N_cells_x, M_cells_y;
    std::vector<InternalCell> cells;
    std::vector<BoundaryFace> bound_west, bound_east,
                              bound_south, bound_north;
    // + bound_*_T (temperature BCs) and obstacle_* bounds for the other scenarios
};
```

There are **three `GridData` instances**: `pMesh` (pressure), `uMesh`
(u-velocity), `vMesh` (v-velocity).

---

## 5. Global algorithm: Chorin's projection method

For each time step `n → n+1`:

```
┌─────────────────────────────────────────────────────────┐
│  Step 1 — PREDICTOR (MomentumSolver::advance)            │
│                                                          │
│  u* = u^n + dt · [ν∇²u - (u·∇)u]   (AB2 or Euler)        │
│  v* = v^n + dt · [ν∇²v - (u·∇)v]                         │
│                                                          │
│  The predicted velocity u* does NOT satisfy ∇·u* = 0     │
└───────────────────────┬─────────────────────────────────┘
                        │
                        ▼
┌─────────────────────────────────────────────────────────┐
│  Step 2a — PRESSURE RHS (PressurePoissonSystem)          │
│                                                          │
│  ∇²p = (ρ/dt) · ∇·u*                                     │
└───────────────────────┬─────────────────────────────────┘
                        │
                        ▼
┌─────────────────────────────────────────────────────────┐
│  Step 2b — SOLVE THE PPE (Conjugate Gradient)            │
│                                                          │
│  A · p = b,   A ≈ −V∇²,   b = −(ρ/dt)·∫∇·u* dV           │
└───────────────────────┬─────────────────────────────────┘
                        │
                        ▼
┌─────────────────────────────────────────────────────────┐
│  Step 2c — VELOCITY CORRECTION (projectVelocity)         │
│                                                          │
│  u^{n+1} = u* - (dt/ρ) · ∂p/∂x                           │
│  v^{n+1} = v* - (dt/ρ) · ∂p/∂y                           │
│                                                          │
│  Now ∇·u^{n+1} ≈ 0  (solenoidal)                         │
└─────────────────────────────────────────────────────────┘
```

The key idea is that the method **decouples** the momentum equation from the
pressure field: the velocity is first advanced ignoring pressure, then
corrected so that its divergence is zero.

---

## 6. Mesh generation — MeshGenerator

File: [src/mesh/MeshGenerator.cpp](../../src/mesh/MeshGenerator.cpp)

### `generate_mesh()` — pressure (P) mesh
Builds the base N×M mesh of pressure cells.

1. Uses the domain parameters read from `data.txt` and `stretch.txt` (by `get_data()` / `get_stretch()`).
2. Generates stretched 1D coordinates with `generate_stretched_coords()` in x and y, zone by zone.
3. Computes the face areas `Se, Sw, Sn, Ss` and volumes `V` of every cell.
4. Adds ghost-cell centres at positions 0 and N+1/M+1.

### `generate_u_mesh(pMesh)` — U velocity mesh
Builds the staggered mesh for the horizontal velocity component.
- Size: (N+1) × M nodes
- `xvc_u[k] = pMesh.x[k]` for k=1..N (the U cell limits are consecutive P centres in x)
- The edges k=1 and k=N+1 coincide with the domain walls

### `generate_v_mesh(pMesh)` — V velocity mesh
Same as above for the vertical component.
- Size: N × (M+1) nodes
- `yvc_v[l] = pMesh.y[l]` for l=1..M

### `generate_stretched_coords(vc, num_cells, length, beta, offset)` — stretched coordinates
Generates a 1D node distribution with optional refinement near the walls via a
`sinh` mapping (see [theory/mesh.md §4](../theory/mesh.md#4-mesh-stretching-refinement)).

### `build_connectivity(GridData&)` — neighbour connectivity
For each cell, finds its 4 neighbours (E, W, N, S), raises the `bc_*` flag when
the neighbour is a boundary, and fills `id_to_i`, `id_to_j`.

### `define_boundaries(GridData&, filename)` — boundary conditions
Reads the boundary-condition file, with the format:
```
NAME     FACE    START   END   TYPE              VALUE
lid      NORTH   0       9999  WALL_FIXED_VALUE   1.0
south    SOUTH   0       9999  WALL               0.0
```
and fills the `bound_west`, `bound_east`, `bound_south`, `bound_north` vectors
of `GridData`. The index range is half-open, `[START, END)`.

### `export_mesh_to_VTK(GridData, filename)` — export for ParaView
Writes the mesh as a VTK rectilinear grid: `output/<scenario>/mesh/mesh_P.vtk`,
`mesh_U.vtk`, `mesh_V.vtk`.

---

## 7. Discrete operators — src/operators/

The four operators used by the algorithm live in **standalone classes** in
`src/operators/`. This makes it possible to:
- Reuse them in the main solver and in the verification suite without duplicating code.
- Test them in isolation with manufactured solutions (Section 13).
- Swap an implementation without touching the rest of the solver.

All of them take `const GridData&` references and work with the alias
`Vector = std::vector<double>`.

---

### 7.1 LaplaceOperator — ∇²p (P → P)

File: [src/operators/LaplaceOperator.hpp](../../src/operators/LaplaceOperator.hpp)

Discretises the Laplace operator on the pressure mesh with a 2nd-order FVM.

**Coefficients per cell (i,j):**
```
a_E = Se[i][j] / (x[i+1] - x[i])    a_W = Sw[i][j] / (x[i] - x[i-1])
a_N = Sn[i][j] / (y[j+1] - y[j])    a_S = Ss[i][j] / (y[j] - y[j-1])
a_P = +(a_E + a_W + a_N + a_S)      (positive diagonal)
```

Positive-definite convention: `(A·φ)[k] = Σ a_nb·(φ_P − φ_nb) ≈ −V·∇²φ`,
hence `(A·φ)/V ≈ −∇²φ`.

**Boundary condition:** homogeneous Neumann (∂p/∂n = 0) — boundary faces are
skipped. An `OUTLET` face instead gets a Dirichlet pressure contribution
(Square Cylinder only).

**Key methods:**
- `apply(p)` → `A·p` (matrix-vector product, used by the MMS tests)
- `assemble()` → `CSRMatrix` (full assembly in CSR format)
- `assembleToTriplets()` → `vector<Triplet>` (passed to `pinReference` or to `CSRMatrix::fromTriplets`)
- `pinReference(triplets, b)` → fixes p[0]=0 to remove the all-Neumann singularity

**Used by:** `PressurePoissonSystem::assembleMatrix()` and
`pinReferencePressure()`, `verify_operators.cpp` (Suites 1, 5, 9, 10)

---

### 7.2 GradientOperator — ∇p (P → U, V)

File: [src/operators/GradientOperator.hpp](../../src/operators/GradientOperator.hpp)

Discretises the pressure gradient from the P mesh onto the staggered U/V
meshes with 2nd-order central finite differences.

**Interior U nodes** (k=2..N, j=1..M):
```
(∂p/∂x)|_{k,j} = (p[k,j] - p[k-1,j]) / (x_P[k] - x_P[k-1])
```

**Interior V nodes** (i=1..N, l=2..M):
```
(∂p/∂y)|_{i,l} = (p[i,l] - p[i,l-1]) / (y_P[l] - y_P[l-1])
```

Boundary nodes (k=1, k=N+1, l=1, l=M+1) are set to zero — the wall velocity is
imposed by the boundary condition.

**Purely geometric operator** — it has no dt/ρ factors; those are applied in
`projectVelocity`.

**Method:** `apply(p, grad_u, grad_v)`

**Used by:** `PressurePoissonSystem::projectVelocity()`, `verify_operators.cpp` (Suites 2, 5, 8)

**Algebraic property:** Gᵀ = −D in the unweighted norm (discrete
integration-by-parts identity, checked in Suite 8).

---

### 7.3 DivergenceOperator — ∇·u (U, V → P)

File: [src/operators/DivergenceOperator.hpp](../../src/operators/DivergenceOperator.hpp)

Discretises the velocity divergence by integrating the net flux over each P
cell with the FVM.

**For cell (i,j):**
```
flux  = u_E·Se - u_W·Sw + v_N·Sn - v_S·Ss

where:  u_E = u[uIdx(i+1, j)],   u_W = u[uIdx(i, j)]
        v_N = v[vIdx(i, j+1)],   v_S = v[vIdx(i, j)]

result[pIdx(i,j)] = flux / V[i][j]     (units: 1/time)
```

The result is divided by the volume so it can be compared directly with the
analytical ∇·u.

**Purely geometric operator** — it has no ρ/dt factors.

**Method:** `apply(u, v)` → `Vector` of size `pMesh.num_active_nodes`

**Boundary sampling note:** boundary U/V nodes must be sampled at the **face**
position (`xvc[k-1]`), not at the centroid, for the divergence to be exact
(see Suite 3).

**Used by:** `PressurePoissonSystem::computeRHS()`, `verify_operators.cpp` (Suites 3, 5, 8)

---

### 7.4 DiffusionOperator — ν∇²u (U, V → U, V)

File: [src/operators/DiffusionOperator.hpp](../../src/operators/DiffusionOperator.hpp)

Discretises the viscous terms ν∇²u and ν∇²v on the staggered meshes with
2nd-order finite differences on a non-uniform mesh.

**Formula at an interior U node (k,j):**
```
d²u/dx² = 2/(dx_e + dx_w) · [(u_E - u_P)/dx_e - (u_P - u_W)/dx_w]
d²u/dy² = 2/(dy_n + dy_s) · [(u_N - u_P)/dy_n - (u_P - u_S)/dy_s]
result  = ν · (d²u/dx² + d²u/dy²)
```

**Ghost cells at the walls (image method):**
For nodes next to the N/S walls (u component) and W/E walls (v component):
```
u_ghost  = 2·u_wall − u_P      →  (u_ghost + u_P)/2 = u_wall  ✓
dy_ghost = 2·(y_P − y_wall)    →  correct distance to the boundary face
```
The wall values are passed to the constructor: `u_ghost_north`,
`u_ghost_south`, `v_ghost_west`, `v_ghost_east`.

**Constructor:**
```cpp
DiffusionOperator(pMesh, uMesh, vMesh, nu,
                  u_ghost_north, u_ghost_south,
                  v_ghost_west,  v_ghost_east);
```

**Key methods:**
- `applyU(u, result)` → ν∇²u at every interior U node (k=2..N, j=1..M)
- `applyV(v, result)` → ν∇²v at every interior V node (i=1..N, l=2..M)
- `diffuseU(u, k, j)` / `diffuseV(v, i, l)` → per-node kernels (used by MomentumSolver)

**Used by:** `MomentumSolver` (as the member `diffOp_` in `uRHS`/`vRHS`),
`verify_operators.cpp` (Suites 4, 6, 7)

---

### Operator summary

| Operator | Domain → Range | Discretisation | Class | Consumers |
|----------|----------------|----------------|-------|-----------|
| **L** Laplacian | P → P | 2nd-order FVM | `LaplaceOperator` | `PressurePoissonSystem`, `verify_operators` |
| **G** Gradient | P → (U,V) | 2nd-order FD | `GradientOperator` | `PressurePoissonSystem`, `verify_operators` |
| **D** Divergence | (U,V) → P | 2nd-order FVM | `DivergenceOperator` | `PressurePoissonSystem`, `verify_operators` |
| **D_ν** Diffusion | (U,V) → (U,V) | 2nd-order FD + ghost cells | `DiffusionOperator` | `MomentumSolver`, `verify_operators` |
| **C(u)** Convection | (U,V) → (U,V) | CDS / Upwind | (inline in MomentumSolver) | `MomentumSolver` |

**Algebraic identities** checked to machine precision (Suites 5 and 8 of
`verify_operators`):
```
D·(G·φ) + (L·φ)/V = 0    (consistency of the Laplacian with D and G)
φᵀ·(D·u) + (G·φ)ᵀ·u = 0  (discrete adjoint Gᵀ = −D)
```

---

## 8. Predictor step — MomentumSolver

Files: [src/MomentumSolver.hpp](../../src/MomentumSolver.hpp) / [src/MomentumSolver.cpp](../../src/MomentumSolver.cpp)

Implements Chorin's **predictor step**: advances the velocities with diffusion
and convection, ignoring the pressure gradient.

### Constructor `MomentumSolver(pMesh, uMesh, vMesh, nu, scheme)`

1. Reads the boundary values from `pMesh.bound_*` through `extractGhost()` / `ghostValue()`.
2. Builds the member `DiffusionOperator diffOp_` with the extracted wall values:
   - `u_ghost_north_` (lid or no-slip), `u_ghost_south_` (no-slip)
   - `v_ghost_west_`, `v_ghost_east_` (no-slip)
3. The convection scheme is chosen at construction: `ConvectionScheme::CDS`
   (2nd order, the default) or `ConvectionScheme::UPWIND` (1st order).

---

### `applyBoundaryConditions(u, v)` — impose the boundary conditions

Assigns the prescribed values to the boundary nodes of the velocity vectors:
- U nodes at k=1 and k=N+1: inlet velocity, zero (wall) or zero-gradient (outlet)
- V nodes at l=1 and l=M+1: zero (no penetration)

---

### `uRHS(u, v, k, j)` — right-hand side of the x-momentum equation

Computes `ν∇²u − (u·∇)u` at the interior U node (k,j):

**Diffusion** (delegated to `DiffusionOperator::diffuseU`):
```
diffuseU(u, k, j)  →  ν·(d²u/dx² + d²u/dy²)
```

**Convection (CDS):**
```
u_P · (u_E − u_W)/(dx_e + dx_w)   +   v_interp · (u_N − u_S)/(dy_n + dy_s)
```

**Convection (UPWIND):**
```
u_P≥0 ? (u_P − u_W)/dx_w : (u_E − u_P)/dx_e     +   the analogue for v
```

**Interpolation of v to the U node** (4-point average):
```
v_interp = 0.25·(v[vIdx(k-1,j)] + v[vIdx(k,j)] + v[vIdx(k-1,j+1)] + v[vIdx(k,j+1)])
```

---

### `vRHS(u, v, i, l)` — right-hand side of the y-momentum equation

Symmetric to `uRHS` for interior V nodes (i,l). Diffusion is delegated to
`diffOp_.diffuseV`; u is interpolated to the V node from the 4 surrounding U
nodes.

---

### `advance(u, v, dt)` — Adams–Bashforth 2 time advance

From the second step on, AB2 gives 2nd order in time:

```
RHS^n     = [uRHS(k,j) for every interior node]
RHS^{n-1} = stored from the previous step in rhs_u_old_

u* = u^n + dt · [(3/2)·RHS^n − (1/2)·RHS^{n-1}]   (AB2)
```

The first step uses explicit Euler (`first_step_ == true`). All increments are
computed **before any update**, always from the `u^n` values.

---

## 9. Corrector step — PressurePoissonSystem

Files: [src/PressurePoissonSystem.hpp](../../src/PressurePoissonSystem.hpp) / [src/PressurePoissonSystem.cpp](../../src/PressurePoissonSystem.cpp)

Solves the **pressure Poisson equation (PPE)** and projects the velocity onto
the divergence-free subspace. The class delegates to the operators of
Section 7: `LaplaceOperator` for the matrix, `DivergenceOperator` for the
right-hand side and `GradientOperator` for the correction.

---

### `assembleMatrix()` → `vector<Triplet>`

Returns `LaplaceOperator::assembleToTriplets()`:

```
a_E = Se[i][j] / (x[i+1] - x[i]),   a_W = Sw[i][j] / (x[i] - x[i-1])
a_N = Sn[i][j] / (y[j+1] - y[j]),   a_S = Ss[i][j] / (y[j] - y[j-1])
a_P = +(a_E + a_W + a_N + a_S)
```

Boundary faces → coefficient 0 (homogeneous Neumann). The resulting matrix is
**SPD** → suitable for CG.

---

### `pinReferencePressure(triplets, b)` — fix the reference pressure

Removes the singularity of the all-Neumann system by fixing `p[0] = 0` through
a symmetric substitution of row/column 0. It also sets `b[0] = 0`. Used for
the closed cavities only: the Square Cylinder has a Dirichlet outlet that
already fixes the pressure level.

---

### `setVelocityFields(u_star, v_star)` — store the predicted velocity

Stores the fields `u*` and `v*` for `computeRHS` and `projectVelocity`.

---

### `computeRHS(dt, rho)` → `Vector`

For every pressure cell (i,j):
```
b[k] = −(ρ/dt) · [u*_E·Se - u*_W·Sw + v*_N·Sn - v*_S·Ss]
     = −(ρ/dt) · DivergenceOperator::apply(u*, v*)[k] · V[i][j]
```
The minus sign matches the positive-definite matrix, A·p ≈ −V·∇²p, so
`A·p = b` is the discrete form of ∇²p = (ρ/dt)·∇·u*.

Staggered index mapping:
- East face of P(i,j)  ↔ U node (k=i+1, j)
- West face of P(i,j)  ↔ U node (k=i, j)
- North face of P(i,j) ↔ V node (i, l=j+1)
- South face of P(i,j) ↔ V node (i, l=j)

`b[0] = 0` is applied explicitly in `main.cpp` every step to keep the
reference pressure (closed cavities only).

---

### `projectVelocity(p_sol, dt, rho, u_next, v_next)` — velocity correction

Subtracts the discrete pressure gradient from the predicted velocity:

**Interior U nodes** (k=2..N, j=1..M):
```
u_next[k,j] = u*[k,j] - (dt/ρ) · (p[k,j] - p[k-1,j]) / (x_P[k] - x_P[k-1])
```

**Interior V nodes** (i=1..N, l=2..M):
```
v_next[i,l] = v*[i,l] - (dt/ρ) · (p[i,l] - p[i,l-1]) / (y_P[l] - y_P[l-1])
```

Boundary nodes are not corrected.

---

## 10. Linear solvers — solver_sparse

Files: [src/solver/sparse/solver_class_sparse.hpp](../../src/solver/sparse/solver_class_sparse.hpp) / [solver_sparse.cpp](../../src/solver/sparse/solver_sparse.cpp)

### `Vector` type and CSR format
```cpp
using Vector = std::vector<double>;   // global alias

struct CSRMatrix {
    int n;
    std::vector<double> values;       // non-zero values
    std::vector<int>    col_indices;  // column of each value
    std::vector<int>    row_ptr;      // start of each row (size n+1)
    std::vector<double> diagonal;     // extracted diagonal, O(1) access for the preconditioner
};
```
Built from triplets with `CSRMatrix::fromTriplets(n, triplets)`.

---

### Available solvers

| Class | Method | Notes |
|-------|--------|-------|
| `Jacobi` | Jacobi iteration | Slow, parallelisable |
| `GS` | Gauss–Seidel | Converges ~2× faster than Jacobi |
| `SOR` | GS + relaxation (default ω=1.5) | Faster than plain GS |
| `CG` | Conjugate Gradient | **Default solver** for the SPD PPE (closed cavities) |
| `PCG` | CG with Jacobi (diagonal) preconditioner | Used for the Square Cylinder |

**Warm-started CG:**
```cpp
Vector solve(const CSRMatrix& A, const Vector& b, const Vector& x0);
```
The previous step's solution (`p^n`) is the initial guess, which cuts the
iteration count to ~5–20 in quasi-steady flow.

---

## 11. Main loop — main.cpp

File: [src/main.cpp](../../src/main.cpp)

### Parameters read from the input files

All files are read from `input/<scenario>/`, where `<scenario>` comes from
`input/scenario.txt`.

| File | Parameters |
|------|------------|
| `data.txt` | `L_domain`, `H_domain`, `W_depth`, zone corner nodes |
| `stretch.txt` | Cell counts and stretching factors in x and y |
| `sim_params.txt` | `U_lid`, `dt`, `rho`, `max_steps`, `ss_tol`, `solve_energy`, `Pr` |
| `reynolds_cases.txt` | Reynolds numbers to sweep (one case per value) |

The kinematic viscosity is derived per case: `nu = U_lid · L_domain / Re`.

### Initialisation (runs once)

```
1. read_scenario_name("input/scenario.txt")  → INPUT_DIR = input/<scenario>/
2. MeshConfig::get_data(), get_stretch(), get_sim_params()
3. generate_mesh()          → pMesh
   generate_u_mesh(pMesh)   → uMesh
   generate_v_mesh(pMesh)   → vMesh
4. define_boundaries(pMesh, INPUT_DIR + "boundaries.txt")
5. build_connectivity(pMesh), build_connectivity(uMesh), build_connectivity(vMesh)
6. export_mesh_to_VTK()  ×3  → output/<scenario>/mesh/mesh_P/U/V.vtk

7. PressurePoissonSystem::assembleMatrix()   → triplets
   pinReferencePressure(triplets, dummy_b)
   CSRMatrix::fromTriplets()                 → A (CSR, assembled once)
```

### Per case (each Reynolds number)

```
8.  MomentumSolver momentum(pMesh, uMesh, vMesh, nu)
9.  u = v = p = 0
10. momentum.applyBoundaryConditions(u, v)
```

### Time loop

```
For step = 0, 1, 2, ... up to max_steps or convergence:

  A) PREDICTOR
     u_old = u_star
     momentum.advance(u_star, v_star, dt)           → u*, v* (AB2)
     momentum.applyBoundaryConditions(u_star, v_star)

  B) PRESSURE RHS
     poisson.setVelocityFields(u_star, v_star)
     b_rhs = poisson.computeRHS(dt, rho)            → b = −(ρ/dt)·∫∇·u* dV
     b_rhs[0] = 0                                   → pressure reference

  C) SOLVE THE PPE
     p_sol = CG_solver.solve(A, b_rhs, p_sol)       → warm start

  D) VELOCITY CORRECTION
     poisson.projectVelocity(p_sol, dt, rho, u_star, v_star)
                                                    → u^{n+1} = u* - (dt/ρ)·∇p

  E) RE-IMPOSE THE BOUNDARY CONDITIONS
     momentum.applyBoundaryConditions(u_star, v_star)

  F) STEADY-STATE CHECK (every 100 steps)
     if ||u^{n+1} - u^n||∞ / dt < ss_tol → CONVERGED, exit
     (a non-finite value aborts the case and nothing is written)
```

### Shipped lid-driven cavity parameters

| Parameter | Value | Description |
|-----------|-------|-------------|
| `Re` | from `reynolds_cases.txt` | Reynolds number |
| `U_lid` | 1.0 | Lid velocity |
| `nu` | U·L/Re | Kinematic viscosity |
| `dt` | 1e-3 | Time step |
| `rho` | 1.0 | Density |
| `max_steps` | 200 000 | Maximum number of time steps |
| `ss_tol` | 0.0 | Steady-state tolerance (0 → always run to `max_steps`) |
| Grid | 120 × 120 | Stretched (`stretch.txt`) |

---

## 12. Call graph

```
main()
│
├── Mesh
│   ├── get_data() / get_stretch() / get_sim_params()
│   ├── generate_mesh() → pMesh
│   │   └── generate_stretched_coords() per zone, in x and y
│   ├── generate_u_mesh(pMesh) → uMesh
│   ├── generate_v_mesh(pMesh) → vMesh
│   ├── define_boundaries(pMesh, ...)
│   ├── build_connectivity(pMesh/uMesh/vMesh)
│   └── export_mesh_to_VTK() ×3
│
├── Pressure linear system (once)
│   ├── PressurePoissonSystem::assembleMatrix()   → triplets
│   │   └── LaplaceOperator::assembleToTriplets()
│   ├── pinReferencePressure(triplets, dummy_b)
│   └── CSRMatrix::fromTriplets()                → A (CSR)
│
└── CASE LOOP (one per Re) ────────────────────────────────────
    │
    ├── MomentumSolver::MomentumSolver(pMesh, uMesh, vMesh, nu)
    │   └── DiffusionOperator::DiffusionOperator(...)  ← member diffOp_
    │
    ├── MomentumSolver::applyBoundaryConditions(u, v)   ← init
    │
    └── TIME LOOP ─────────────────────────────────────────────
        │
        ├── MomentumSolver::advance(u, v, dt)
        │   ├── [j=1..M, k=2..N]  uRHS(u, v, k, j)
        │   │   ├── diffOp_.diffuseU(u, k, j)        ← DiffusionOperator
        │   │   │   ├── d²u/dx²  (2nd-order FD)
        │   │   │   └── d²u/dy²  (ghost cell at j=1 and j=M)
        │   │   └── CDS/Upwind convection + 4-point v interpolation
        │   ├── [l=2..M, i=1..N]  vRHS(u, v, i, l)
        │   │   ├── diffOp_.diffuseV(v, i, l)        ← DiffusionOperator
        │   │   └── Convection + 4-point u interpolation
        │   └── u += dt·(3/2·RHS^n − 1/2·RHS^{n-1});  v += … (AB2)
        │
        ├── PressurePoissonSystem::setVelocityFields(u*, v*)
        │
        ├── PressurePoissonSystem::computeRHS(dt, rho)
        │   └── −(ρ/dt) · DivergenceOperator::apply(u*, v*) · V
        │
        ├── b_rhs[0] = 0
        │
        ├── CG::solve(A, b_rhs, p_previous)          ← warm start
        │   ├── r = b - A·p_prev
        │   └── [iter] α, x, r, β, p  (converged when ||r|| < 1e-6·||b||)
        │
        ├── PressurePoissonSystem::projectVelocity(p_sol, dt, rho, u*, v*)
        │   └── GradientOperator::apply(p) → u* -= (dt/ρ)·∂p/∂x,  v* -= (dt/ρ)·∂p/∂y
        │
        ├── MomentumSolver::applyBoundaryConditions(u*, v*)
        │
        └── [every 100 steps]  ||u*−u_old||∞/dt < ss_tol → BREAK
```

---

## 13. Operator verification

Executable: `./build/operator_verify`  
Source: [src/verification/verify_operators.cpp](../../src/verification/verify_operators.cpp)  
Detailed documentation: [operators_verification.md](operators_verification.md)

The verification program uses the **Method of Manufactured Solutions (MMS)**:
a smooth analytical function is chosen, sampled on the mesh, the discrete
operator is applied and the result is compared with the analytical one. This is
repeated on uniform meshes N ∈ {8, 16, 32, 64} and the convergence rate is
estimated as `r = log₂(e_N / e_{2N})`.

When it finishes it writes `output/verification/data/convergence.csv` with all
the errors per suite and N.

### Implemented suites

| Suite | Operator | Manufactured function | Expected result |
|-------|----------|-----------------------|-----------------|
| **1** | `LaplaceOperator` | φ=cos(πx)cos(πy), Neumann | rate ~2.0 |
| **2** | `GradientOperator` | φ=cos(πx)cos(πy), interior nodes | rate ~2.0 (∂p/∂x and ∂p/∂y) |
| **3** | `DivergenceOperator` | u=sin(2πx)cos(πy), face sampling | rate ~2.0 |
| **4** | `DiffusionOperator::applyU` | u=sin(πy) (const. in x), ν=1 | rate ~2.0 |
| **5** | Identity D·G + L/V = 0 | φ=cos(πx)cos(πy) | error ~1e-14..1e-12 |
| **6** | `DiffusionOperator::applyV` | v=sin(πx) (const. in y), ν=1 | rate ~2.0 |
| **7** | Diffusion with a non-zero BC | u=y (linear), u_N=1, u_S=0 | error exactly 0 |
| **8** | Adjoint Gᵀ = −D | random interior vectors | residual ~1e-11 (machine precision) |
| **9** | LaplaceOperator SPD | random vectors | symmetry ~1e-14, φᵀAφ > 0 |
| **10** | Laplacian on a tanh mesh (γ=2) | φ=cos(πx)cos(πy) | rate ~1.0 (*) |

(*) Suite 10 exposes a **limitation of the current stencil**: on non-uniform
meshes (`h_e ≠ h_w`), the cross term `φ'·(h_e−h_w)/dx_i = O(1/N)` reduces the
accuracy to first order. The fix is to replace `1/dx_i` with `2/(h_e+h_w)` in
the FVM coefficients.

### Convergence-rate figure

```bash
./build/operator_verify                # writes output/verification/data/convergence.csv
python3 plotters/plot_convergence.py   # writes output/verification/figures/convergence_rates.{pdf,png}
```

The figure has 3 panels:
- (a) Uniform mesh — every operator, slope ~2.0
- (b) Mesh stretching — uniform vs stretched Laplacian, slope 2.0 vs 1.0
- (c) Discrete identities — residuals of Suites 5, 8 and 9 at machine precision

---

## 14. Post-processing and validation

### Exporting fields to VTK

At the end of each case, `main.cpp` writes three files with the `writeFieldVTK`
lambda:

```
output/<scenario>/Re_<N>/fields/pressure.vtk     — pressure field (P mesh)
output/<scenario>/Re_<N>/fields/u_velocity.vtk   — u component (U mesh)
output/<scenario>/Re_<N>/fields/v_velocity.vtk   — v component (V mesh)
```

Format: ASCII VTK rectilinear grid, viewable in **ParaView**.

### Centreline profiles

```
output/lid_driven_cavity/Re_<N>/data/centerline_u.csv  — u(y) at x = L/2 (U column closest to the centre)
output/lid_driven_cavity/Re_<N>/data/centerline_v.csv  — v(x) at y = H/2 (V row closest to the centre)
```

Columns: `y, y_norm, u, u_norm` and `x, x_norm, v, v_norm`.

### Validation figures

```bash
python3 plotters/plot_results.py
```

For every case in `output/lid_driven_cavity/` the script draws the flow fields
(stream function, vorticity, pressure) and the centreline profiles against
Ghia et al. (1982), then the cross-case figures (profile comparison,
stream-function overview) and `validation_summary.csv`. The Ghia column is
selected from the case's Re (available: 100, 400, 1000, 3200, 5000, 7500,
10000). Results are summarised in
[lid_driven_cavity.md §8](../scenarios/lid_driven_cavity.md#8-validation-results).

---

## 15. Numerical parameters and stability

### Time-step limit (explicit diffusion)
```
dt ≤ C · dx² / ν
```
The explicit scheme also has a convective limit (CFL, `dt ≤ C·dx/|u|`); on
stretched grids the smallest cell sets both limits. See
[theory/solver.md §6](../theory/solver.md#6-stability-this-is-a-fully-explicit-scheme)
and the heated-cavity caveat in
[scenarios/differentially_heated_cavity.md](../scenarios/differentially_heated_cavity.md).

### Linear system convergence (PPE)
- Relative tolerance: `||r|| < 1e-6 · ||b||`
- Warm-started CG: ~5–20 iterations per step in quasi-steady flow
- Condition number: κ(A) ~ O(N²)

### Convergence to steady state
- Criterion: `||u^{n+1} - u^n||∞ / dt < ss_tol`, checked every 100 steps
- With `ss_tol = 0` (the shipped value) the run always goes to `max_steps`

### AB2 vs Euler time integration
- The first step uses Euler (no history yet) — a small initial transient
- From the second step on, AB2 gives 2nd order in time
- `rhs_u_old_` / `rhs_v_old_` in `MomentumSolver` hold the previous step's RHS
