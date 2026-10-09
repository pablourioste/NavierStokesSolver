# Differentially Heated Cavity (DHC) — Implementation Documentation

This guide explains, end to end, how the **Differentially Heated Cavity** (DHC)
is implemented in the solver, which equations it solves, and **how information
flows** from the input files → the mesh and its discretization (including
*stretching*) → the operators → the momentum, pressure and energy solvers.

The code fully reuses the *Lid-Driven Cavity* solver. The DHC is a **mode**
activated by a single flag, adding two modular pieces: the energy equation
([EnergySolver.cpp](../../src/EnergySolver.cpp)) and the Boussinesq buoyancy term
([Buoyancy.hpp](../../src/Buoyancy.hpp)). **No core file (momentum, pressure,
operators) was modified** to support the DHC.

---

## 1. The physical equations

The DHC solves the incompressible Navier–Stokes equations **coupled with the
energy equation** under the **Boussinesq approximation**, in the **Ra–Pr**
non-dimensional scaling (de Vahl Davis, 1983).

### 1.1 Mass conservation (continuity)

$$
\nabla \cdot \mathbf{u} = 0
$$

Enforced exactly via the **projection method** (pressure step).

### 1.2 Momentum (with Boussinesq buoyancy)

$$
\frac{\partial \mathbf{u}}{\partial t} + (\mathbf{u}\cdot\nabla)\mathbf{u}
= -\nabla p + \nu\,\nabla^2 \mathbf{u} + \mathbf{f}_b
$$

In the Ra–Pr scaling, the **momentum diffusion** equals `nu = Pr` and the
buoyancy force acts only in the vertical direction:

$$
\mathbf{f}_b = (0,\; Ra\,Pr\,\theta)
$$

where `θ` is the non-dimensional temperature. This appears directly in the code:

- In [main.cpp:196](../../src/main.cpp#L196): `const double nu = energy_on ? Pr : ...`
- In [Buoyancy.hpp:34](../../src/Buoyancy.hpp#L34): `const double coeff = dt * Ra * Pr;`

### 1.3 Energy (temperature advection–diffusion)

$$
\frac{\partial \theta}{\partial t} + (\mathbf{u}\cdot\nabla)\theta
= \alpha\,\nabla^2 \theta
$$

With thermal diffusivity `α = 1` in the Ra–Pr scaling
([main.cpp:221](../../src/main.cpp#L221), `alpha=1.0` is passed to `EnergySolver`).
The equation is stated literally in the header of
[EnergySolver.hpp:14](../../src/EnergySolver.hpp#L14):

```
d(theta)/dt + (u . grad) theta = alpha * laplacian(theta)
```

### 1.4 Non-dimensional numbers

| Symbol | Meaning | Where it enters |
|--------|---------|-----------------|
| `Ra` (Rayleigh) | strength of natural convection | swept in `rayleigh_cases.txt` → `coeff = dt·Ra·Pr` |
| `Pr` (Prandtl)  | momentum / thermal diffusivity ratio | `sim_params.txt` → `nu = Pr` |
| `α = 1`         | non-dimensional thermal diffusivity | fixed in `main.cpp` |

---

## 2. Problem selection: how the DHC is activated

The switch between **Lid-Driven Cavity** and **Differentially Heated Cavity** is
a single integer in [input/sim_params.txt](../../input/sim_params.txt):

```text
solve_energy = 1     # 0 -> Lid cavity (momentum + pressure)
                     # 1 -> Differentially Heated Cavity (+ energy + buoyancy)
Pr = 0.71            # Prandtl number (used only if solve_energy = 1)
```

Flag flow:

1. [`MeshConfig::get_sim_params`](../../src/mesh/MeshGenerator.cpp#L686) reads the
   file and `update_parameter` assigns `solve_energy` and `Pr` to the class
   members ([MeshGenerator.cpp:732](../../src/mesh/MeshGenerator.cpp#L732)).
2. In [main.cpp:129](../../src/main.cpp#L129):
   `const bool energy_on = (myConfig.solve_energy != 0);`
3. If `energy_on`:
   - The Rayleigh numbers to sweep are read from `rayleigh_cases.txt`
     ([main.cpp:135](../../src/main.cpp#L135)).
   - **Temperature** boundary conditions are loaded from `boundaries_T.txt`
     ([main.cpp:136](../../src/main.cpp#L136)).
   - The case loop iterates over **Rayleigh** instead of Reynolds.

---

## 3. The input files and what each controls

Everything in `input/`. `main()` reads them at startup in this order
([main.cpp:82-93](../../src/main.cpp#L82-L93)):

| File | Read by | Controls |
|------|---------|----------|
| [data.txt](../../input/data.txt) | `get_data` | Geometry: `L_domain`, `H_domain`, depth `W_depth`, and the NW/NE/SW/SE nodes that split the domain into zones |
| [stretch.txt](../../input/stretch.txt) | `get_stretch` | Cell count per zone and the **stretching β** (refinement) |
| [sim_params.txt](../../input/sim_params.txt) | `get_sim_params` | `dt`, `rho`, `U_lid`, `nu`, `max_steps`, `ss_tol`, **`solve_energy`**, **`Pr`** |
| [boundaries.txt](../../input/boundaries.txt) | `define_boundaries` | **Velocity** boundary conditions |
| [boundaries_T.txt](../../input/boundaries_T.txt) | `define_temperature_boundaries` | **Temperature** boundary conditions (DHC only) |
| [rayleigh_cases.txt](../../input/rayleigh_cases.txt) | `read_rayleigh_cases` | List of Rayleigh numbers to sweep |

### 3.1 `data.txt` — geometry

```text
L_domain = 1.0    # domain width
H_domain = 1.0    # domain height
W_depth  = 0.1    # depth (truly 2D, but areas are multiplied by W_depth)
# NW/NE/SW/SE nodes = 1.0  -> single-zone domain (no obstacle)
```

For the standard square cavity, the four nodes are placed so that the
**"west/south" zone covers the whole domain** and the "center"/"east" zones have
zero length (with `N_cells = 0` they are not copied into the mesh).

### 3.2 `stretch.txt` — discretization and refinement

```text
beta_x_west   = 0.8     # stretching in X
beta_y_south  = 0.8     # stretching in Y
beta_x_center = 0.0     # 0.0 -> uniform mesh in that zone
N_cells_west  = 50      # 50 cells in X
M_cells_west  = 50      # 50 cells in Y
N_cells_center = 0      # empty zones
```

The total number of cells is the sum of the three zones:
`N_cells = N_cells_west + N_cells_center + N_cells_east`
([MeshGenerator.cpp:158](../../src/mesh/MeshGenerator.cpp#L158)).

### 3.3 `boundaries_T.txt` — thermal boundaries (the DHC's key)

```text
# BOCO_NAME  FACE   IDX_START IDX_END  TYPE             VALUE
hot_west     WEST   0   9999   WALL_ISOTHERMAL   1.0    # hot wall  θ = 1
cold_east    EAST   0   9999   WALL_ISOTHERMAL   0.0    # cold wall θ = 0
top          NORTH  0   9999   WALL_ADIABATIC    0.0    # zero flux
bottom       SOUTH  0   9999   WALL_ADIABATIC    0.0    # zero flux
```

This is the physical configuration that **defines** the differentially heated
cavity: hot west wall, cold east wall, adiabatic top and bottom. The horizontal
temperature gradient generates the buoyancy that drives the fluid.

> Note: the `IDX_START IDX_END` range is **half-open** `[start, end)` (the loop
> uses `k < idx_end`), which is why `0 9999` covers the entire wall.

---

## 4. The mesh: generation, staggered grid and stretching

### 4.1 Staggered grid (MAC)

**Three meshes** are generated from the pressure mesh
([main.cpp:86-88](../../src/main.cpp#L86-L88)):

```
pMesh = generate_mesh();           // pressure and temperature θ (cell center)
uMesh = generate_u_mesh(pMesh);    // velocity u  (vertical faces,   N+1 in X)
vMesh = generate_v_mesh(pMesh);    // velocity v  (horizontal faces, M+1 in Y)
```

- **p and θ** live at the **center** of each cell.
- **u** lives on the vertical faces (the U mesh has `N+1` columns and shares P's
  Y faces, [MeshGenerator.cpp:343](../../src/mesh/MeshGenerator.cpp#L343)).
- **v** lives on the horizontal faces (the V mesh has `M+1` rows and shares P's
  X faces, [MeshGenerator.cpp:400](../../src/mesh/MeshGenerator.cpp#L400)).

This layout is why `θ` is stored just like the pressure: the staggered
velocities **are directly** the face-normal velocities of each P cell, with no
interpolation needed for the convective flux (see §6.2).

### 4.2 Stretching (refinement)

Refinement is done with a `sinh` function in
[`generate_stretched_coords`](../../src/mesh/MeshGenerator.cpp#L79). Idea:

1. If `|β| < 1e-6` → **uniform mesh** (`Δx = L/num_cells`).
2. If `β > 0` → refine **at the end** of the segment (`x = L`).
3. If `β < 0` → refine **at the start** (`x = 0`).

The algorithm:

```
xi      = -1 + 2*i/num_cells          # uniform logical coordinate [-1,1]
g_xi    = (sinh(β(xi - xi_focus)) - val_min) / val_range * 2 - 1   # stretched
x       = (L/2)*(1 + g_xi) + offset   # mapping to physical coordinate
```

The larger `|β|`, the more cells concentrate near the focus. For the DHC this
allows **refining near the hot/cold walls**, where thermal boundary layers are
thin and where the Nusselt number is measured.

The `offset` positions each zone (west/center/east) within the global domain, and
then the segments are **merged** into `data.xvc` / `data.yvc`
([MeshGenerator.cpp:208-237](../../src/mesh/MeshGenerator.cpp#L208-L237)).

### 4.3 From coordinates to the geometry the solver uses

After generating the faces (`xvc`, `yvc`), `generate_mesh` computes and stores in
`GridData` everything the operators need
([MeshGenerator.cpp:239-261](../../src/mesh/MeshGenerator.cpp#L239-L261)):

- `x[i]`, `y[j]` — cell **centroids** (includes ghost cells 0 and N+1).
- `Se, Sw, Sn, Ss` — **face areas** (= `Δy·W_depth` or `Δx·W_depth`).
- `V` — cell **volume** (`Δx·Δy·W_depth`).

These structures are **the contract** between the mesh and the solvers: the
`EnergySolver` reads `pMesh.Se/Sw/Sn/Ss/V` and the centroids directly to build
its fluxes and Laplacian (see §6).

### 4.4 Connectivity and temperature boundaries

- [`build_connectivity`](../../src/mesh/MeshGenerator.cpp#L448) builds the neighbor
  table (E/W/N/S) and flags which cells touch a boundary.
- [`define_temperature_boundaries`](../../src/mesh/MeshGenerator.cpp#L993) creates
  the parallel `bound_*_T` lists (west/east/north/south). It **copies the
  geometry** from the velocity lists `bound_*` and then overwrites only the
  `type` and `value` with what `boundaries_T.txt` says. Faces not mentioned
  default to **adiabatic**. This leaves the velocity BCs untouched.

---

## 5. The time loop: how the solver is fed step by step

The heart is in [main.cpp:225-258](../../src/main.cpp#L225-L258). Each `caseVal` of
the sweep is a Rayleigh number; inside, the time loop is:

```cpp
for (int step = 0; step < max_steps; ++step) {
    u_old = u_star;

    // (1) Momentum PREDICTOR (AB2): u* = u^n + dt(3/2 R^n - 1/2 R^{n-1})
    momentum.advance(u_star, v_star, dt);

    // (2) BUOYANCY: injects dt·Ra·Pr·θ into v* (DHC only)
    if (energy_on)
        addBuoyancy(T, v_star, dt, Ra, Pr, pMesh, vMesh);

    // (3) PRESSURE: builds the Poisson RHS from the divergence of u*
    poisson.setVelocityFields(u_star, v_star);
    b_rhs = poisson.computeRHS(dt, rho);
    b_rhs[0] = 0.0;                       // fixes the reference pressure (pinning)

    // (4) Solves the Poisson system with Conjugate Gradient (warm-start)
    p_sol = CG_solver.solve(A, b_rhs, p_sol);

    // (5) PROJECTION: corrects u* -> divergence-free velocity field
    poisson.projectVelocity(p_sol, dt, rho, u_star, v_star);
    momentum.applyBoundaryConditions(u_star, v_star);

    // (6) ENERGY: transports θ with the already corrected field (DHC only)
    if (energy_on)
        energy.advance(T, u_star, v_star, dt);

    // (7) Steady-state check every 100 steps
}
```

**Order matters.** Buoyancy is injected into `v*` **after** the momentum
predictor and **before** the pressure projection: this way the projection makes
the final field incompressible **including** the buoyant forcing. And the energy
is advanced **last**, using the already corrected, divergence-free velocity field
(`u_star`, `v_star`), which keeps θ transport consistent. This is exactly why the
`MomentumSolver` and `PressurePoissonSystem` **did not need to be touched**: the
DHC couples in via external injection.

### 5.1 Initial condition

- `u_star, v_star, p_sol = 0` ([main.cpp:211-213](../../src/main.cpp#L211-L213)).
- `T = 0.5` in all cells (mid-value between the hot wall θ=1 and the cold θ=0)
  ([main.cpp:223](../../src/main.cpp#L223)).

---

## 6. The EnergySolver in detail

File: [src/EnergySolver.cpp](../../src/EnergySolver.cpp). Solves
`dθ/dt + (u·∇)θ = α∇²θ` on the **pressure mesh**. Its structure mirrors the
`MomentumSolver`.

### 6.1 Construction: reading the temperature BCs

In the constructor ([EnergySolver.cpp:17-34](../../src/EnergySolver.cpp#L17-L34)),
`extractWallT` reads the type and value of each wall from `pMesh.bound_*_T`
(`bc_west_`, `Tw_west_`, …). By default, an empty list → adiabatic.

### 6.2 The convective term (conservative finite volumes)

In [`tRHS`](../../src/EnergySolver.cpp#L48), the face-normal velocities of the P
cell are **directly** the staggered velocities (no interpolation):

```cpp
const double u_e = u[uIdx(i+1, j)];  // east face
const double u_w = u[uIdx(i,   j)];  // west face
const double v_n = v[vIdx(i, j+1)];  // north face
const double v_s = v[vIdx(i, j  )];  // south face
```

Fluxes use first-order **UPWIND** by default (recommended for high-Rayleigh
stability) or **CDS** (central differences). With UPWIND, the value at each face
is taken from the *upstream* cell
([EnergySolver.cpp:100-106](../../src/EnergySolver.cpp#L100-L106)):

```cpp
Fe = u_e * ((u_e >= 0.0) ? T_P : T_E) * Se;   // and similarly Fw, Fn, Fs
conv = (Fe - Fw + Fn - Fs) / Vol;             // divergence per unit volume
```

Since the velocity **normal to the walls is zero** (no penetration), the fluxes
on boundary faces vanish and convection needs no temperature ghost.

### 6.3 The diffusive term (Laplacian with ghost-cell image)

Second-order finite differences on the non-uniform mesh. At the walls the
**image-cell method** is used
([EnergySolver.cpp:59-78](../../src/EnergySolver.cpp#L59-L78)):

- **Isothermal (Dirichlet):** `θ_ghost = 2·θ_wall − θ_P`
- **Adiabatic (Neumann):** `θ_ghost = θ_P` (zero normal gradient)

```cpp
d2T_dx2 = 2/(dx_e+dx_w) * ((T_E - T_P)/dx_e - (T_P - T_W)/dx_w);
d2T_dy2 = 2/(dy_n+dy_s) * ((T_N - T_P)/dy_n - (T_P - T_S)/dy_s);
diff    = alpha * (d2T_dx2 + d2T_dy2);
```

The total RHS is `return diff - conv;`.

### 6.4 Time integration (Adams-Bashforth 2)

In [`advance`](../../src/EnergySolver.cpp#L118): the **entire** RHS is computed from
`θ^n` into a buffer (deferred-Euler style, so updating one cell does not
contaminate its neighbors within the step), and then:

- **First step:** Euler — `θ += dt·rhs`.
- **Rest:** AB2 — `θ += dt·(1.5·rhs − 0.5·rhs_old)`.

> **Note — why does temperature use the same weights as momentum, but pressure
> has no weights?**
>
> The time integration of `θ` is **identical, line for line**, to the
> `MomentumSolver`: same AB2 weights (`3/2`, `−1/2`), same Euler start and the
> same deferred buffer. The `EnergySolver` was deliberately designed as a mirror
> of the `MomentumSolver` for the scalar `θ`. Compare
> [EnergySolver.cpp:128-134](../../src/EnergySolver.cpp#L128-L134) with
> [MomentumSolver.cpp:289-300](../../src/MomentumSolver.cpp#L289-L300).
>
> **Pressure is a different case**: it is not an evolution equation (no
> `∂p/∂t`), but an **elliptic constraint**. At each step the Poisson equation
> `∇²p = (ρ/dt)·∇·u*` is solved with Conjugate Gradient
> ([main.cpp:236](../../src/main.cpp#L236)) to enforce `∇·u = 0`. That is why it has
> **no AB2 weights and no `rhs_old`**: it is an implicit *solve* yielding the
> instantaneous pressure, not a time march.
>
> | Field | Equation type | Time advance |
> |-------|---------------|--------------|
> | `u, v` (momentum) | evolution (parabolic) | explicit AB2 `3/2, −1/2` |
> | `θ` (energy)      | evolution (parabolic) | explicit AB2 `3/2, −1/2` ← same as momentum |
> | `p` (pressure)    | constraint (elliptic) | not integrated; Poisson via CG each step |

---

## 7. The buoyancy coupling (Buoyancy.hpp)

[`addBuoyancy`](../../src/Buoyancy.hpp#L24) is header-only and applied as an
explicit external increment on `v*`:

```cpp
const double coeff = dt * Ra * Pr;
for (int l = 2; l <= M; ++l)          // interior V nodes only
    for (int i = 1; i <= N; ++i) {
        // θ lives in P cells; V node (i,l) sits between P(i,l-1) and P(i,l)
        const double theta_v = 0.5 * (T[pIdx(i, l-1)] + T[pIdx(i, l)]);
        v[vIdx(i, l)] += coeff * theta_v;
    }
```

Important details:

- Since `θ` is on P and `v` is on the horizontal faces, **θ is linearly
  interpolated** to the V node (average of the two adjacent P cells).
- Only **interior** V nodes (`l = 2..M`) are touched; wall nodes (no
  penetration) are left untouched.
- It only adds a **vertical** force (the `v` column), consistent with
  `f_b = (0, Ra·Pr·θ)`.

---

## 8. Validation and output

When `energy_on`, at the end of each case
([main.cpp:406-461](../../src/main.cpp#L406-L461)) the **de Vahl Davis (1983)**
benchmark metrics are computed:

- `u_max` on the vertical mid-plane (`x = L/2`).
- `v_max` on the horizontal mid-plane (`y = H/2`).
- **Average Nusselt number** on the hot (west) wall:

$$
\overline{Nu} = \frac{1}{H}\sum_j \frac{\theta_{wall} - \theta_P}{x_P - x_{wall}}\,\Delta y_j
$$

```cpp
Nu_avg += (Twall - T[k]) / dist * dy;   // local flux * cell height
Nu_avg /= H_domain;
```

Outputs per case (in `output/Ra_<value>/`):

- `temperature.vtk`, `u_velocity.vtk`, `v_velocity.vtk`, `pressure.vtk` (fields in
  rectilinear VTK format, [main.cpp:266-270](../../src/main.cpp#L266-L270)).
- `dhc_metrics.csv` with `Ra, Pr, u_max, v_max, Nu_avg`.

---

## 9. CRITICAL stability warning

Both energy and momentum are **fully explicit**, so the DHC needs a **much
smaller** `dt` than the lid cavity. The explicit diffusion limit is:

$$
dt \lesssim \frac{\Delta x^2}{2\,\max(\alpha,\nu)}, \qquad \alpha = 1
$$

On a production 150² mesh this means `dt ~ 1e-5` or less. Using `dt = 0.001` (the
lid cavity default) **diverges to NaN** (CG fails to converge first). The current
`sim_params.txt` is already set to `dt = 0.00001` for the DHC.

> Verified at Ra=1e4 on a 40² mesh with `dt=2e-5`: `u_max≈16.7`, `v_max≈20.1`,
> `Nu≈2.31` — matches de Vahl Davis (1983) within a few percent.

---

## 10. Information flow summary

```
data.txt ─────────┐
stretch.txt ──────┤→ MeshConfig ──→ generate_mesh ──→ pMesh (xvc,yvc,x,y,Se..Ss,V)
sim_params.txt ───┘                      │                 ├─→ generate_u_mesh → uMesh
   (solve_energy, Pr, dt, rho)           │                 └─→ generate_v_mesh → vMesh
                                         │
boundaries.txt ───→ define_boundaries ──→ bound_*  (velocity BC, ghosts)
boundaries_T.txt ─→ define_temperature_boundaries → bound_*_T (temperature BC)
rayleigh_cases.txt → read_rayleigh_cases → list of Ra

                         ┌──────────────── TIME LOOP (per Ra) ─────────────────────────┐
 pMesh,uMesh,vMesh,nu=Pr │ 1. momentum.advance(u*,v*)          (ν∇²u − (u·∇)u, AB2)    │
 ───────────────────────→│ 2. addBuoyancy → v* += dt·Ra·Pr·θ                            │
 T (=θ), Ra, Pr          │ 3-5. Poisson(CG) + projectVelocity  (∇·u = 0)               │
 ───────────────────────→│ 6. energy.advance(θ, u*, v*)        (α∇²θ − (u·∇)θ, AB2)    │
                         └──────────────────────────────────────────────────────────────┘
                                         │
                                         └─→ output/Ra_*/{temperature,u,v,pressure}.vtk
                                                        dhc_metrics.csv (u_max,v_max,Nu)
```

---

### Key files

| Piece | File |
|-------|------|
| Main loop and mode selection | [src/main.cpp](../../src/main.cpp) |
| Mesh generation + stretching | [src/mesh/MeshGenerator.cpp](../../src/mesh/MeshGenerator.cpp) |
| Mesh structures and config | [src/mesh/Mesh_config.h](../../src/mesh/Mesh_config.h) |
| Energy equation (θ) | [src/EnergySolver.cpp](../../src/EnergySolver.cpp) · [.hpp](../../src/EnergySolver.hpp) |
| Boussinesq buoyancy | [src/Buoyancy.hpp](../../src/Buoyancy.hpp) |
| Momentum predictor | [src/MomentumSolver.cpp](../../src/MomentumSolver.cpp) |
| Pressure / projection | [src/PressurePoissonSystem.cpp](../../src/PressurePoissonSystem.cpp) |
```
</content>
