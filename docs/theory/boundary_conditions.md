# Boundary Conditions — Types, Input Format and Implementation

This page is the single reference for **every boundary condition the solver
supports**: how it is written in the input files, what it imposes on each
variable ($u$, $v$, $p$, $\theta$), and where in the code it is enforced.
The scenario pages ([../scenarios/](../scenarios/)) only list *which*
conditions each problem uses and link back here.

Related pages: [equations.md](equations.md) (the equations these conditions
close), [mesh.md](mesh.md) (where the variables live on the staggered grid),
[solver.md](solver.md) (the order of operations in a time step).

---

## Contents

1. [Where boundary conditions are specified](#1-where-boundary-conditions-are-specified)
2. [Boundary-condition types](#2-boundary-condition-types)
3. [Velocity: normal component](#3-velocity-normal-component-nodes-on-the-boundary)
4. [Velocity: tangential component (ghost-cell image method)](#4-velocity-tangential-component-ghost-cell-image-method)
5. [Pressure](#5-pressure)
6. [Temperature](#6-temperature-dhc-only)
7. [Immersed solid obstacle](#7-immersed-solid-obstacle-square-cylinder-only)
8. [When the conditions are applied in a time step](#8-when-the-conditions-are-applied-in-a-time-step)
9. [Summary per scenario](#9-summary-per-scenario)
10. [Limitations](#10-limitations-what-is-not-implemented)

---

## 1. Where boundary conditions are specified

| File | Controls | Read by |
|---|---|---|
| `input/<scenario>/boundaries.txt` | velocity (and, implicitly, pressure) on the four outer sides | [`MeshConfig::define_boundaries`](../../src/mesh/MeshGenerator.cpp#L996) |
| `input/differentially_heated_cavity/boundaries_T.txt` | temperature on the four outer sides | [`MeshConfig::define_temperature_boundaries`](../../src/mesh/MeshGenerator.cpp#L1086) |
| `input/square_cylinder*/obstacle.txt` | position of the solid obstacle (its surface is a no-slip wall) | `read_obstacle_bounds` in [main.cpp](../../src/main.cpp#L54) → [`MeshConfig::mark_solid_block`](../../src/mesh/MeshGenerator.cpp#L533) |

Both `boundaries*.txt` files use the same line format:

```
# BOCO_NAME   FACE    IDX_START  IDX_END   TYPE               VALUE
lid           NORTH   0          9999      WALL_FIXED_VALUE   1.0
```

| Column | Meaning |
|---|---|
| `BOCO_NAME` | free label (printed in the log, stored in `BoundaryFace::boco`) |
| `FACE` | `WEST` ($x = 0$), `EAST` ($x = L$), `SOUTH` ($y = 0$), `NORTH` ($y = H$) |
| `IDX_START IDX_END` | half-open range `[start, end)` of boundary faces along that side, 0-based (cell rows for WEST/EAST, cell columns for SOUTH/NORTH). `0 9999` = the whole side, whatever the mesh size (the end is clamped) |
| `TYPE` | keyword from §2 |
| `VALUE` | optional number (default 0.0); its meaning depends on `TYPE` |

Lines starting with `#` and empty lines are ignored. Each face of the side
gets its own `BoundaryFace` entry (`type`, `value`, position, area) in
`GridData::bound_west/east/south/north` (velocity) or `bound_*_T`
(temperature) — see [Mesh_config.h](../../src/mesh/Mesh_config.h).

---

## 2. Boundary-condition types

Keywords accepted in the `TYPE` column
([`string_to_bc_enum`](../../src/mesh/MeshGenerator.cpp#L984)) and the
`BCType` enum value they map to
([Mesh_config.h](../../src/mesh/Mesh_config.h)):

| Keyword | `BCType` | Meaning of `VALUE` | Allowed faces | Used by |
|---|---|---|---|---|
| `WALL` | `WALL_ADIABATIC` | ignored | all | all scenarios (no-slip walls) |
| `WALL_FIXED_VALUE` | `WALL_FIXED_VALUE` | tangential wall velocity (moving wall) | NORTH/SOUTH in practice | lid-driven cavity (the lid) |
| `INLET` | `INLET` | uniform normal inflow velocity $U_\infty$ | WEST/EAST only | `square_cylinder` |
| `INLET_PARABOLIC` | `INLET` (per-face value) | centreline velocity $U_{\max}$ | WEST/EAST only | `square_cylinder_breuer` |
| `OUTLET` | `OUTLET` | must be 0 (see §5) | WEST/EAST only (velocity side) | both cylinder scenarios |
| `SYMMETRY` | `SYMMETRY` | — | — | **not implemented: behaves as a no-slip wall** (§10) |
| `WALL_ISOTHERMAL` | `WALL_ISOTHERMAL` | wall temperature $\theta_w$ | all | `boundaries_T.txt` only |
| `WALL_ADIABATIC` | `WALL_ADIABATIC` | ignored | all | `boundaries_T.txt` (and as an alias of `WALL`) |

What each type imposes on each variable:

| Type | $u_n$ (normal velocity) | $u_t$ (tangential velocity) | $p$ | $\theta$ |
|---|---|---|---|---|
| `WALL` | $0$ | $0$ (no-slip) | $\partial p/\partial n = 0$ | — |
| `WALL_FIXED_VALUE` | $0$ | `VALUE` | $\partial p/\partial n = 0$ | — |
| `INLET` | `VALUE` | $0$ | $\partial p/\partial n = 0$ | — |
| `INLET_PARABOLIC` | $4U_{\max}(y-y_0)(y_1-y)/H^2$ | $0$ | $\partial p/\partial n = 0$ | — |
| `OUTLET` | $\partial u_n/\partial n = 0$ | $0$ (see §4) | $p = 0$ | — |
| `WALL_ISOTHERMAL` | — | — | — | $\theta = $ `VALUE` |
| `WALL_ADIABATIC` (in `boundaries_T.txt`) | — | — | — | $\partial\theta/\partial n = 0$ |

---

## 3. Velocity: normal component (nodes on the boundary)

On the staggered grid the **normal** velocity of each outer side is stored
*on* the side itself: $u$-nodes $k = 1$ and $k = N{+}1$ lie on the WEST and
EAST sides, $v$-nodes $l = 1$ and $l = M{+}1$ on the SOUTH and NORTH sides
([mesh.md §2](mesh.md#2-three-meshes-one-geometry)). These nodes are not
unknowns: they are overwritten by
[`MomentumSolver::applyBoundaryConditions`](../../src/MomentumSolver.cpp#L84)
and the projection never corrects them.

| Side | Node | Value set |
|---|---|---|
| WEST | $u_{1,j}$ | `INLET`: `bound_west[j-1].value`; `OUTLET`: $u_{2,j}$ (zero gradient); any wall: $0$ |
| EAST | $u_{N+1,j}$ | `INLET`: `bound_east[j-1].value`; `OUTLET`: $u_{N,j}$ (zero gradient); any wall: $0$ |
| SOUTH | $v_{i,1}$ | always $0$ |
| NORTH | $v_{i,M+1}$ | always $0$ |

- The WEST/EAST values are taken **face by face** (`bound_*[j-1]`), which is
  what makes the parabolic inlet possible: `define_boundaries` stores
  $4U_{\max}(y_f-y_0)(y_1-y_f)/H^2$ at each face centre $y_f$
  ([MeshGenerator.cpp:1051-1065](../../src/mesh/MeshGenerator.cpp#L1051-L1065)).
- SOUTH/NORTH are **always impermeable**: an `INLET` or `OUTLET` on those
  sides is parsed but has no effect on $v$ (§10).
- The zero-gradient outlet copies the last interior value
  ($u_{N+1,j} = u_{N,j}$); it is refreshed every time the BCs are applied
  (§8), so it follows the interior flow.

---

## 4. Velocity: tangential component (ghost-cell image method)

The **tangential** velocity is not stored on the side: the last $u$-node row
lies half a cell inside the SOUTH/NORTH walls, and the last $v$-node column
half a cell inside the WEST/EAST walls. Their diffusion and convection
stencils need a value beyond the wall, supplied by an **image (ghost) node**
mirrored across the wall:

$$
\phi_{ghost} = 2\,\phi_{wall} - \phi_P, \qquad
\Delta_{P\to ghost} = 2\,(\text{distance from } P \text{ to the wall})
$$

so that the linear interpolation between $P$ and its ghost gives exactly
$\phi_{wall}$ at the wall. The second derivative at the wall-adjacent node
then reads, with $h$ = distance to the wall,

$$
\frac{\partial^2\phi}{\partial y^2}\bigg|_P \approx
\frac{2}{\Delta_n + 2h}\left[\frac{\phi_N-\phi_P}{\Delta_n} -
\frac{\phi_P-(2\phi_{wall}-\phi_P)}{2h}\right].
$$

> Using $\phi_{ghost} = \phi_{wall}$ (without the factor 2) gives the wrong
> wall condition; this was one of the early bugs, see
> [../reference/development_log.md](../reference/development_log.md).

The wall value $\phi_{wall}$ is chosen once per case, in the
`MomentumSolver` constructor, by
[`ghostValue`](../../src/MomentumSolver.cpp#L35):

| Type of the side | $\phi_{wall}$ (tangential) |
|---|---|
| `WALL_FIXED_VALUE` | `VALUE` (e.g. $U_{lid} = 1$) |
| `WALL`, `INLET`, `INLET_PARABOLIC`, `OUTLET`, anything else | $0$ |

| Ghost | Built from | Used in |
|---|---|---|
| $u$ above NORTH / below SOUTH | `bound_north[0]`, `bound_south[0]` | `uRHS` ([MomentumSolver.cpp:200-222](../../src/MomentumSolver.cpp#L200-L222)), `DiffusionOperator::diffuseU` |
| $v$ left of WEST / right of EAST | `bound_west[0]`, `bound_east[0]` | `vRHS` ([MomentumSolver.cpp:300-319](../../src/MomentumSolver.cpp#L300-L319)), `DiffusionOperator::diffuseV` |

Consequences worth knowing:

- The tangential wall value is **uniform along each side** (taken from the
  first face of the list): a side cannot be part moving, part fixed.
- At an `INLET` the tangential velocity is $v = 0$ — correct for a uniform or
  parabolic inflow, which has no cross-stream component.
- At an `OUTLET` the tangential velocity is also imposed as $v = 0$ (a
  Dirichlet condition), not as a zero-gradient condition. This is a modelling
  choice of the present code; the cylinder outlet is placed 37 $D$ downstream
  so that the wake reaches it nearly parallel.

---

## 5. Pressure

Pressure needs no input line of its own: its condition follows from the type
of each side, inside
[`LaplaceOperator::assembleToTriplets`](../../src/operators/LaplaceOperator.cpp)
(called by `PressurePoissonSystem::assembleMatrix`).

**Walls and inlets — homogeneous Neumann, $\partial p/\partial n = 0$.**
In the projection method the normal velocity on these sides is already
prescribed and is never corrected, so the pressure gradient across them must
not act: the boundary face is simply left out of the stencil (no
off-diagonal entry, no diagonal contribution).

**Outlet — Dirichlet, $p = 0$ on the outlet face.** The face still
contributes $a_{face} = S/\Delta$ to the diagonal but gets no off-diagonal
entry ([LaplaceOperator.cpp:66-95](../../src/operators/LaplaceOperator.cpp#L66-L95)).
The right-hand-side term $-a_{face}\,p_{bc}$ is **not** added, i.e. the
outlet pressure is always 0 whatever `VALUE` says — keep it at `0.0`.

**Reference pressure.** With Neumann conditions on every side (closed
cavities) the matrix is singular — $p$ is defined up to a constant. The
system is regularised by fixing $p = 0$ in cell 0 (the SOUTH-WEST corner
cell) — zeroing row *and* column 0 to keep the matrix symmetric for CG
([`LaplaceOperator::pinReference`](../../src/operators/LaplaceOperator.cpp#L113),
called from [main.cpp:197-201](../../src/main.cpp#L197-L201), with
`b[0] = 0` re-imposed every step at
[main.cpp:378](../../src/main.cpp#L378)). The cylinder scenarios skip the
pinning: the Dirichlet outlet already fixes the constant.

| Scenario | Pressure on the sides | Reference | Linear solver |
|---|---|---|---|
| Lid-driven cavity | Neumann ×4 | pinned, $p_{(1,1)} = 0$ | CG |
| Differentially heated cavity | Neumann ×4 | pinned, $p_{(1,1)} = 0$ | CG |
| Square cylinder (both) | Neumann (inlet, walls, obstacle); Dirichlet $p = 0$ (outlet) | outlet | Jacobi-PCG |

---

## 6. Temperature (DHC only)

Temperature $\theta$ lives at P-cell centres, like pressure. Its conditions
are read from `boundaries_T.txt` into the parallel lists `bound_*_T`, which
`define_temperature_boundaries` first fills with `WALL_ADIABATIC` copies of
the velocity faces — **any side not listed is adiabatic**.

[`EnergySolver::tRHS`](../../src/EnergySolver.cpp#L48) applies them with the
same image method as the velocity
([EnergySolver.cpp:59-75](../../src/EnergySolver.cpp#L59-L75)):

| Type | Ghost value | Imposes |
|---|---|---|
| `WALL_ISOTHERMAL` | $\theta_{ghost} = 2\theta_w - \theta_P$, distance $2h$ | $\theta = \theta_w$ on the wall (Dirichlet) |
| `WALL_ADIABATIC` | $\theta_{ghost} = \theta_P$ | $\partial\theta/\partial n = 0$ (no heat flux) |

The convective flux through every outer face is zero because the wall-normal
velocity is zero there, so convection needs no ghost. As for the velocity,
the type and value are taken from the first face of each side (uniform along
the side).

---

## 7. Immersed solid obstacle (Square Cylinder only)

The obstacle is not a boundary of the domain: it is a block of P-cells
flagged `is_solid` inside the mesh
([`mark_solid_block`](../../src/mesh/MeshGenerator.cpp#L533)). Its bounds
from `obstacle.txt` are snapped to the nearest existing face coordinates, so
the surface always coincides with cell faces; a warning is printed if the
snap moves a bound by more than $10^{-6}$. The obstacle surface is a
**no-slip, impermeable, adiabatic-for-pressure wall**:

| Variable | Treatment | Where |
|---|---|---|
| $u$, $v$ on the surface (normal) | $0$; nodes with exactly one solid neighbour cell ([`uNodeOnObstacleSurface`](../../src/mesh/ObstacleMask.hpp)) | `applyBoundaryConditions` |
| $u$, $v$ inside the solid | held at $0$, never evaluated | `applyBoundaryConditions`, `advance` skips them |
| $u$, $v$ tangential | image ghost $\phi_{ghost} = -\phi_P$ (wall value 0) at distance $2h$ to the **facing** obstacle side (a solid neighbour to the east means the wall is the obstacle's west face, `obstacle_xmin`) | `uRHS`/`vRHS`, `DiffusionOperator` |
| $p$ | Neumann: `mark_solid_block` sets `bc_* = true` on both sides of the interface, so the Laplacian skips those faces; solid cells get an identity row ($p = 0$, decoupled) | `LaplaceOperator` ([line 27](../../src/operators/LaplaceOperator.cpp#L27)) |
| $\nabla\cdot u^*$ | 0 in solid cells | [DivergenceOperator.cpp:18](../../src/operators/DivergenceOperator.cpp#L18) |

The helpers deciding, at each stencil lookup, whether a neighbour is fluid,
on the surface, or inside the solid are in
[src/mesh/ObstacleMask.hpp](../../src/mesh/ObstacleMask.hpp). For the two
cavity scenarios no cell is solid, so all of these branches are inactive.

---

## 8. When the conditions are applied in a time step

```
for each step (src/main.cpp, time loop):
    momentum.advance(u, v, dt)                 # predictor: interior nodes only
    momentum.applyBoundaryConditions(u*, v*)   # (1) BCs on u* BEFORE the PPE
    addBuoyancy(...)                           # DHC only, interior v* nodes
    b = -(rho/dt) * div(u*) * V                # uses the boundary values set in (1)
    p = CG/PCG solve(A, b)
    projectVelocity(p, ...)                    # corrects interior nodes only
    momentum.applyBoundaryConditions(u, v)     # (2) BCs on u^{n+1}
    energy.advance(T, u, v, dt)                # DHC only, temperature ghosts inside
```

- Call **(1)** ([main.cpp:372](../../src/main.cpp#L372)) makes sure the
  divergence that drives the pressure equation sees the correct wall,
  inlet and outlet fluxes. Without it, the predictor's values on the
  obstacle surface leaked into the divergence of the neighbouring cells
  (a bug fixed during the cylinder validation, see
  [../reference/development_log.md](../reference/development_log.md)).
  For the cavities it changes nothing (the wall values are fixed).
- Call **(2)** ([main.cpp:388](../../src/main.cpp#L388)) resets the
  boundary nodes after the projection and refreshes the zero-gradient outlet.
- The tangential ghosts (§4) and temperature ghosts (§6) are not stored in
  the fields; they are evaluated inside the stencils every time.

---

## 9. Summary per scenario

| | WEST | EAST | SOUTH | NORTH | Obstacle |
|---|---|---|---|---|---|
| **Lid-driven cavity** | `WALL` | `WALL` | `WALL` | `WALL_FIXED_VALUE 1.0` | — |
| **Diff. heated cavity**, velocity | `WALL` | `WALL` | `WALL` | `WALL` | — |
| **Diff. heated cavity**, temperature | `WALL_ISOTHERMAL 1.0` | `WALL_ISOTHERMAL 0.0` | `WALL_ADIABATIC` | `WALL_ADIABATIC` | — |
| **Square cylinder** (`square_cylinder`) | `INLET 1.0` | `OUTLET` | `WALL` | `WALL` | no-slip |
| **Square cylinder, Breuer** (`square_cylinder_breuer`) | `INLET_PARABOLIC 1.0` | `OUTLET` | `WALL` | `WALL` | no-slip |

The scenario pages give the same conditions in mathematical form, with a
sketch of the domain.

---

## 10. Limitations (what is not implemented)

| Limitation | Effect | Where to change it |
|---|---|---|
| `SYMMETRY` (slip wall) is parsed but `ghostValue` falls through to 0 and the normal velocity is set to 0 → it behaves exactly as a **no-slip wall** | blocks the unconfined cylinder benchmark (Sohankar et al. 1998), which needs slip side walls | `MomentumSolver::ghostValue`, `applyBoundaryConditions` |
| `INLET`/`OUTLET` only act on WEST/EAST; $v$ on SOUTH/NORTH is always 0 | no vertical inflow/outflow | `applyBoundaryConditions` |
| Tangential velocity ghost and temperature condition are uniform per side (first face) | no partially moving walls or partially heated walls | `MomentumSolver` constructor, `EnergySolver::extractWallT` |
| Outlet pressure is always 0 (`VALUE` ignored) and the outlet tangential velocity is $v = 0$ | fine for a far outlet; a convective outlet is not available | `LaplaceOperator`, `ghostValue` |
| Only one rectangular obstacle, aligned with the mesh | no multiple bodies or curved shapes | `mark_solid_block` |
| No energy equation with an obstacle or open boundaries | DHC only | `EnergySolver` |
