# The Staggered (MAC) Grid

Theory notes on the mesh: why it is staggered, how it is built, how
stretching/refinement works, and how boundary conditions attach to it. For
the equations being discretized, see [equations.md](equations.md); for how
the mesh is consumed by the solver classes, see [solver.md](solver.md).

---

## 1. Why a staggered grid?

- On a **collocated** grid (all variables at cell centers), the discrete
  pressure gradient and divergence operators have a null space of
  oscillatory "checkerboard" pressure modes — pressure can oscillate
  node-to-node without producing any force, so it is never damped
  (checkerboarding).
- The **Marker-and-Cell (MAC) staggered grid** avoids this by storing
  different variables at different locations within each cell:
  - **Pressure** $p$ (and, for the heated cavity, **temperature** $\theta$) at
    the cell **center**.
  - **$u$-velocity** at the vertical **east/west faces**.
  - **$v$-velocity** at the horizontal **north/south faces**.
- This placement is also the natural one for finite-volume fluxes: the face
  velocities ARE the normal velocities needed to compute the flux through
  that face — no interpolation required for the divergence or convective
  terms.

---

## 2. Three meshes, one geometry

The solver builds **three separate logical grids** from one set of face
coordinates, all sharing the same physical domain:

```cpp
GridData pMesh = myConfig.generate_mesh();        // pressure / temperature
GridData uMesh = myConfig.generate_u_mesh(pMesh);  // u-velocity
GridData vMesh = myConfig.generate_v_mesh(pMesh);  // v-velocity
```
([src/main.cpp:85-87](../../src/main.cpp#L85-L87))

| Grid | Lives at | Size (N×M pressure cells) | Role |
|---|---|---|---|
| `pMesh` | cell centers | $N \times M$ | pressure $p$, temperature $\theta$ |
| `uMesh` | vertical (east/west) faces | $(N{+}1) \times M$ | $u$-velocity |
| `vMesh` | horizontal (north/south) faces | $N \times (M{+}1)$ | $v$-velocity |

- **`uMesh` construction**
  ([src/mesh/MeshGenerator.cpp:333-385](../../src/mesh/MeshGenerator.cpp#L333-L385)):
  shares `pMesh`'s Y faces exactly, and places U-cell boundaries at the
  **centers** of the P cells in X — so U node $k$ sits exactly on the shared
  face between P cells $(k-1)$ and $k$. This is what makes the staggering
  exact and consistent: no extra interpolation is needed to know where a U
  node lies relative to its neighboring P cells.
- **`vMesh` construction** is the symmetric analogue in Y
  ([src/mesh/MeshGenerator.cpp:391+](../../src/mesh/MeshGenerator.cpp#L391)).

### Indexing convention (used identically in every solver class)

1-based logical indices, mapped to 0-based flat storage, row-major
(j/l-major, i/k-minor):

```cpp
pIdx(i, j) = (j-1)*N_cells_x + (i-1)     // P cell  (i=1..N,   j=1..M)
uIdx(k, j) = (j-1)*N_cells_x + (k-1)     // U node  (k=1..N+1, j=1..M)
vIdx(i, l) = (l-1)*N_cells_x + (i-1)     // V node  (i=1..N,   l=1..M+1)
```

- Appears verbatim and consistently in
  [PressurePoissonSystem.hpp](../../src/PressurePoissonSystem.hpp),
  [MomentumSolver.hpp](../../src/MomentumSolver.hpp), and
  [EnergySolver.hpp](../../src/EnergySolver.hpp) — this shared convention is
  what lets all the solver pieces address the same physical location
  consistently without passing coordinates around.
- Staggered face mapping for a P cell $(i,j)$:
  - East face $\leftrightarrow$ U node $k = i+1$; West face $\leftrightarrow$ U node $k=i$.
  - North face $\leftrightarrow$ V node $l = j+1$; South face $\leftrightarrow$ V node $l=j$.

---

## 3. Geometric data stored per cell

`generate_mesh()` computes and stores, for every cell
([src/mesh/MeshGenerator.cpp:152-260](../../src/mesh/MeshGenerator.cpp#L152-L260)):

- `xvc[i]`, `yvc[j]` — **face** coordinates (cell boundaries).
- `x[i]`, `y[j]` — cell **centroids**, including one layer of **ghost
  cells** (index 0 and N+1 / M+1) mirrored just outside the domain, used by
  the ghost-cell image method (see [equations.md §4](equations.md#4-boundary-conditions--the-ghost-cell-image-method)).
- `Se, Sw, Sn, Ss` — face **areas** (`Δy·W_depth` for E/W faces, `Δx·W_depth`
  for N/S faces) — the depth `W_depth` makes a "2D" simulation a true 3D
  finite-volume problem with one cell in the third direction.
- `V` — cell **volume** (`Δx·Δy·W_depth`).

These fields are **the contract** between the mesh and every operator/solver:
`DiffusionOperator`, `GradientOperator`, `DivergenceOperator`,
`LaplaceOperator`, `EnergySolver` all read centroids and face areas directly
from `GridData` rather than recomputing geometry themselves.

---

## 4. Mesh stretching (refinement)

Refinement towards walls (where gradients/boundary layers are steepest) is
done with a `sinh`-based mapping in
[`generate_stretched_coords`](../../src/mesh/MeshGenerator.cpp#L79):

- **Uniform case:** if $|\beta| < 10^{-6}$, cells are evenly spaced,
  $\Delta x = L/N$.
- **Refine at the end** ($x=L$) if $\beta > 0$.
- **Refine at the start** ($x=0$) if $\beta < 0$.
- **Algorithm**, for `num_cells` cells over `[0, length]`:
  1. Generate a uniform logical coordinate $\xi_i \in [-1, 1]$.
  2. Shift and evaluate a $\sinh$ function centered on the focus point,
     then re-normalize it back to $[-1, 1]$:
     $$
     g(\xi) = \frac{\sinh(\beta(\xi - \xi_{focus})) - v_{min}}{v_{max}-v_{min}}
     \cdot 2 - 1
     $$
  3. Map the normalized coordinate to physical space:
     $x = \tfrac{L}{2}(1 + g(\xi)) + \text{offset}$.
- Larger $|\beta|$ concentrates more cells near the focus point. In the
  differentially heated cavity this lets the mesh resolve the thin thermal
  boundary layers at the hot/cold walls, where the Nusselt number is measured.
- Each of `beta_x_west`, `beta_x_east`, `beta_y_north`, `beta_y_south`
  (`input/stretch.txt`) controls stretching **independently per zone/side**.

> **Known accuracy caveat:** `LaplaceOperator`'s current FVM stencil uses
> `a_E = Se / (x[i+1] - x[i])`, which is only first-order accurate when the
> mesh is non-uniform (stretched) — confirmed by the MMS convergence test
> ("Suite 10" in [docs/OPERATORS_REPORT.md](../OPERATORS_REPORT.md)), which
> measures rate ≈ 1.0 instead of the ≈ 2.0 achieved on uniform grids. The fix
> (replacing `1/dx` with `2/(h_e+h_w)` in the stencil) is known but not yet
> applied — see [solver.md §3](solver.md#3-pressurepoissonsystem-the-elliptic-constraint)
> for the implication on production runs.

---

## 5. Zones: the domain is split into a 3×3 grid of regions

- The domain is partitioned in X into **west / center / east** zones and in Y
  into **south / center / north** zones, each stretched independently and
  then merged into one set of global face coordinates
  ([src/mesh/MeshGenerator.cpp:195-229](../../src/mesh/MeshGenerator.cpp#L195-L229)).
- Zone boundaries are set by four corner nodes in
  [input/data.txt](../../input/data.txt): `node_NW`, `node_NE`, `node_SW`,
  `node_SE`.
- **For both scenarios currently implemented** (lid-driven cavity,
  differentially heated cavity), the corner nodes are placed at the domain
  corners and `center`/`east` zone cell counts are set to 0 in
  `input/stretch.txt`, which collapses the partition to a **single** west/south
  block covering the whole domain — i.e. a plain rectangular cavity, just
  built through a more general 3-zone machinery.
- This 3×3 zone layout is scaffolding that would allow concentrating cells
  around a **non-degenerate center zone** (e.g. around an interior obstacle).
  It is not yet paired with any interior/immersed-boundary condition
  handling — see
  [../scenarios/square_cylinder.md](../scenarios/square_cylinder.md) for what
  such a scenario would need.

---

## 6. Connectivity and boundary definition

- [`build_connectivity`](../../src/mesh/MeshGenerator.cpp#L448) builds, for
  every cell, its East/West/North/South neighbor id (or `-1` with a boundary
  flag set if the cell touches a domain edge). Used mainly for inspection
  and future general-stencil needs; the production solvers address
  neighbors directly via the `pIdx`/`uIdx`/`vIdx` formulas in §2 rather than
  walking this table.
- [`define_boundaries`](../../src/mesh/MeshGenerator.cpp#L915) reads
  [input/boundaries.txt](../../input/boundaries.txt) and populates
  `mesh.bound_west/east/north/south` — one `BoundaryFace` entry per boundary
  cell on that side, carrying a `BCType` and a prescribed value.
  - File format: `BOCO_NAME  FACE  IDX_START  IDX_END  TYPE  VALUE`, where
    the index range is **half-open** `[start, end)` — `0 9999` covers an
    entire side regardless of mesh resolution.
- [`define_temperature_boundaries`](../../src/mesh/MeshGenerator.cpp#L993)
  does the same for temperature, reading
  [input/boundaries_T.txt](../../input/boundaries_T.txt) into the parallel
  `bound_*_T` lists. It **copies geometry** (position, area) from the
  velocity boundary lists and only overwrites type/value — faces not
  mentioned default to adiabatic (`WALL_ADIABATIC`).
- **`BCType`** enum
  ([src/mesh/Mesh_config.h:18-26](../../src/mesh/Mesh_config.h#L18-L26)):
  `WALL_ADIABATIC`, `WALL_FIXED_VALUE`, `INLET`, `OUTLET`, `SYMMETRY`,
  `WALL_ISOTHERMAL`, `UNDEFINED`. Velocity BCs mainly use `WALL_FIXED_VALUE`
  (e.g. the moving lid) and plain `WALL` (no-slip, value 0); temperature BCs
  use `WALL_ISOTHERMAL` (Dirichlet $\theta$) and `WALL_ADIABATIC` (zero flux).

---

## 7. VTK export (for ParaView)

- `export_mesh_to_VTK` and the field-writer in
  [src/main.cpp](../../src/main.cpp) write a **rectilinear grid** VTK file
  (`DATASET RECTILINEAR_GRID`) with one scalar per cell
  (`CELL_DATA` / `SCALARS` / `LOOKUP_TABLE`).
- Because the grid is rectilinear (not a general unstructured mesh), only
  the 1D face coordinate arrays (`X_COORDINATES`, `Y_COORDINATES`) need to be
  written — ParaView reconstructs the full 2D grid from the tensor product.
- Outputs per case land in `output/Re_<value>/` or `output/Ra_<value>/`:
  `pressure.vtk`, `u_velocity.vtk`, `v_velocity.vtk`, and (DHC only)
  `temperature.vtk`.
