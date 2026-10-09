# Scenario: Square Cylinder Problem (Planned / Not Yet Implemented)

> **Status: not implemented.** This document describes a physical problem the
> mesh infrastructure is partially scaffolded for, but which has no working
> solver support yet. It is included so contributors know what exists, what's
> missing, and what the physical validation target would be.

---

## 1. The physical problem

Uniform flow past a square cylinder (a solid square obstacle) placed inside
a channel or open domain — a standard bluff-body benchmark:

- At low Reynolds number, flow separates smoothly around the obstacle with
  two steady, symmetric recirculation zones behind it.
- Above a critical Reynolds number ($Re \approx 45$–60 depending on blockage
  ratio), the wake becomes unstable and starts **periodically shedding
  vortices** (a von Kármán vortex street), producing oscillating lift and
  drag forces on the cylinder.
- The canonical validation metric is the **Strouhal number**
  $St = f L / U_\infty$ (shedding frequency vs. inflow velocity and obstacle
  size) as a function of $Re$, together with time-averaged drag/lift
  coefficients — compared against well-established numerical/experimental
  literature values.

## 2. What already exists that's relevant

The mesh generator supports a **3×3 zone domain partition** (west/center/east
× south/center/north), controlled by four corner nodes in
[input/data.txt](../../input/data.txt) (`node_NW/NE/SW/SE`) and per-zone cell
counts/stretching in [input/stretch.txt](../../input/stretch.txt) — see
[../theory/mesh.md §5](../theory/mesh.md#5-zones-the-domain-is-split-into-a-33-grid-of-regions).

- For the two currently implemented scenarios, the `center`/`east` zones are
  deliberately collapsed to zero cells, degenerating the partition into one
  uniform/stretched block.
- With non-zero `center` zone dimensions, this machinery **could** carve out
  and independently refine a rectangular region in the middle of the domain
  — structurally, exactly the footprint a square obstacle would occupy, with
  extra mesh density concentrated around it via the existing `beta_*_center`
  stretching parameters.

## 3. What is missing

None of the following exists in the codebase today:

- **Interior/immersed boundary conditions.** `define_boundaries` and
  `build_connectivity`
  ([src/mesh/MeshGenerator.cpp](../../src/mesh/MeshGenerator.cpp)) only
  recognize boundary faces on the four **outer** sides of the domain
  (`WEST`/`EAST`/`NORTH`/`SOUTH`). There is no concept of an interior solid
  region whose cells should be excluded from the fluid solve, nor of
  no-slip conditions applied on an interior face.
- **Inlet/outlet treatment for an open channel.** `BCType::INLET` and
  `BCType::OUTLET` exist in the enum
  ([src/mesh/Mesh_config.h:18-26](../../src/mesh/Mesh_config.h#L18-L26)) but
  are not exercised by either current scenario (both are closed cavities)
  — their handling in `PressurePoissonSystem` (Dirichlet pressure at an
  outlet, prescribed inflow profile) is only partially sketched in comments,
  not implemented/validated.
- **Force/Strouhal post-processing.** No code currently integrates pressure
  and viscous stress over an obstacle surface to produce drag/lift
  coefficients, nor an FFT or peak-detection routine to extract a shedding
  frequency.
- **Unsteady-flow-appropriate validation tooling.** The existing plotting
  scripts (`plotters/plot_results.py`, `plotters/plot_dhc.py`) are built
  around steady-state centerline profiles / Nusselt numbers; a vortex-shedding
  case would need time-series force plots and frequency-domain analysis.

## 4. What implementing this would require (rough scope)

1. Extend `define_boundaries`/`build_connectivity` (or add a parallel
   mechanism) to mask out an interior rectangular block of P/U/V cells and
   apply no-slip conditions on its four faces, reusing the existing
   ghost-cell image method conceptually but now for an **interior** wall
   instead of only exterior ones.
2. Implement and validate inlet (prescribed/uniform $u$) and outlet
   (zero-gradient or convective) boundary conditions in
   `PressurePoissonSystem` and `MomentumSolver` for an open channel, since a
   closed cavity has no meaningful vortex shedding.
3. Add surface-integration post-processing for drag/lift, and a
   frequency-domain (or peak-counting) post-processor for the Strouhal
   number.
4. Validate against known literature values (e.g. $St \approx 0.13$–0.15 at
   $Re \sim 100$–200 for a square cylinder, blockage-ratio dependent) the
   same way the other two scenarios are validated against Ghia et al. and
   de Vahl Davis.

## 5. See also

- [../theory/mesh.md](../theory/mesh.md) — the zone/stretching system this
  would build on.
- [lid_driven_cavity.md](lid_driven_cavity.md),
  [differentially_heated_cavity.md](differentially_heated_cavity.md) — the
  two scenarios that are implemented and validated today.
