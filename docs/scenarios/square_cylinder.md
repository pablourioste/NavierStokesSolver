# Scenario: Square Cylinder Problem

> **Status: implemented, validation pending.** The solver runs this scenario
> end to end — immersed solid obstacle, inlet/outlet channel flow, drag/lift
> recording, and a Strouhal-number estimate — without NaN or divergence on
> the shipped mesh. **The Strouhal number has not yet been calibrated against
> the literature** (~0.13–0.15 at Re ~ 100–200); matching it is future work
> (mesh/`dt`/domain-length tuning), not something this round attempted. Treat
> any `St` value this solver reports as provisional.

---

## 1. The physical problem

Uniform flow past a square cylinder (a solid square obstacle) placed inside
a channel — a standard bluff-body benchmark:

- At low Reynolds number, flow separates smoothly around the obstacle with
  two steady, symmetric recirculation zones behind it.
- Above a critical Reynolds number ($Re \approx 45$–60 depending on blockage
  ratio), the wake becomes unstable and starts **periodically shedding
  vortices** (a von Kármán vortex street), producing oscillating lift and
  drag forces on the cylinder.
- The canonical validation metric is the **Strouhal number**
  $St = f L / U_\infty$ (shedding frequency vs. inflow velocity and obstacle
  size) as a function of $Re$, together with time-averaged drag/lift
  coefficients.

## 2. How to run it

1. Set `scenario = square_cylinder` in [input/scenario.txt](../../input/scenario.txt).
2. Geometry lives in [input/square_cylinder/](../../input/square_cylinder/):
   - `data.txt` — channel domain (default: length $3.125$, height $1.0$,
     blockage ratio $1/8$) and the zone-corner placement of the obstacle.
   - `obstacle.txt` — the obstacle's exact bounds
     (`obstacle_x0/x1/y0/y1`), read by `main.cpp` and passed to
     `MeshConfig::mark_solid_block()`.
   - `stretch.txt` — per-zone cell counts/stretching, concentrating
     resolution around the obstacle via the (now non-degenerate) center zone.
   - `sim_params.txt` — `U_inf` (free-stream inlet speed), `dt`, `max_steps`
     (large — shedding never reaches steady state, so `ss_tol = 0.0`).
   - `boundaries.txt` — `WEST = INLET`, `EAST = OUTLET`, `NORTH`/`SOUTH = WALL`.
   - `reynolds_cases.txt` — the Re sweep.
3. Build and run:
   ```bash
   bash run.sh mesh_gen
   ```
4. Visualize:
   ```bash
   python3 plotters/plot_square_cylinder.py
   ```

## 3. What the implementation does

- **Immersed solid obstacle**: `InternalCell::is_solid` (set by
  `MeshConfig::mark_solid_block()`, called only for this scenario) marks
  P-cells inside the obstacle footprint. `LaplaceOperator` decouples solid
  cells with a trivial identity pressure row; `DivergenceOperator` zeroes
  their divergence; `MomentumSolver`/`DiffusionOperator`'s stencils treat a
  blocked neighbor exactly like an outer wall — same ghost-image, no-slip
  method, now triggered by the obstacle surface in addition to the domain
  edge. See [../theory/solver.md](../theory/solver.md) and
  [src/mesh/ObstacleMask.hpp](../../src/mesh/ObstacleMask.hpp) for the shared
  helpers this relies on.
- **Inlet/outlet**: `WEST` prescribes a uniform normal velocity (`U_inf`);
  `EAST` uses a zero-gradient (extrapolated) outflow velocity and a Dirichlet
  pressure row ($p=0$ gauge) instead of the closed-cavity homogeneous-Neumann
  treatment — this alone regularizes the pressure system, so (unlike
  Lid-Driven Cavity / DHC) `pinReferencePressure` is **not** called for this
  scenario.
- **Linear solver**: the near-obstacle mesh is highly stretched (high cell
  aspect ratio), which raises the pressure Laplacian's condition number well
  beyond what plain CG handles in a reasonable iteration budget — this
  scenario uses Jacobi-preconditioned CG (`PCG`) instead. Lid-Driven Cavity
  and DHC are unaffected (they still use plain `CG`, unchanged).
- **Force/Strouhal post-processing**:
  [`ObstacleForces`](../../src/postprocessing/ObstacleForces.hpp)
  surface-integrates pressure and viscous shear over the 4 obstacle faces
  into `Cd`/`Cl`, recorded every step to `forces.csv`.
  `estimate_shedding_frequency()` does peak-detection on `Cl(t)` (discarding
  the first half of the series as the startup transient) to produce a
  provisional Strouhal number in `strouhal.csv`.

## 4. Outputs

Per Reynolds case, written to `output/square_cylinder/Re_<value>/`:
- `pressure.vtk`, `u_velocity.vtk`, `v_velocity.vtk`.
- `forces.csv` — `step, time, Cd, Cl`.
- `strouhal.csv` — `Re, St` (provisional estimate).
- `square_cylinder_plot.png` (from the plotter) — vorticity field with the
  obstacle outlined, plus the `Cd(t)`/`Cl(t)` time series.

A one-off `output/obstacle_mask.vtk` (cell-data scalar, 1 = solid) lets you
confirm the obstacle footprint visually in ParaView independent of any flow
case. `output/square_cylinder/strouhal_summary.png` plots `St` vs `Re` across
the sweep against the literature reference band — explicitly annotated as
unvalidated.

## 5. Validation status (read before trusting any number here)

What has been checked (this round):
- The obstacle footprint exports correctly and sits where configured
  (`output/obstacle_mask.vtk`, checkpoint 1).
- A short run produces no NaNs/divergence at Re=100 on the shipped mesh
  (checkpoint 2) and stays numerically well-behaved over hundreds of steps
  (residuals and `Cd`/`Cl` both settle smoothly, no oscillatory blow-up).
- `Cd` converges to a plausible order of magnitude during the initial
  transient; the qualitative wake structure (recirculation behind the
  obstacle) looks as expected in the vorticity plot.

What has **not** been checked:
- Whether the shipped domain size, blockage ratio, mesh resolution, and `dt`
  actually reproduce the literature Strouhal band (~0.13–0.15 at Re~100–200)
  — this requires running the full configured step count (tens of thousands
  of steps) to capture enough shedding periods, and likely several rounds of
  mesh/`dt`/domain-length tuning. Not attempted this round by explicit scope
  decision.
- Drag/lift coefficient values against published steady-state or
  time-averaged references.
- Mesh-independence (convergence study) of any of the above.

## 6. Known limitations

- The Laplacian's FVM stencil is only first-order accurate on stretched
  meshes (see [../theory/mesh.md](../theory/mesh.md)) — this scenario's mesh
  is stretched near the obstacle, so this caveat applies here too.
- The ghost-image obstacle-surface stencils assume the obstacle boundary
  coincides exactly with mesh faces (enforced by snapping in
  `mark_solid_block`) — an obstacle size/position that doesn't align well
  with the configured mesh zones will be silently snapped to the nearest
  face, which a warning on stdout flags if the snap distance is non-trivial.
- Only a single rectangular obstacle is supported (no multiple bodies, no
  non-rectangular shapes).

## 7. See also

- [../theory/mesh.md §5](../theory/mesh.md#5-zones-the-domain-is-split-into-a-33-grid-of-regions) — the zone/stretching system the obstacle placement builds on.
- [../theory/solver.md](../theory/solver.md) — how the obstacle-aware stencils extend the existing ghost-image method.
- [lid_driven_cavity.md](lid_driven_cavity.md),
  [differentially_heated_cavity.md](differentially_heated_cavity.md) — the
  other two scenarios, both fully validated.
