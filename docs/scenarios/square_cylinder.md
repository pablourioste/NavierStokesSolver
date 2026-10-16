# Scenario: Square Cylinder Problem

> **Status: validated against Breuer et al. (2000)** in the dedicated
> `square_cylinder_breuer` scenario (blockage 1/8, parabolic inflow):
> $C_D$ and $L_r$ at $Re = 40$ within 0.5 %, onset between $Re = 50$ and
> 70, $St$ within 1.4 % at $Re = 100$ (grid-converged). See
> [square_cylinder_validation.md](square_cylinder_validation.md) for the
> references, the four solver bugs fixed on the way, and what remains. The
> shipped `square_cylinder` scenario (uniform inflow between no-slip walls)
> matches no published benchmark.

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
- `fields/pressure.vtk`, `fields/u_velocity.vtk`, `fields/v_velocity.vtk`.
- `data/forces.csv` — `step, time, Cd, Cl`.
- `data/strouhal.csv` — `Re, St` (provisional estimate).
- `figures/vorticity`, `figures/force_coefficients`, `figures/lift_spectrum`
  (from the plotter) — vorticity field with the obstacle drawn in,
  $C_D(t^*)$ / $C_L(t^*)$, and the amplitude spectrum of $C_L$ over the second
  half of the record (its peak gives an FFT-based $St$, independent of the
  solver's zero-crossing estimate).

A one-off `output/square_cylinder/mesh/obstacle_mask.vtk` (cell-data scalar,
1 = solid) lets you confirm the obstacle footprint visually in ParaView
independent of any flow case. `output/square_cylinder/figures/strouhal` plots
`St` vs `Re` across the sweep against Sohankar et al. (1998, $\beta = 5\,\%$) —
explicitly annotated as unvalidated; it is only drawn when at least two
Reynolds numbers have been run — and `force_summary.csv` lists mean $C_D$,
rms $C_L$, `St` (solver) and `St_fft` per case.

## 5. Validation status (read before trusting any number here)

What has been checked (this round):
- The obstacle footprint exports correctly and sits where configured
  (`output/square_cylinder/mesh/obstacle_mask.vtk`, checkpoint 1).
- A short run produces no NaNs/divergence at Re=100 on the shipped mesh
  (checkpoint 2) and stays numerically well-behaved over hundreds of steps
  (residuals and `Cd`/`Cl` both settle smoothly, no oscillatory blow-up).
- `Cd` converges to a plausible order of magnitude during the initial
  transient; the qualitative wake structure (recirculation behind the
  obstacle) looks as expected in the vorticity plot.

Shipped scenario, $Re = 150$, coarse 126×88 grid, `dt = 5e-4`,
$tU_\infty/D = 200$, re-run after the solver fixes (results in
`output/square_cylinder/Re_150/`): periodic shedding from
$tU_\infty/D \approx 130$, $St = 0.196$, $\overline{C_D} = 1.84$,
$C_{L,\mathrm{rms}} = 0.30$, $\overline{C_L} \approx 0$. Not comparable with a
published value (uniform inflow between no-slip walls); the earlier numbers
in this section ($\overline{C_D} = 2.57$, $\overline{C_L} \approx -0.4$) came
from the bugs listed in [square_cylinder_validation.md §4](square_cylinder_validation.md)
and are void.

What has **not** been checked:
- Whether the shipped domain size, blockage ratio, mesh resolution, and `dt`
  actually reproduce the reference values in
  [square_cylinder_validation.md](square_cylinder_validation.md)
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
