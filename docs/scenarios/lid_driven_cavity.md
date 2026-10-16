# Scenario: Lid-Driven Cavity

The classic incompressible-flow benchmark: a square cavity with three
stationary no-slip walls and a lid moving at constant horizontal velocity,
driving a primary recirculating vortex plus secondary corner eddies.

For the governing equations and fractional-step method, see
[../theory/equations.md](../theory/equations.md); for the mesh and solver
classes involved, see [../theory/mesh.md](../theory/mesh.md) and
[../theory/solver.md](../theory/solver.md).

---

## 1. Physical setup

- Unit square domain, $L = H = 1$.
- **North wall (lid):** $u = U_{lid}$ (typically 1), $v = 0$ — `WALL_FIXED_VALUE`.
- **South / West / East walls:** no-slip, $u = v = 0$ — `WALL`.
- Governing non-dimensional number: Reynolds number
  $Re = U_{lid}\,L/\nu$.
- Only the momentum + pressure equations are solved — no energy equation,
  no buoyancy (`solve_energy = 0`).

## 2. How to run it

1. Set `scenario = lid_driven_cavity` in [input/scenario.txt](../../input/scenario.txt).
2. `solve_energy = 0` in
   [input/lid_driven_cavity/sim_params.txt](../../input/lid_driven_cavity/sim_params.txt)
   (already set).
3. Choose which Reynolds numbers to sweep in
   [input/lid_driven_cavity/reynolds_cases.txt](../../input/lid_driven_cavity/reynolds_cases.txt)
   (one or more integers, space/comma separated). If the file is
   empty/missing, the solver falls back to $Re = 100, 200, \dots, 1000$.
4. Boundary conditions are read from
   [input/lid_driven_cavity/boundaries.txt](../../input/lid_driven_cavity/boundaries.txt)
   — it already ships configured for this scenario (moving lid on `NORTH`,
   no-slip elsewhere).
5. Build and run:
   ```bash
   bash run.sh mesh_gen
   ```
5. Visualize and validate:
   ```bash
   python3 plotters/plot_results.py
   ```

## 3. What the solver computes

For each Reynolds number in the sweep:
- `nu` is derived from `Re` as $\nu = U_{lid}L/Re$
  ([src/main.cpp:302](../../src/main.cpp#L302)).
- The time loop runs the standard predictor → pressure-projection pipeline
  (see [../theory/solver.md §0](../theory/solver.md#0-overview-the-per-timestep-pipeline))
  with no buoyancy/energy steps.
- Steady state is detected when the max velocity change per step,
  $\lVert \Delta u\rVert_\infty/\Delta t$, drops below `ss_tol`
  ([src/main.cpp:392](../../src/main.cpp#L392)).

## 4. Outputs

Per case, written to `output/lid_driven_cavity/Re_<value>/`:
- `fields/pressure.vtk`, `fields/u_velocity.vtk`, `fields/v_velocity.vtk` —
  rectilinear VTK fields, open directly in ParaView.
- `data/centerline_u.csv` — $u(y)$ along the vertical mid-plane $x=L/2$.
- `data/centerline_v.csv` — $v(x)$ along the horizontal mid-plane $y=H/2$.
- `figures/` — written by `plotters/plot_results.py` (see §5).

The meshes go to `output/lid_driven_cavity/mesh/`.

## 5. Validation — Ghia, Ghia & Shin (1982)

The solver's centerline velocity profiles are compared against the
widely-used reference data from Ghia et al. (1982)
([input/lid_driven_cavity/ghia_reference.txt](../../input/lid_driven_cavity/ghia_reference.txt),
also embedded directly in [src/main.cpp](../../src/main.cpp) for the console
printout), available at $Re = 100, 400, 1000, 3200, 5000, 7500, 10000$.

| $Re$ | Expected flow pattern |
|---|---|
| 100 | single primary vortex, tiny corner eddies |
| 400–1000 | growing corner eddies, vortex center shifts |
| 3200+ | tertiary corner eddies appear, resolution-sensitive |

Run `python3 plotters/plot_results.py`. Per case it writes
`figures/flow_fields` (stream function with Ghia's contour levels, vorticity,
pressure) and `figures/centerlines` (profiles vs. Ghia, only when the case's
$Re$ is tabulated there). Across cases it writes
`output/lid_driven_cavity/figures/{centerline_comparison,centerline_trend,streamfunction_overview}`
and `output/lid_driven_cavity/validation_summary.csv` (primary-vortex
$\psi_{\min}$ and centre, RMS/max centreline errors vs. Ghia). Every figure is
saved as a vector PDF (for LaTeX) and a 300 dpi PNG.

![Centreline profiles vs. Ghia et al.](../figures/lid_centerline_comparison.png)

## 6. Known caveats

- **Suspect reference value, Re = 400, $v$ table, $x/L = 0.9063$.** The
  embedded table has $-0.23827$; the solution ($-0.3816$) matches the other
  14 points of the column to 0.001 and the tabulated value breaks the
  monotonic trend of its own column. Likely a transcription error — check
  against Ghia et al. (1982), Table II, before quoting the Re = 400 $v$
  error (0.037 with it, ≈ 0.003 without).

- At higher $Re$ (≥3200), resolving the secondary/tertiary corner eddies
  needs a correspondingly finer mesh (`stretch.txt`) — a coarse uniform grid
  will visibly disagree with Ghia near the corners while still matching the
  bulk primary vortex.
- CDS convection (the default) can show mild oscillations at high cell
  Peclet numbers on coarse grids; switch to `ConvectionScheme::UPWIND` in
  `main.cpp`'s `MomentumSolver` construction if this becomes an issue (not
  exposed as an input-file parameter — see
  [README "Future improvements"](../../README.md) for the related UX note).
