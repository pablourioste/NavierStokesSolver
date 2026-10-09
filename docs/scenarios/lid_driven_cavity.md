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

1. Set `solve_energy = 0` in [input/sim_params.txt](../../input/sim_params.txt).
2. Choose which Reynolds numbers to sweep in
   [input/reynolds_cases.txt](../../input/reynolds_cases.txt) (one or more
   integers, space/comma separated). If the file is empty/missing, the
   solver falls back to $Re = 100, 200, \dots, 1000$.
3. Boundary conditions are read from
   [input/boundaries.txt](../../input/boundaries.txt) — it already ships
   configured for this scenario (moving lid on `NORTH`, no-slip elsewhere).
   [input/boundaries_lid_cavity.txt](../../input/boundaries_lid_cavity.txt) is
   an equivalent, more heavily commented reference copy of the same file.
4. Build and run:
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
  ([src/main.cpp:191](../../src/main.cpp#L191)).
- The time loop runs the standard predictor → pressure-projection pipeline
  (see [../theory/solver.md §0](../theory/solver.md#0-overview-the-per-timestep-pipeline))
  with no buoyancy/energy steps.
- Steady state is detected when the max velocity change per step,
  $\lVert \Delta u\rVert_\infty/\Delta t$, drops below `ss_tol`
  ([src/main.cpp:252](../../src/main.cpp#L252)).

## 4. Outputs

Per case, written to `output/Re_<value>/`:
- `pressure.vtk`, `u_velocity.vtk`, `v_velocity.vtk` — rectilinear VTK fields,
  open directly in ParaView.
- `centerline_u.csv` — $u(y)$ along the vertical mid-plane $x=L/2$.
- `centerline_v.csv` — $v(x)$ along the horizontal mid-plane $y=H/2$.

## 5. Validation — Ghia, Ghia & Shin (1982)

The solver's centerline velocity profiles are compared against the
widely-used reference data from Ghia et al. (1982)
([input/ghia_reference.txt](../../input/ghia_reference.txt), also embedded
directly in [src/main.cpp](../../src/main.cpp) for the console printout),
available at $Re = 100, 400, 1000, 3200, 5000, 7500, 10000$.

| $Re$ | Expected flow pattern |
|---|---|
| 100 | single primary vortex, tiny corner eddies |
| 400–1000 | growing corner eddies, vortex center shifts |
| 3200+ | tertiary corner eddies appear, resolution-sensitive |

Run `python3 plotters/plot_results.py` to overlay the computed $u(y)$/$v(x)$
profiles on the Ghia reference points — see
`no_Stretched/validation.png` and `no_Stretched/validation_all_re_u.png`
for pre-computed examples at multiple Reynolds numbers.

## 6. Known caveats

- At higher $Re$ (≥3200), resolving the secondary/tertiary corner eddies
  needs a correspondingly finer mesh (`stretch.txt`) — a coarse uniform grid
  will visibly disagree with Ghia near the corners while still matching the
  bulk primary vortex.
- CDS convection (the default) can show mild oscillations at high cell
  Peclet numbers on coarse grids; switch to `ConvectionScheme::UPWIND` in
  `main.cpp`'s `MomentumSolver` construction if this becomes an issue (not
  exposed as an input-file parameter — see
  [README "Future improvements"](../../README.md) for the related UX note).
