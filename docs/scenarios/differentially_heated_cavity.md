# Scenario: Differentially Heated Cavity (DHC)

A square cavity with a **hot** vertical wall and a **cold** vertical wall
(adiabatic top/bottom), driving buoyancy-induced convection rolls — the
classic de Vahl Davis (1983) natural-convection benchmark.

This page is a quick-start summary. For the full line-by-line implementation
walkthrough (equations, information flow, file-by-file trace through the
code), see the in-depth companion document:
**[differentially_heated_cavity_EN.md](differentially_heated_cavity_EN.md)**.

For the underlying equations (Boussinesq buoyancy, energy equation) see
[../theory/equations.md](../theory/equations.md); for the solver classes
involved (`EnergySolver`, `Buoyancy.hpp`) see
[../theory/solver.md §5](../theory/solver.md#5-energysolver--buoyancyhpp-bolting-on-a-new-transported-scalar).

---

## 1. Physical setup

- Unit square domain.
- **West wall:** hot, $\theta = 1$ (`WALL_ISOTHERMAL`).
- **East wall:** cold, $\theta = 0$ (`WALL_ISOTHERMAL`).
- **North/South walls:** adiabatic, zero heat flux (`WALL_ADIABATIC`).
- All four velocity walls: no-slip (closed cavity, no lid).
- Governing non-dimensional numbers: Rayleigh $Ra$ (convection strength) and
  Prandtl $Pr$ (momentum/thermal diffusivity ratio), de Vahl Davis scaling.

## 2. How to run it

1. Set `solve_energy = 1` and choose `Pr` (typically `0.71`, air) in
   [input/sim_params.txt](../../input/sim_params.txt).
2. **Use a small `dt`** — see the stability warning below. The shipped
   `sim_params.txt` already has `dt = 1e-5`, validated for this scenario.
3. List the Rayleigh numbers to sweep in
   [input/rayleigh_cases.txt](../../input/rayleigh_cases.txt) (defaults to
   $Ra = 10^3, 10^4, 10^5, 10^6$ if empty).
4. Temperature boundary conditions come from
   [input/boundaries_T.txt](../../input/boundaries_T.txt) (hot west / cold
   east / adiabatic top-bottom, already configured).
5. Build and run:
   ```bash
   bash run.sh mesh_gen
   ```
6. Visualize:
   ```bash
   python3 plotters/plot_dhc.py
   ```

## 3. Outputs

Per Rayleigh case, written to `output/Ra_<value>/`:
- `pressure.vtk`, `u_velocity.vtk`, `v_velocity.vtk`, `temperature.vtk`.
- `dhc_metrics.csv` — `Ra, Pr, u_max, v_max, Nu_avg`.

## 4. Validation — de Vahl Davis (1983)

| $Ra$ | Reference $\overline{Nu}$ (hot wall) |
|---|---|
| $10^3$ | 1.118 |
| $10^4$ | 2.243 |
| $10^5$ | 4.519 |
| $10^6$ | 8.800 |

The solver matches the $Ra=10^4$ benchmark within a few percent
(measured $\overline{Nu}\approx 2.31$ on a 40² mesh, `dt=2e-5`). See
[differentially_heated_cavity_EN.md §8](differentially_heated_cavity_EN.md#8-validation-and-output)
for the Nusselt-number formula and how it's computed from the discrete
temperature field.

## 5. CRITICAL stability warning

Both the momentum and energy updates are **fully explicit**, so this
scenario needs a **much smaller** time step than the lid-driven cavity:
$$
\Delta t \lesssim \frac{\Delta x^2}{2\,\max(\alpha,\nu)}, \qquad \alpha=1,\ \ \nu = Pr
$$
On a production 150² mesh this means `dt ~ 1e-5` or smaller — using the
lid-cavity-scale `dt = 1e-3` **diverges to NaN** (the pressure CG solve
fails to converge first). See
[../theory/solver.md §6](../theory/solver.md#6-stability-this-is-a-fully-explicit-scheme)
for the general explanation.
