# NavierStokesSolver — Project Summary

## Mesh & Solver Core

A 2D incompressible Navier-Stokes solver built from scratch in C++ on a staggered MAC grid, using the **fractional-step (projection) method**.

The fractional-step method splits each timestep into stages so that the velocity and pressure fields don't have to be solved for simultaneously (as a fully coupled system would require):

1. **Predictor step — momentum equation:** solve the momentum equation *without* the pressure gradient term, using only convection and diffusion. This gives an intermediate, not-yet-divergence-free velocity field **u\*, v\*** (`MomentumSolver`, explicit Adams–Bashforth 2 with an Euler first step). Boundary conditions are enforced via a ghost-cell image method.
2. **Pressure correction step:** the intermediate field u*, v* generally does not satisfy continuity (∇·u = 0), so a pressure field is solved for that, when its gradient is subtracted from u*/v*, restores a divergence-free velocity. This is the pressure Poisson equation (PPE), assembled by `PressurePoissonSystem` and solved with a warm-started Conjugate Gradient (CG) solver (Jacobi-preconditioned CG for the Square Cylinder).
3. **Projection:** u, v = u*, v* − (dt/ρ)·∇p, giving the final divergence-free velocity field for the timestep.

**Supporting infrastructure:**
- `MeshConfig` (`src/mesh/MeshGenerator.cpp`) builds the staggered grid.
- Gradient / divergence / Laplacian / diffusion operators are factored into their own classes (`src/operators/`).
- `operator_verify` executable — checks the operators independently (10 MMS suites).
- `reynolds_sweep` executable — the same driver as `mesh_gen` under a separate name; both sweep the active scenario's Re/Ra list.

**Key bugs resolved during development:**
- Sign conventions in the PPE RHS and matrix assembly.
- Pinning the reference pressure (zeroing row *and* column, not just diagonal).
- Ghost-cell wall distance needing a factor-of-2 correction.
- Square Cylinder: swapped viscous force components, obstacle-surface flux leaking into the pressure equation, inverted wall distances at the obstacle, and a NaN mistaken for convergence (see [docs/scenarios/square_cylinder_validation.md](docs/scenarios/square_cylinder_validation.md)).

The core layer is stable and validated.

## Scenario selection

The active scenario is chosen in `input/scenario.txt` (`lid_driven_cavity`, `differentially_heated_cavity`, `square_cylinder` or `square_cylinder_breuer`). Every other input file is then read from `input/<scenario>/`.

## Scenario 1 — Lid-Driven Cavity

The classic benchmark: top wall moving at U = 1, other walls no-slip.

- Validated against **Ghia et al. (1982)** at Re = 100, 400, 1000 and 3200 on a 120×120 stretched grid: ψ_min within 0.05–0.3 % up to Re = 1000 and 2.5 % at Re = 3200.
- **Result:** a primary recirculating vortex plus secondary corner eddies; the u/v centreline velocity profiles overlay the Ghia reference data.

## Scenario 2 — Differentially Heated Cavity (DHC)

The core solver above (fractional step, momentum predictor u*/v* + pressure correction) solves the **momentum equation only**. For the Differentially Heated Cavity, it is extended with the **energy equation** on top, without modifying the momentum/pressure core:

- `EnergySolver` — advection-diffusion of temperature θ, solved as an additional step each timestep.
- `Buoyancy.hpp` — a source term computed from θ and injected into the momentum predictor step (between computing u*/v* and the pressure correction), coupling temperature back into the flow.

Non-dimensionalized with Rayleigh/Prandtl numbers, following **de Vahl Davis (1983)** scaling. Hot/cold side walls, adiabatic top/bottom.

- **Result:** buoyancy-driven convection rolls, and an average Nusselt number (Nu) characterizing heat transfer.
- Validated at Ra = 10³–10⁶: Nu within 1.2 / 2.2 / 4.5 / 2.9 % of the reference (e.g. Nu = 2.292 vs 2.243 at Ra = 10⁴).
- **Stability caveat:** the scheme is fully explicit, so DHC requires a much smaller timestep (dt ~ 1e-5) than the lid-driven case to remain stable — this is a known/expected constraint, not a bug.

The energy equation runs only when `scenario = differentially_heated_cavity` **and** `solve_energy = 1` in `input/differentially_heated_cavity/sim_params.txt`.

## Scenario 3 — Square Cylinder

A square obstacle immersed in a channel (inlet/outlet), with drag, lift and Strouhal post-processing (`src/postprocessing/ObstacleForces.*`).

- Validated against **Breuer et al. (2000)** (blockage 1/8, parabolic inflow): C_D +0.5 % and recirculation length −0.4 % at Re = 40; shedding onset between Re = 50 and 70 (reference ≈ 60); Strouhal number +1.4 % at Re = 100.

Full results and references are in the [README](README.md#results) and the per-scenario pages in [docs/scenarios/](docs/scenarios/).
