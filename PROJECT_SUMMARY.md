# NavierStokesSolver — Project Summary

## Mesh & Solver Core

A 2D incompressible Navier-Stokes solver built from scratch in C++ on a staggered MAC grid, using the **fractional-step (projection) method**.

The fractional-step method splits each timestep into two stages so that the velocity and pressure fields don't have to be solved for simultaneously (as a fully coupled system would require):

1. **Predictor step — momentum equation:** solve the momentum equation *without* the pressure gradient term, using only convection and diffusion. This gives an intermediate, not-yet-divergence-free velocity field **u\*, v\*** (`MomentumSolver`, explicit Euler). Boundary conditions are enforced via a ghost-cell image method.
2. **Pressure correction step:** the intermediate field u*, v* generally does not satisfy continuity (∇·u = 0), so a pressure field is solved for that, when its gradient is subtracted from u*/v*, restores a divergence-free velocity. This is the pressure Poisson equation (PPE), assembled and solved by `PressurePoissonSystem` with a warm-started Conjugate Gradient (CG) solver.
3. **Projection:** u, v = u*, v* − (dt/ρ)·∇p, giving the final divergence-free velocity field for the timestep.

**Supporting infrastructure:**
- `MeshGenerator` builds the staggered grid.
- Gradient / divergence / Laplacian operators are factored into their own classes.
- `operator_verify` executable — checks the operators independently.
- `reynolds_sweep` executable — runs batches of simulations across parameter values (Re/Ra).

**Supporting infrastructure:**
- Gradient / divergence / Laplacian operators are factored into their own classes.
- `operator_verify` executable — checks the operators independently.
- `reynolds_sweep` executable — runs batches of simulations across parameter values (Re/Ra).

**Key bugs resolved during development:**
- Sign conventions in the PPE RHS and matrix assembly.
- Pinning the reference pressure (zeroing row *and* column, not just diagonal).
- Ghost-cell wall distance needing a factor-of-2 correction.

The core layer is now stable and validated.

## Scenario 1 — Lid-Driven Cavity

The classic benchmark: top wall moving at U = 1, other walls no-slip.

- Validated directly against **Ghia et al. (1982)** reference velocity profiles at Re = 100 / 400.
- **Expected result:** a primary recirculating vortex plus secondary corner eddies; u/v centerline velocity profiles should closely overlay the Ghia reference data.

## Scenario 2 — Differentially Heated Cavity (DHC)

The core solver above (fractional-step, momentum predictor u*/v* + pressure correction) solves the **momentum equation only**. For the Differentially Heated Cavity, we extend it by adding the **energy equation** on top, without modifying the momentum/pressure core:

- `EnergySolver` — advection-diffusion of temperature θ, solved as an additional step each timestep.
- `Buoyancy` — a source term computed from θ and injected into the momentum predictor step (between computing u*/v* and the pressure correction), coupling temperature back into the flow.

Non-dimensionalized with Rayleigh/Prandtl numbers, following **de Vahl Davis (1983)** scaling. Hot/cold side walls, adiabatic top/bottom.

- **Expected result:** buoyancy-driven convection rolls, and a Nusselt number (Nu) output characterizing heat transfer.
- Validated at Ra = 1e4: Nu ≈ 2.31, matching literature within a few percent.
- **Stability caveat:** the scheme is fully explicit, so DHC requires a much smaller timestep (dt ~ 1e-5) than the lid-driven case to remain stable — this is a known/expected constraint, not a bug.

Both scenarios are toggled via `solve_energy` (0 = Lid-Driven, 1 = DHC) in `input/sim_params.txt`.
