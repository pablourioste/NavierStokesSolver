# Governing Equations & the Fractional-Step Method

Theory notes on the physics this solver implements and the time-integration
method used to advance the discrete equations. For the spatial grid, see
[mesh.md](mesh.md); for how these equations map onto C++ classes, see
[solver.md](solver.md). For how to run each physical case, see
[../scenarios/](../scenarios/).

---

## 1. Governing equations

### 1.1 Continuity (mass conservation)

- The fluid is incompressible: $\nabla \cdot \mathbf{u} = 0$ everywhere, at
  every instant.
- This is not integrated in time — it is a **constraint**, enforced at every
  step by the pressure projection (§3).

### 1.2 Momentum (incompressible Navier–Stokes)

$$
\frac{\partial \mathbf{u}}{\partial t} + (\mathbf{u}\cdot\nabla)\mathbf{u}
= -\nabla p + \nu\,\nabla^2 \mathbf{u} + \mathbf{f}_b
$$

- $\mathbf{u} = (u, v)$: velocity field.
- $p$: pressure (kinematic, i.e. already divided by density in this
  formulation — see §1.4).
- $\nu$: kinematic viscosity.
- $\mathbf{f}_b$: body-force term, zero for the lid-driven cavity, and equal
  to the Boussinesq buoyancy term for the differentially heated cavity (§1.3).

### 1.3 Energy equation (only when buoyancy/heat transfer is modeled)

$$
\frac{\partial \theta}{\partial t} + (\mathbf{u}\cdot\nabla)\theta
= \alpha\,\nabla^2 \theta
$$

- $\theta$: non-dimensional temperature, transported like any scalar by the
  (already divergence-free) velocity field.
- $\alpha$: thermal diffusivity (= 1 in the Ra–Pr scaling used here, see §1.5).
- This equation is **decoupled from the core solver** and only switched on
  for the Differentially Heated Cavity scenario — see
  [solver.md §5](solver.md#5-energysolver--buoyancyhpp-bolting-on-a-new-transported-scalar).

### 1.4 Boussinesq buoyancy coupling

- Under the Boussinesq approximation, density variations are neglected
  everywhere **except** in a buoyancy source term added to the vertical
  momentum equation:
  $$
  \mathbf{f}_b = (0,\; Ra\,Pr\,\theta)
  $$
- This couples temperature back into the flow: hotter fluid is lighter and
  rises, colder fluid sinks — the mechanism that drives convection rolls in
  the differentially heated cavity.
- Implementation: [src/Buoyancy.hpp](../../src/Buoyancy.hpp) injects this term
  directly into the momentum predictor (see §2 below) — no other file needed
  to change to add it.

### 1.5 Non-dimensionalization

Two different non-dimensional scalings are used depending on the scenario
(both are already built into `nu` and the forcing term, not applied manually):

| Scenario | Governing number | How it enters |
|---|---|---|
| Lid-Driven Cavity | Reynolds number $Re = U_{lid} L / \nu$ | `nu` is set directly from `Re` in [sim_params.txt](../../input/sim_params.txt) / `reynolds_cases.txt` sweep |
| Differentially Heated Cavity | Rayleigh $Ra$, Prandtl $Pr$ (de Vahl Davis 1983 scaling) | momentum diffusion `nu = Pr`; buoyancy coefficient `dt·Ra·Pr` — see [Buoyancy.hpp:34](../../src/Buoyancy.hpp#L34) |

In the Ra–Pr scaling, thermal diffusivity is fixed at $\alpha = 1$, and all
of the Rayleigh-number dependence shows up solely in the buoyancy term — the
energy equation itself is unchanged by $Ra$.

---

## 2. The fractional-step (projection) method

Solving the coupled velocity-pressure system directly (fully implicit, all
unknowns at once) is expensive. **Chorin's fractional-step method** splits
each time step into two decoupled sub-steps so pressure and velocity are
never solved for simultaneously:

### Step 1 — Predictor (momentum without pressure)

Advance the momentum equation **ignoring the pressure gradient**:

$$
\mathbf{u}^* = \mathbf{u}^n + \Delta t \left[\nu\nabla^2\mathbf{u}^n -
(\mathbf{u}^n\cdot\nabla)\mathbf{u}^n + \mathbf{f}_b\right]
$$

- Produces an intermediate velocity field $\mathbf{u}^* = (u^*, v^*)$ that
  satisfies the boundary conditions but is **not** divergence-free in general.
- Implemented by `MomentumSolver::advance` — see
  [solver.md §2](solver.md#2-momentumsolver-the-predictor-step).
- For the Differentially Heated Cavity, the buoyancy term $\mathbf{f}_b$ is
  injected **after** this predictor step, directly into $v^*$ (§2 continued
  below) — the momentum solver itself never needs to know buoyancy exists.

### Step 2 — Pressure Poisson equation (PPE)

We want the **final** velocity to be divergence-free:
$$
\mathbf{u}^{n+1} = \mathbf{u}^* - \frac{\Delta t}{\rho}\nabla p
\qquad\text{such that}\qquad \nabla\cdot\mathbf{u}^{n+1} = 0
$$

Taking the divergence of the correction formula and imposing
$\nabla\cdot\mathbf{u}^{n+1}=0$ gives a Poisson equation for pressure:

$$
\nabla^2 p = \frac{\rho}{\Delta t}\,\nabla\cdot\mathbf{u}^*
$$

- This is an **elliptic constraint**, not a time-evolution equation — there
  is no $\partial p/\partial t$. It is solved anew (via Conjugate Gradient)
  at every time step from the current $\mathbf{u}^*$.
- Assembled and solved by `PressurePoissonSystem` — see
  [solver.md §3](solver.md#3-pressurepoissonsystem-the-elliptic-constraint).
- With homogeneous Neumann boundary conditions everywhere (closed cavity, no
  pressure outlet), the discrete system is singular (constant pressure is in
  the null space) and must be regularized by pinning one reference value —
  see [solver.md §3](solver.md#3-pressurepoissonsystem-the-elliptic-constraint).

### Step 3 — Projection (correction)

$$
\mathbf{u}^{n+1} = \mathbf{u}^* - \frac{\Delta t}{\rho}\nabla p
$$

- Subtracting the pressure gradient from the predictor field restores
  $\nabla\cdot\mathbf{u}=0$ exactly (up to the linear solver's tolerance).
- Boundary velocity nodes are **not** touched by this correction — they keep
  their prescribed values.

### Why split it this way?

- The predictor step is **fully explicit and local** (no linear solve needed,
  cheap per cell).
- The only implicit/global solve per time step is the pressure Poisson
  equation, which is linear, symmetric positive-definite, and well suited to
  Conjugate Gradient — far cheaper than solving the full coupled
  Navier–Stokes system every step.

---

## 3. Time integration

| Field | Equation type | Time advance |
|---|---|---|
| $u, v$ (momentum) | evolution (parabolic) | explicit Adams–Bashforth 2 (AB2), weights $3/2, -1/2$ |
| $\theta$ (energy, DHC only) | evolution (parabolic) | explicit AB2, same weights as momentum |
| $p$ (pressure) | constraint (elliptic) | **not** time-integrated — solved fresh via CG each step |

- **Adams–Bashforth 2 (AB2):**
  $$
  \phi^{n+1} = \phi^n + \Delta t\left[\tfrac{3}{2}R(\phi^n) -
  \tfrac12 R(\phi^{n-1})\right]
  $$
  where $R$ is the right-hand side (diffusion − convection, or diffusion −
  convection + buoyancy). Second-order accurate in time, using only
  already-computed data (no implicit solve).
- **First step exception:** at $n=0$ there is no $R(\phi^{-1})$ yet, so the
  very first step falls back to explicit Euler, $\phi^1 = \phi^0 + \Delta t\,
  R(\phi^0)$.
- **Why pressure has no AB2 weights:** $p$ is not advanced from its previous
  value at all — each step re-solves $\nabla^2 p = (\rho/\Delta t)\nabla\cdot
  \mathbf{u}^*$ from scratch (though warm-started from the previous solution
  for performance, see [solver.md §4](solver.md#4-the-linear-solver-warm-started-conjugate-gradient)).
  It is an algebraic constraint on the current state, not a field with its
  own dynamics.

---

## 4. Boundary conditions — the ghost-cell (image) method

All solid-wall and prescribed-value boundary conditions are enforced the same
way throughout the code (momentum, energy): by defining a **ghost value**
just outside the domain so that the standard interior finite-difference
stencil, applied at the last interior node, automatically reproduces the
correct wall condition.

- **Dirichlet (fixed value $\phi_{wall}$)**, e.g. no-slip velocity or an
  isothermal wall:
  $$
  \phi_{ghost} = 2\,\phi_{wall} - \phi_P
  $$
  so that the average of the ghost and the last interior node
  $(\phi_{ghost}+\phi_P)/2$ equals $\phi_{wall}$ exactly.
- **Neumann (zero flux, adiabatic)**:
  $$
  \phi_{ghost} = \phi_P
  $$
  i.e. zero normal gradient across the wall.
- The ghost values are computed **once per case**, from the boundary
  condition tables (`input/boundaries.txt`, `input/boundaries_T.txt`), and
  reused at every time step — not recomputed per iteration.
- Geometric details (why the wall-distance needs a factor of 2, how this
  interacts with the non-uniform/stretched grid) are covered in
  [mesh.md](mesh.md).

---

## Summary: one umbrella, two scenarios

| | Lid-Driven Cavity | Differentially Heated Cavity |
|---|---|---|
| Equations solved | continuity + momentum | continuity + momentum (+ buoyancy) + energy |
| Driving mechanism | moving lid (Dirichlet velocity BC) | horizontal temperature gradient → buoyancy |
| Non-dimensional number | $Re$ | $Ra$, $Pr$ |
| Core solver touched? | — | **no** — energy/buoyancy are added externally |
| Switch | `solve_energy = 0` in `input/sim_params.txt` | `solve_energy = 1` |

See [../scenarios/lid_driven_cavity.md](../scenarios/lid_driven_cavity.md) and
[../scenarios/differentially_heated_cavity.md](../scenarios/differentially_heated_cavity.md)
for how to run each one.
