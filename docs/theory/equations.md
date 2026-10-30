# Governing Equations and the Fractional-Step Method

The physics the solver implements, the non-dimensional form used by each
scenario, and the time-integration method. For the grid see
[mesh.md](mesh.md); for the boundary conditions see
[boundary_conditions.md](boundary_conditions.md); for how the equations map
onto C++ classes see [solver.md](solver.md).

---

## Contents

1. [Governing equations](#1-governing-equations)
2. [Non-dimensional form used by each scenario](#2-non-dimensional-form-used-by-each-scenario)
3. [The fractional-step (projection) method](#3-the-fractional-step-projection-method)
4. [Time integration](#4-time-integration)
5. [Spatial discretisation (summary)](#5-spatial-discretisation-summary)
6. [Stability of the explicit scheme](#6-stability-of-the-explicit-scheme)
7. [Derived quantities](#7-derived-quantities-used-for-validation)

---

## 1. Governing equations

Two-dimensional, incompressible, constant-property flow.

**Continuity**

$$
\nabla\cdot\mathbf{u} = 0
$$

It is not integrated in time: it is a **constraint**, enforced every step by
the pressure projection (§3).

**Momentum (Navier–Stokes)**

$$
\frac{\partial \mathbf{u}}{\partial t} + (\mathbf{u}\cdot\nabla)\mathbf{u}
= -\frac{1}{\rho}\nabla p + \nu\,\nabla^2 \mathbf{u} + \mathbf{f}_b
$$

with $\mathbf{u} = (u, v)$, pressure $p$, density $\rho$ (= 1 in every
input), kinematic viscosity $\nu$ and body force $\mathbf{f}_b$ — zero except
for the buoyancy of the heated cavity.

**Energy (differentially heated cavity only)**

$$
\frac{\partial \theta}{\partial t} + (\mathbf{u}\cdot\nabla)\theta
= \alpha\,\nabla^2 \theta
$$

with the non-dimensional temperature $\theta = (T - T_c)/(T_h - T_c)$ and
thermal diffusivity $\alpha$. Under the **Boussinesq approximation** density
is constant except in the buoyancy force, which in dimensional form is
$\mathbf{f}_b = g\beta_T (T - T_{ref})\,\hat{\mathbf{e}}_y$. The constant
part of the buoyancy is absorbed by the pressure, so the code applies the
force with $\theta$ directly.

---

## 2. Non-dimensional form used by each scenario

All inputs are already non-dimensional ($\rho = 1$, reference length and
velocity equal to 1). The swept number (Re or Ra) sets the coefficients
below; the code derives them per case in [main.cpp:302-312](../../src/main.cpp#L302-L312).

### 2.1 Lid-driven cavity

Scales: cavity side $L$, lid velocity $U_{lid}$, time $L/U_{lid}$, pressure
$\rho U_{lid}^2$.

$$
\nabla\cdot\mathbf{u}=0,\qquad
\frac{\partial \mathbf{u}}{\partial t} + (\mathbf{u}\cdot\nabla)\mathbf{u}
= -\nabla p + \frac{1}{Re}\nabla^2\mathbf{u},
\qquad Re = \frac{U_{lid}L}{\nu}
$$

Code: `nu = U_lid * L_domain / Re`.

### 2.2 Differentially heated cavity (de Vahl Davis 1983 scaling)

Scales: cavity side $L$, velocity $\alpha/L$, time $L^2/\alpha$, pressure
$\rho\alpha^2/L^2$.

$$
\nabla\cdot\mathbf{u}=0,\qquad
\frac{\partial \mathbf{u}}{\partial t} + (\mathbf{u}\cdot\nabla)\mathbf{u}
= -\nabla p + Pr\,\nabla^2\mathbf{u} + Ra\,Pr\,\theta\,\hat{\mathbf{e}}_y,
\qquad
\frac{\partial \theta}{\partial t} + (\mathbf{u}\cdot\nabla)\theta = \nabla^2\theta
$$

$$
Ra = \frac{g\beta_T\,\Delta T\,L^3}{\nu\,\alpha},\qquad Pr = \frac{\nu}{\alpha}
$$

Code: `nu = Pr`, `alpha = 1`, buoyancy increment `dt * Ra * Pr * theta`
added to $v^*$ ([Buoyancy.hpp](../../src/Buoyancy.hpp)). The Rayleigh number
appears **only** in the buoyancy term. Velocities come out in units of
$\alpha/L$, the units in which de Vahl Davis tabulates $u_{\max}$ and
$v_{\max}$.

### 2.3 Square cylinder (both variants)

Scales: cylinder side $D$, inflow velocity $U$ (uniform $U_\infty$, or the
centreline velocity $U_{\max}$ of the parabolic inflow in the Breuer set-up).

$$
\nabla\cdot\mathbf{u}=0,\qquad
\frac{\partial \mathbf{u}}{\partial t} + (\mathbf{u}\cdot\nabla)\mathbf{u}
= -\nabla p + \frac{1}{Re}\nabla^2\mathbf{u},
\qquad Re = \frac{U D}{\nu}
$$

Code: `nu = U_inf * D / Re`, with $D$ = `obstacle_xmax - obstacle_xmin`
after snapping. Lengths in the input files are not divided by $D$
($D = 0.125$, $H = 1$), so the convective time unit is
$t^* = tU/D = 8\,t$ for $D = 0.125$, $U = 1$.

### 2.4 Summary

| | Lid-driven cavity | Differentially heated cavity | Square cylinder |
|---|---|---|---|
| Equations | continuity + momentum | continuity + momentum + buoyancy + energy | continuity + momentum |
| Driving mechanism | moving lid | temperature difference → buoyancy | inflow past an obstacle |
| Swept number | $Re$ | $Ra$ (fixed $Pr$) | $Re$ |
| Viscosity in the code | $U_{lid}L/Re$ | $Pr$ | $U D/Re$ |
| Steady or unsteady | steady | steady | steady ($Re \lesssim 50$) / periodic |

---

## 3. The fractional-step (projection) method

Velocity and pressure are not solved together. Chorin's projection method
splits every time step into three sub-steps.

**Step 1 — Predictor.** Advance momentum without the pressure gradient:

$$
\mathbf{u}^* = \mathbf{u}^n + \Delta t\left[\tfrac32 \mathbf{R}^n - \tfrac12
\mathbf{R}^{n-1}\right],\qquad
\mathbf{R} = \nu\nabla^2\mathbf{u} - (\mathbf{u}\cdot\nabla)\mathbf{u}
$$

($\mathbf{R}$ evaluated on the current velocity; AB2, §4). For the heated
cavity the buoyancy increment $\Delta t\,Ra\,Pr\,\theta^n$ is then added to
$v^*$. $\mathbf{u}^*$ satisfies the boundary conditions but is not
divergence-free.

**Step 2 — Pressure Poisson equation (PPE).** Requiring
$\nabla\cdot\mathbf{u}^{n+1} = 0$ in Step 3 gives

$$
\nabla^2 p^{n+1} = \frac{\rho}{\Delta t}\,\nabla\cdot\mathbf{u}^*
$$

an elliptic problem solved from scratch every step (warm-started CG). The
code assembles the positive-definite form $A\,p = b$ with
$A = -\nabla^2_h$ and $b = -(\rho/\Delta t)\,(\nabla\cdot\mathbf{u}^*)\,V$
per cell (finite-volume integral). The pressure conditions
(Neumann at walls, Dirichlet at an outlet, reference pinning in closed
cavities) are in [boundary_conditions.md §5](boundary_conditions.md#5-pressure).

**Step 3 — Projection.**

$$
\mathbf{u}^{n+1} = \mathbf{u}^* - \frac{\Delta t}{\rho}\nabla p^{n+1}
$$

applied at interior velocity nodes only; boundary nodes keep their
prescribed values. The discrete divergence of $\mathbf{u}^{n+1}$ is zero to
the tolerance of the linear solver, because the discrete Laplacian is
exactly the discrete divergence of the discrete gradient on the staggered
grid ([../reference/operators_verification.md](../reference/operators_verification.md),
Suite 5).

**Step 4 (heated cavity only) — Energy.** Advance $\theta$ with the
divergence-free $\mathbf{u}^{n+1}$ (AB2).

This is the *non-incremental* projection: the pressure of the previous step
is not used in the predictor. Its splitting error makes steady results depend
weakly on $\Delta t$ (≈ 1 % in the cylinder drag when $\Delta t$ is halved,
see [../scenarios/square_cylinder.md](../scenarios/square_cylinder.md)).

---

## 4. Time integration

| Field | Equation type | Time advance | Code |
|---|---|---|---|
| $u, v$ | parabolic (evolution) | explicit Adams–Bashforth 2 | [`MomentumSolver::advance`](../../src/MomentumSolver.cpp#L359) |
| $\theta$ | parabolic (evolution) | explicit Adams–Bashforth 2 | [`EnergySolver::advance`](../../src/EnergySolver.cpp#L118) |
| $p$ | elliptic (constraint) | not time-integrated, solved each step | `PressurePoissonSystem` + `CG`/`PCG` |

$$
\phi^{n+1} = \phi^n + \Delta t\left[\tfrac32 R(\phi^n) - \tfrac12 R(\phi^{n-1})\right]
$$

Second-order in time; the first step uses explicit Euler, since
$R(\phi^{-1})$ does not exist. The whole right-hand side is computed from
$\phi^n$ before any node is updated.

---

## 5. Spatial discretisation (summary)

| Term | Scheme | Order (uniform grid) |
|---|---|---|
| Diffusion $\nu\nabla^2\mathbf{u}$, $\alpha\nabla^2\theta$ | central differences on the non-uniform grid, image ghosts at walls | 2 |
| Convection $(\mathbf{u}\cdot\nabla)\mathbf{u}$ | central differences (CDS, default); first-order upwind available | 2 (CDS) |
| Convection $\nabla\cdot(\mathbf{u}\theta)$ | conservative finite volume, first-order upwind (default) | 1 |
| Pressure Laplacian | finite volume, 5-point | 2 (1 on stretched grids, see [mesh.md §4](mesh.md#4-mesh-stretching-refinement)) |
| Gradient, divergence | staggered central differences / face fluxes | 2 |

The convection schemes are chosen in [main.cpp](../../src/main.cpp)
(`MomentumSolver(..., nu)` uses the CDS default; `EnergySolver(..., UPWIND)`),
not in the input files. Details of each stencil:
[solver.md](solver.md) and
[../reference/operators_verification.md](../reference/operators_verification.md).

---

## 6. Stability of the explicit scheme

Both convection and diffusion are explicit, so the time step is limited by

$$
\Delta t \lesssim \frac{h_{\min}^2}{8\,\nu_{\max}}\quad\text{(AB2, 2D diffusion)},
\qquad
\Delta t \lesssim \frac{h_{\min}}{|\mathbf{u}|_{\max}}\quad\text{(convection)},
$$

with $\nu_{\max} = \max(\nu, \alpha)$ and $h_{\min}$ the smallest cell size.
The diffusion limit is the binding one for the cavities at low $Re$ and for
the heated cavity, where $\nu = Pr$ and $\alpha = 1$ are $O(1)$: on a
uniform 40×40 grid $h^2/8 \approx 8\times10^{-5}$, and the validated runs
used $\Delta t = 2\times10^{-5}$. Refining or stretching the grid shrinks
$h_{\min}$ and tightens the limit quadratically. A diverging run stops with
an error and writes no fields ([main.cpp:415-420](../../src/main.cpp#L415-L420)).

---

## 7. Derived quantities used for validation

| Quantity | Definition | Computed in |
|---|---|---|
| Stream function $\psi$ | $u = \partial\psi/\partial y$, $v = -\partial\psi/\partial x$, $\psi = 0$ at the SW corner, integrated exactly on the MAC grid | `stream_function` in [plotters/nsplot.py](../../plotters/nsplot.py) |
| Vorticity $\omega$ | $\partial v/\partial x - \partial u/\partial y$ at cell corners | `vorticity` in [plotters/nsplot.py](../../plotters/nsplot.py) |
| Mean Nusselt number (hot wall) | $\overline{Nu} = \frac{1}{H}\sum_j \frac{\theta_w - \theta_{1,j}}{x_{1} - x_w}\,\Delta y_j$ | [main.cpp:628-640](../../src/main.cpp#L628-L640) |
| $u_{\max}$, $v_{\max}$ (DHC) | max $\lvert u\rvert$ on $x = L/2$, max $\lvert v\rvert$ on $y = H/2$ | [main.cpp:606-626](../../src/main.cpp#L606-L626) |
| Drag / lift coefficients | $C_D = F_x / (\tfrac12\rho U^2 D\,W)$, $C_L = F_y / (\tfrac12\rho U^2 D\,W)$, pressure + viscous shear on the four obstacle faces | [ObstacleForces.cpp](../../src/postprocessing/ObstacleForces.cpp) |
| Strouhal number | $St = f D/U$, $f$ from the lift signal (zero crossings in the solver; FFT and zero crossings over the saturated window in the plotter) | `estimate_shedding_frequency`, [plotters/plot_square_cylinder.py](../../plotters/plot_square_cylinder.py) |
| Recirculation length $L_r/D$ | from the rear face to the first point downstream where $u$ on the centreline becomes positive | `recirculation_length` in [plotters/plot_square_cylinder.py](../../plotters/plot_square_cylinder.py) |

$W$ = `W_depth` is the out-of-plane depth used for face areas and volumes; it
cancels in every non-dimensional result.
