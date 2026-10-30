# Square Cylinder — Benchmarks, Reference Data and Solver Changes

Companion to [square_cylinder.md](square_cylinder.md). This page collects
the literature values the Square Cylinder scenario is compared with, the
exact configuration each reference uses, and the solver changes (including
four bug fixes) that were needed to reproduce them. Every number below was
checked against the cited source; where a value is only shown graphically in
the original, that is said explicitly.

**The results of the present solver** (tables and figures) are in
[square_cylinder.md §8](square_cylinder.md#8-validation-results-square_cylinder_breuer).

## Contents

1. [Why a separate benchmark scenario](#1-why-a-separate-benchmark-scenario)
2. [Benchmark A — confined cylinder, Breuer et al. (2000)](#2-benchmark-a--confined-cylinder-breuer-et-al-2000)
3. [Benchmark B — quasi-unconfined cylinder, Sohankar et al. (1998)](#3-benchmark-b--quasi-unconfined-cylinder-sohankar-et-al-1998)
4. [Solver changes made for the validation](#4-solver-changes-made-for-the-validation)
5. [Sources](#sources)

## 1. Why a separate benchmark scenario

A reference value is only meaningful for the configuration it was computed
in. Strouhal number and forces depend on the blockage $\beta = D/H$, on the
inflow profile, on the side-wall condition and on the velocity used to define
$Re$. The `square_cylinder` variant matches **neither** standard
benchmark, which is why the validation uses a separate scenario:

| | `square_cylinder` (variant) | `square_cylinder_breuer` | Breuer et al. (2000) | Sohankar et al. (1998) |
|---|---|---|---|---|
| Type | confined channel | confined channel | confined channel | quasi-unconfined |
| Blockage $\beta = D/H$ | 1/8 | 1/8 | **1/8** | 5 % (and 2.5 %) |
| Inflow | uniform $U_\infty$ | parabolic | **parabolic**, $U_{\max}$ on the axis | uniform $U_\infty$ |
| Side walls | no-slip | no-slip | no-slip | **slip** (symmetry) |
| $Re$ defined with | $U_\infty$ (= mean) | $U_{\max}$ | **$U_{\max}$** | $U_\infty$ |
| Inlet → cylinder | 6 $D$ | 12 $D$ (face), 12.5 $D$ (centre) | 12.5 $D$ ($L/4$) | 10 $D$ |
| Cylinder → outlet | 18 $D$ | 37 $D$ ($L = 50D$) | 37 $D$ ($L = 50D$) | 10–26 $D$ |
| Outlet condition | zero gradient | zero gradient | — | convective (or Neumann, $X_d = 26D$) |

With a uniform inflow between no-slip walls (`square_cylinder` variant), the velocity
profile develops a boundary layer along the 6 $D$ before the cylinder, so
the cylinder sees neither a uniform nor a parabolic stream and the effective
Reynolds number is undefined.

## 2. Benchmark A — confined cylinder, Breuer et al. (2000)

M. Breuer, J. Bernsdorf, T. Zeiser, F. Durst, *Accurate computations of the
laminar flow past a square cylinder based on two different methods:
lattice-Boltzmann and finite-volume*, Int. J. Heat Fluid Flow 21 (2000)
186–196, doi:10.1016/S0142-727X(99)00081-8.

**Configuration** (as restated by Lamura et al. 2001, who reproduce it):
$\beta = D/H = 1/8$; channel length $L = 50D$; cylinder centred in the
channel at $L/4$ from the inlet; parabolic inflow
$u(y) = 4U_{\max}\,y(H-y)/H^2$; no-slip channel and cylinder walls;
$Re = U_{\max}D/\nu$, $St = fD/U_{\max}$, coefficients normalised with
$U_{\max}$.

| Quantity | Reference value | Source of the number |
|---|---|---|
| $L_r/D$ (steady, $Re = 40$) | 2.151 | Breuer et al. via Korukçu (2020), Table 2 |
| $C_D$ (steady, $Re = 40$) | 1.708 | Breuer et al. via Korukçu (2020), Table 2 |
| $L_r/D$ vs $Re$, steady range | $-0.065 + 0.0554\,Re$ | Breuer et al. via Lamura et al. (2001), Eq. (4) |
| Onset of shedding | $Re_{cr} \approx 60$ | Breuer et al. via Lamura et al. (2001), §4 |
| $St$ at $Re = 100$ | ≈ 0.137 | Chen, Xu & Li (2015), §IV.4: "about 0.137, which is identical to the reference data" (Breuer) |
| $C_D$, $St$ at other $Re$ | figures only | Breuer et al. Figs.; values not available in open sources — **read from the original paper** |

The two steady-state values at $Re = 40$ are consistent with each other
($-0.065 + 0.0554 \times 40 = 2.151$).

## 3. Benchmark B — quasi-unconfined cylinder, Sohankar et al. (1998)

A. Sohankar, C. Norberg, L. Davidson, *Low-Reynolds-number flow around a
square cylinder at incidence: study of blockage, onset of vortex shedding and
outlet boundary condition*, Int. J. Numer. Meth. Fluids 26 (1998) 39–56.

**Configuration:** uniform inflow $U_\infty$; slip (frictionless) side
boundaries; $H = 20D$ ($\beta = 5\,\%$) or $40D$ ($\beta = 2.5\,\%$); inlet
$X_u = 10D$ upstream; convective outlet $X_d = 10D$ downstream (or Neumann
outlet with $X_d = 26D$); $Re = U_\infty D/\nu$; QUICK convection,
Crank–Nicolson, $\Delta t = 0.025\,D/U_\infty$.

Table III of the paper (zero incidence; $C_L'$ is the rms lift):

| $Re$ | $\beta$ | $St$ | $\overline{C_D}$ | $C_L'$ |
|---|---|---|---|---|
| 100 | 5 % | 0.146 | 1.460 | 0.139 |
| 100 | 2.5 % | 0.144 | 1.444 | 0.141 |
| 150 | 5 % | 0.162 | 1.417 | 0.160 |
| 150 | 2.5 % | 0.160 | 1.400 | 0.166 |
| 200 | 5 % | 0.165 | 1.427 | 0.232 |

Grid-refined values at $Re = 100$ (Table IV, finest case): $St = 0.146$,
$\overline{C_D} = 1.478$, $C_L' = 0.153$. Onset of shedding:
$Re_{cr} = 51.2 \pm 1.0$ at $\beta = 5\,\%$ (§4.4). Also from the paper:
increasing the blockage increases $St$ and $C_D$, "similar to an increase in
the Reynolds number".

Other values at $Re = 100$ for cross-checking (as tabulated in Haque 2023,
Table III, unconfined uniform inflow): Sharma & Eswaran (2004)
$\overline{C_D} = 1.494$, $C_L' = 0.1922$, $St = 0.1488$; Zhao (2013)
$1.452$ / $0.1908$ / $0.1447$; Berrone (2011) $1.48$ / — / $0.145$;
experiments: Okajima (1982) $St = 0.14$–$0.144$ at $Re = 120$, Norberg
(1996) $St = 0.142$ ($Re = 90$) and $0.147$ ($Re = 110$).

## 4. Solver changes made for the validation

| Change | Status | Where |
|---|---|---|
| Parabolic inlet `INLET_PARABOLIC` (value = $U_{\max}$), per-face $u(y) = 4U_{\max}(y-y_0)(y_1-y)/H^2$ | done | `MeshConfig::define_boundaries`, [MeshGenerator.cpp](../../src/mesh/MeshGenerator.cpp) |
| Scenario names `square_cylinder*` run the cylinder solver, so benchmark set-ups get their own `input/` and `output/` folders | done | [main.cpp](../../src/main.cpp) |
| **Bug: viscous force components swapped.** Shear $\mu\,\partial v/\partial x$ on the front/rear faces was added to the drag and $\mu\,\partial u/\partial y$ on the top/bottom faces to the lift. Symptom: non-zero mean lift on a symmetric body ($C_L = 9.3$ at step 1) and drag missing its friction part | fixed | [ObstacleForces.cpp](../../src/postprocessing/ObstacleForces.cpp) |
| **Bug: obstacle-surface velocities entered the pressure equation.** `advance()` updated the nodes on the obstacle faces and the BCs were only re-applied after the projection, so a spurious wall flux entered the divergence of the adjacent cells ($\lvert\nabla\cdot u\rvert \approx U/D$ there). BCs are now re-applied before the PPE and `advance()` skips surface/solid nodes | fixed | [main.cpp](../../src/main.cpp), [MomentumSolver.cpp](../../src/MomentumSolver.cpp) |
| **Bug: inverted obstacle wall distances.** For a node whose neighbour lies in the solid, the image-ghost distance used the *opposite* face of the obstacle ($D + h$ instead of $h$), so the cylinder walls were almost free-slip. Symptom: friction drag that decreased, and turned negative, under grid refinement, and a spurious leading-edge separation at $Re = 40$ | fixed | [DiffusionOperator.cpp](../../src/operators/DiffusionOperator.cpp), [MomentumSolver.cpp](../../src/MomentumSolver.cpp) |
| **Bug: a diverged run was reported as converged.** `std::max(du, NaN)` drops the NaN, so the steady-state test passed and NaN fields were written. Now the case stops with an error and writes nothing | fixed | [main.cpp](../../src/main.cpp) |
| Velocity nodes strictly inside the solid are held at zero (they were never used, but polluted the exported fields) | done | [MomentumSolver.cpp](../../src/MomentumSolver.cpp) |
| OpenMP-parallel pressure solver (Jacobi-PCG), ≈ 3× faster on 8 threads, identical results | done | [solver_sparse.cpp](../../src/solver/sparse/solver_sparse.cpp), `CMakeLists.txt` |
| Real slip wall for Benchmark B: `SYMMETRY` is parsed but `ghostValue()` falls through to the no-slip default | **to do** | [MomentumSolver.cpp](../../src/MomentumSolver.cpp) |

The lid-driven cavity output is byte-identical before and after these
changes, and `operator_verify` reproduces the same convergence data: they
only affect scenarios with an immersed obstacle. **All cylinder results
computed before these fixes are invalid.**

## Sources

- Breuer et al. (2000), the primary confined benchmark:
  [ScienceDirect](https://www.sciencedirect.com/science/article/abs/pii/S0142727X99000818),
  [FAU CRIS record](https://cris.fau.de/publications/119688624/).
- Lamura, Gompper, Ihle, Kroll (2001), *Multi-particle-collision dynamics:
  flow around a circular and a square cylinder*, reproduces the Breuer
  set-up and quotes its $L_r$ fit and $Re_{cr}$:
  [arXiv:cond-mat/0110148](https://arxiv.org/html/cond-mat/0110148).
- Korukçu (2020), Uludağ Univ. J. Eng. 25(1), doi:10.17482/uumfd.630535,
  Table 2 (Breuer values at $Re = 40$, $\beta = 1/8$):
  [PDF](https://acikerisim.uludag.edu.tr/server/api/core/bitstreams/08a47aee-d09a-482a-908a-fee3bbf6fc8e/content).
- Chen, Xu, Li (2015), *Cartesian grid method for gas kinetic scheme*, §IV.4
  ($St$ at $Re = 100$ in the Breuer set-up):
  [arXiv:1502.04494](https://arxiv.org/pdf/1502.04494).
- Sohankar, Norberg, Davidson (1998), Tables III and IV:
  [PDF](https://www.cfd-sweden.se/lada/postscript_files/Sohankar_num-fluids.pdf),
  [Wiley](https://onlinelibrary.wiley.com/doi/abs/10.1002/%28SICI%291097-0363%2819980115%2926%3A1%3C39%3A%3AAID-FLD623%3E3.0.CO%3B2-P).
- Haque (2023), AIP Advances 13, 085219, Table III (compilation at
  $Re = 100$):
  [AIP](https://pubs.aip.org/aip/adv/article/13/8/085219/2906520/Flow-simulation-around-square-cylinder-with).
