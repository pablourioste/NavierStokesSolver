# Square Cylinder — Validation References and Plan

This page collects the literature values against which the Square Cylinder
scenario must be validated, states the exact configuration each reference
uses, and lists what the solver must change to reproduce it. Every number
below was checked against the cited source; where a value is only shown
graphically in the original, that is said explicitly.

## 1. Status

**Benchmark A (Breuer et al. 2000) is reproduced in the steady regime.** The
validation set-up lives in its own scenario, `square_cylinder_breuer`
(see §5); at $Re = 40$ the solver gives $C_D = 1.717$ and $L_r/D = 2.150$
against 1.708 and 2.151, grid-independent to 0.1 % and 0.5 %
respectively. Getting there required four solver fixes (§4). The unsteady
cases ($Re = 50$–$100$) and Benchmark B are not yet validated.

A reference value is only meaningful for the configuration it was computed
in. Strouhal number and forces depend on the blockage $\beta = D/H$, on the
inflow profile, on the side-wall condition and on the velocity used to define
$Re$. The shipped `square_cylinder` scenario matches **neither** standard
benchmark, which is why the validation uses a separate scenario:

| | `square_cylinder` (shipped) | `square_cylinder_breuer` | Breuer et al. (2000) | Sohankar et al. (1998) |
|---|---|---|---|---|
| Type | confined channel | confined channel | confined channel | quasi-unconfined |
| Blockage $\beta = D/H$ | 1/8 | 1/8 | **1/8** | 5 % (and 2.5 %) |
| Inflow | uniform $U_\infty$ | parabolic | **parabolic**, $U_{\max}$ on the axis | uniform $U_\infty$ |
| Side walls | no-slip | no-slip | no-slip | **slip** (symmetry) |
| $Re$ defined with | $U_\infty$ (= mean) | $U_{\max}$ | **$U_{\max}$** | $U_\infty$ |
| Inlet → cylinder | 6 $D$ | 12 $D$ (face), 12.5 $D$ (centre) | 12.5 $D$ ($L/4$) | 10 $D$ |
| Cylinder → outlet | 18 $D$ | 37 $D$ ($L = 50D$) | 37 $D$ ($L = 50D$) | 10–26 $D$ |
| Outlet condition | zero gradient | zero gradient | — | convective (or Neumann, $X_d = 26D$) |

With a uniform inflow between no-slip walls (shipped scenario), the velocity
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
| **Bug: obstacle-surface velocities entered the pressure equation.** `advance()` updated the nodes on the obstacle faces and the BCs were only re-applied after the projection, so a spurious wall flux entered the divergence of the adjacent cells ($|\nabla\cdot u| \approx U/D$ there). BCs are now re-applied before the PPE and `advance()` skips surface/solid nodes | fixed | [main.cpp](../../src/main.cpp), [MomentumSolver.cpp](../../src/MomentumSolver.cpp) |
| **Bug: inverted obstacle wall distances.** For a node whose neighbour lies in the solid, the image-ghost distance used the *opposite* face of the obstacle ($D + h$ instead of $h$), so the cylinder walls were almost free-slip. Symptom: friction drag that decreased, and turned negative, under grid refinement, and a spurious leading-edge separation at $Re = 40$ | fixed | [DiffusionOperator.cpp](../../src/operators/DiffusionOperator.cpp), [MomentumSolver.cpp](../../src/MomentumSolver.cpp) |
| **Bug: a diverged run was reported as converged.** `std::max(du, NaN)` drops the NaN, so the steady-state test passed and NaN fields were written. Now the case stops with an error and writes nothing | fixed | [main.cpp](../../src/main.cpp) |
| Velocity nodes strictly inside the solid are held at zero (they were never used, but polluted the exported fields) | done | [MomentumSolver.cpp](../../src/MomentumSolver.cpp) |
| OpenMP-parallel pressure solver (Jacobi-PCG), ≈ 3× faster on 8 threads, identical results | done | [solver_sparse.cpp](../../src/solver/sparse/solver_sparse.cpp), `CMakeLists.txt` |
| Real slip wall for Benchmark B: `SYMMETRY` is parsed but `ghostValue()` falls through to the no-slip default | **to do** | [MomentumSolver.cpp](../../src/MomentumSolver.cpp) |

The lid-driven cavity output is byte-identical before and after these
changes, and `operator_verify` reproduces the same convergence data: they
only affect scenarios with an immersed obstacle. **All cylinder results
computed before these fixes are invalid.**

## 5. Results — Benchmark A (`square_cylinder_breuer`)

Run it with `scenario = square_cylinder_breuer` in `input/scenario.txt`, then
`python3 plotters/plot_square_cylinder.py square_cylinder_breuer`. Default
grid: 208×80 cells, $D/16$ across the cylinder, neighbouring cells differing
by < 7 %; $\Delta t = 10^{-3}$ in solver units, i.e. $0.008\,D/U$ per step.
Note the explicit AB2 diffusion limit $\Delta t \lesssim h^2/(8\nu)$ when
refining.

**Steady regime** (D/16, run to $\lVert\Delta u\rVert_\infty/\Delta t < 10^{-4}$):

| $Re$ | $C_D$ | $L_r/D$ | Breuer, $L_r/D$ | Breuer, $C_D$ |
|---|---|---|---|---|
| 20 | 2.366 | 1.017 | 1.043 (fit) | — |
| 30 | 1.944 | 1.586 | 1.597 (fit) | — |
| 40 | 1.716 | 2.143 | 2.151 | 1.708 |

**Grid study at $Re = 40$** (`output/square_cylinder_breuer/grid_study_Re40.csv`):

| Grid | Cells | $C_D$ | of which friction | $L_r/D$ |
|---|---|---|---|---|
| $D/12$ | 156×60 | 1.7175 | 0.320 | 2.132 |
| $D/16$ | 208×80 | 1.7155 | 0.308 | 2.143 |
| $D/20$ | 260×100 | 1.7173 | 0.299 | 2.150 |
| Breuer et al. | | 1.708 | | 2.151 |

$C_D$ is grid-independent to 0.1 % and 0.5–0.6 % above the reference;
$L_r/D$ converges monotonically to the reference (0.05 % on $D/20$). Halving
$\Delta t$ changes $C_D$ by ≈ 1 % (splitting error of the non-incremental
projection) and $L_r$ by < 0.5 %. Mean lift is zero to round-off
($|C_L| < 10^{-6}$).

**Onset of shedding** (D/16): at $Re = 50$ (run to $tU/D = 256$) the lift
perturbation stays at round-off level ($\approx 3\times10^{-6}$) — steady,
$C_D = 1.571$, $L_r/D = 2.681$ (fit: 2.705). At $Re = 70$ it grows and
saturates (run to $tU/D = 648$; saturated from $tU/D \approx 412$, 30
periods), with $St = 0.126$. The onset therefore lies between 50 and 70,
bracketing Breuer's $Re_{cr} \approx 60$.

**Periodic shedding** (D/16 unless stated; statistics over the saturated
part of the record, from `t_star_stats_from` in `force_summary.csv`):

| $Re$ | Window $tU/D$ | $St$ (FFT) | $St$ (zero crossings) | $\overline{C_D}$ | $C_{L,\mathrm{rms}}$ | $\overline{C_L}$ | Breuer $St$ |
|---|---|---|---|---|---|---|---|
| 70 | 412–648 | 0.1260 | 0.1260 | 1.437 | 0.0702 | $-3\times10^{-4}$ | — |
| 100 | 187–400 | **0.1389** | 0.1389 | 1.372 | 0.1377 | $9\times10^{-5}$ | **0.137** |
| 100 ($D/20$) | 186–400 | 0.1391 | 0.1391 | 1.359 | 0.1369 | | |
| 150 | 140–300 | 0.1486 | 0.1487 | 1.355 | 0.2260 | $-2\times10^{-3}$ | — |
| 200 | 138–300 | 0.1474 | 0.1475 | 1.386 | 0.3602 | $2\times10^{-3}$ | — |

At $Re = 100$ the Strouhal number is 1.4 % above the reference and changes by
0.2 % between $D/16$ and $D/20$; the mean lift is zero to $10^{-4}$.
$\overline{C_D}$ at $Re = 100$ differs by 1 % between the two grids. Breuer
et al. give $St(Re)$, $C_D$ and $C_L$ at $Re > 40$ only graphically, so no
numerical comparison is possible for $Re = 70, 150, 200$ without the original
paper. The FFT resolution is $\Delta St = 1/T_{\mathrm{window}}$ (≈ 0.005
at $Re = 100$); the zero-crossing estimate agrees with it to $10^{-4}$.

## 6. Next steps

1. Obtain the original Breuer et al. (2000) paper (tables/figures of $St(Re)$
   and $C_D(Re)$) to compare $Re = 70, 150, 200$ quantitatively.
2. Implement the slip wall and set up Benchmark B (Sohankar et al. 1998,
   $\beta = 5\,\%$) at $Re = 100$ and 150 for $St$, $\overline{C_D}$ and
   $C_L'$, with a grid and time-step study.
3. Consider an incremental pressure-correction scheme (removes the ≈ 1 %
   $\Delta t$ dependence of the steady drag) and a stronger preconditioner
   for the pressure solver (the Jacobi-PCG hits its iteration cap on the
   finer grids).

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
