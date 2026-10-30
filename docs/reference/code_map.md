# Code Map — Which File Does What

Every file of the repository that is part of the solver, its tools or its
post-processing, with its role, its main functions and the scenarios that use
it. For the order in which these pieces run in a time step see
[../theory/solver.md §0](../theory/solver.md#0-overview-the-per-timestep-pipeline);
for a function-by-function walkthrough see
[algorithm_walkthrough.md](algorithm_walkthrough.md).

Scenario columns: **LDC** lid-driven cavity, **DHC** differentially heated
cavity, **CYL** square cylinder (both variants).

---

## 1. C++ solver — `src/`

### Driver

| File | Role | Main contents | LDC | DHC | CYL |
|---|---|---|:-:|:-:|:-:|
| [main.cpp](../../src/main.cpp) | the program: reads the scenario, builds the mesh, runs the case sweep and the time loop, writes results | `read_scenario_name`, `read_obstacle_bounds`, `read_reynolds_cases`, `read_rayleigh_cases`, `main()` (§0 scenario, §1 mesh, pressure matrix, case loop, time loop, VTK/CSV output, Ghia / de Vahl Davis metrics) | ✓ | ✓ | ✓ |
| [reynolds_sweep.cpp](../../src/reynolds_sweep.cpp) | `#include "main.cpp"` — the same program under another executable name | — | ✓ | ✓ | ✓ |
| [mainUtils.hpp](../../src/mainUtils.hpp) | interactive mesh inspector (`run_mesh_inspector`); its call at the end of `main()` is commented out | `run_mesh_inspector` | | | |

### Mesh — `src/mesh/`

| File | Role | Main contents | LDC | DHC | CYL |
|---|---|---|:-:|:-:|:-:|
| [Mesh_config.h](../../src/mesh/Mesh_config.h) | data structures shared by every class | `InternalCell` (centroid, volume, neighbours, `is_solid`), `BCType`, `BoundaryFace`, `GridData` (coordinates, face areas, volumes, `bound_*`, `bound_*_T`, obstacle bounds), class `MeshConfig` | ✓ | ✓ | ✓ |
| [MeshGenerator.cpp](../../src/mesh/MeshGenerator.cpp) | implementation of `MeshConfig`: input parsing, grid generation, boundary conditions, obstacle, VTK export | `generate_stretched_coords` (L79), `generate_mesh` (L152), `generate_u_mesh` (L333), `generate_v_mesh` (L391), `build_connectivity` (L448), `mark_solid_block` (L533), `export_mesh_to_VTK` (L647), `get_data`/`get_stretch`/`get_sim_params` (L708–787), `update_parameter` (L789, list of accepted keys), `inspect_mesh` (L879), `string_to_bc_enum` (L984), `define_boundaries` (L996), `define_temperature_boundaries` (L1086) | ✓ | ✓ | ✓ |
| [ObstacleMask.hpp](../../src/mesh/ObstacleMask.hpp) | header-only helpers: is a P-cell solid, is a U/V node inside the solid or on its surface | `isSolidP`, `uNodeBlocked`, `vNodeBlocked`, `uNodeOnObstacleSurface`, `vNodeOnObstacleSurface` | (inactive) | (inactive) | ✓ |

### Discrete operators — `src/operators/`

| File | Operator | Used by | Verified by |
|---|---|---|---|
| [LaplaceOperator.hpp](../../src/operators/LaplaceOperator.hpp) / [.cpp](../../src/operators/LaplaceOperator.cpp) | $-\nabla^2$ on P (FV, SPD form), Neumann walls, Dirichlet outlet, identity rows in the solid; `pinReference` | `PressurePoissonSystem` | Suites 1, 5, 9, 10 |
| [GradientOperator.hpp](../../src/operators/GradientOperator.hpp) / [.cpp](../../src/operators/GradientOperator.cpp) | $\nabla p$, P → U, V | `PressurePoissonSystem::projectVelocity` | Suites 2, 5, 8 |
| [DivergenceOperator.hpp](../../src/operators/DivergenceOperator.hpp) / [.cpp](../../src/operators/DivergenceOperator.cpp) | $\nabla\cdot\mathbf{u}$ per unit volume, U, V → P; 0 in solid cells | `PressurePoissonSystem::computeRHS` | Suites 3, 5, 8 |
| [DiffusionOperator.hpp](../../src/operators/DiffusionOperator.hpp) / [.cpp](../../src/operators/DiffusionOperator.cpp) | $\nu\nabla^2 u$, $\nu\nabla^2 v$ with image ghosts at walls and at the obstacle | `MomentumSolver` | Suites 4, 6, 7 |

Stencils and the verification suites:
[operators_verification.md](operators_verification.md).

### Solvers

| File | Role | Main contents | LDC | DHC | CYL |
|---|---|---|:-:|:-:|:-:|
| [MomentumSolver.hpp](../../src/MomentumSolver.hpp) / [.cpp](../../src/MomentumSolver.cpp) | predictor $u^*$, $v^*$ (AB2, CDS or upwind convection) and velocity boundary conditions | `ghostValue` (L35), constructor (L58), `applyBoundaryConditions` (L84), `uRHS` (L156), `vRHS` (L257), `advance` (L359); enum `ConvectionScheme` | ✓ | ✓ | ✓ |
| [PressurePoissonSystem.hpp](../../src/PressurePoissonSystem.hpp) / [.cpp](../../src/PressurePoissonSystem.cpp) | pressure equation and projection | `assembleMatrix` (L105), `pinReferencePressure` (L119), `computeRHS` (L147), `projectVelocity` (L190); index helpers `pIdx`/`uIdx`/`vIdx` | ✓ | ✓ | ✓ |
| [EnergySolver.hpp](../../src/EnergySolver.hpp) / [.cpp](../../src/EnergySolver.cpp) | temperature transport on the P grid (AB2, upwind FV convection, image-ghost diffusion) | `extractWallT`, `tRHS` (L48), `advance` (L118) | | ✓ | |
| [Buoyancy.hpp](../../src/Buoyancy.hpp) | Boussinesq coupling: adds `dt·Ra·Pr·θ` to interior $v^*$ | `addBuoyancy` | | ✓ | |

### Linear algebra — `src/solver/`

| File | Role | Built? |
|---|---|---|
| [sparse/solver_class_sparse.hpp](../../src/solver/sparse/solver_class_sparse.hpp) | `Vector`, `Triplet`, `CSRMatrix` (`fromTriplets`), solver classes `Jacobi`, `GS`, `SOR`, `CG`, `PCG` (all with a warm-start overload) | yes |
| [sparse/solver_sparse.cpp](../../src/solver/sparse/solver_sparse.cpp) | implementations; `CG::solve` (L195, used by LDC/DHC), `PCG::solve` (L293, Jacobi-preconditioned, OpenMP, used by CYL); tolerance $10^{-6}\lVert b\rVert$ | yes |
| [sparse/main_solver_sparse.cpp](../../src/solver/sparse/main_solver_sparse.cpp) | stand-alone test program for the sparse solvers | no |
| [general/](../../src/solver/general/) (`solver_class.hpp`, `solver_coding.cpp`, `main_solver.cpp`) | earlier dense-matrix versions of the same solvers | no |

### Post-processing and verification

| File | Role | Used by |
|---|---|---|
| [postprocessing/ObstacleForces.hpp](../../src/postprocessing/ObstacleForces.hpp) / [.cpp](../../src/postprocessing/ObstacleForces.cpp) | `computeForces`: pressure + shear on the four obstacle faces → $F_x$, $F_y$; `estimate_shedding_frequency`: zero crossings of $C_L$ | CYL |
| [verification/verify_operators.cpp](../../src/verification/verify_operators.cpp) | executable `operator_verify`: 10 MMS suites, writes `output/verification/data/convergence.csv` | operator verification |

---

## 2. Build and run scripts

| File | Role |
|---|---|
| [CMakeLists.txt](../../CMakeLists.txt) | builds `mesh_gen`, `reynolds_sweep` and `operator_verify` into `build/`; links OpenMP when found |
| [run.sh](../../run.sh) | `bash run.sh [mesh_gen\|verify\|sweep\|build]`: configure, build and optionally run |
| [tools/run_case.sh](../../tools/run_case.sh) | runs one case in an isolated copy of `input/` with overrides, stores fields, log and inputs in `output/` |
| [tools/run_pending_cases.sh](../../tools/run_pending_cases.sh) | the batch that produced the stored validation runs, then regenerates every figure |
| [tools/status.sh](../../tools/status.sh) | progress of running `run_case.sh` jobs and list of finished/failed ones |
| [requirements.txt](../../requirements.txt) | Python packages for the plotters |

Usage: [../getting_started.md](../getting_started.md).

---

## 3. Python post-processing — `plotters/`

| File | Reads | Writes | Scenario |
|---|---|---|---|
| [nsplot.py](../../plotters/nsplot.py) | VTK fields (`read_vtk`, `load_staggered`) | — (shared style `apply_style`/`save`, `stream_function`, `vorticity`, `vortex_centre`, case discovery) | all |
| [plot_results.py](../../plotters/plot_results.py) | `output/lid_driven_cavity/Re_*/` | per-case `flow_fields`, `centerlines`; `figures/centerline_comparison`, `centerline_trend`, `streamfunction_overview`; `validation_summary.csv` | LDC |
| [plot_dhc.py](../../plotters/plot_dhc.py) | `output/differentially_heated_cavity/Ra_*/` | per-case `flow_fields`; `figures/fields_overview`, `benchmark`; `validation_summary.csv` | DHC |
| [plot_square_cylinder.py](../../plotters/plot_square_cylinder.py) `[scenario]` | `output/<scenario>/Re_*/` | per-case `vorticity`, `force_coefficients`, `lift_spectrum`, `recirculation`; `figures/strouhal`, `recirculation_length`; `force_summary.csv` | CYL |
| [plot_convergence.py](../../plotters/plot_convergence.py) | `output/verification/data/convergence.csv` | `output/verification/figures/convergence_rates` | verification |

Every figure is saved as a vector PDF (for LaTeX, 6.3 in wide, 9 pt Computer
Modern) and a 300 dpi PNG. Output details: [../outputs.md](../outputs.md).

---

## 4. Inputs and documentation

| Path | Content |
|---|---|
| [input/scenario.txt](../../input/scenario.txt) | active scenario |
| `input/<scenario>/` | per-scenario configuration — every key in [../input_reference.md](../input_reference.md) |
| [docs/](../) | this documentation; index in [docs/README.md](../README.md) |
| [docs/figures/](../figures/) | curated figures shown in the README and the scenario pages |
| [docs/legacy/Lid_driven_structured.cpp](../legacy/Lid_driven_structured.cpp) | original single-file lid-driven cavity prototype, used as reference for the ghost-cell treatment; not built |
