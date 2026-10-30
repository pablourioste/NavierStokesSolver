# NavierStokesSolver — Documentation

Start here. Every question below points to the one page that answers it.

## Where do I find …?

| I want to … | Page |
|---|---|
| build the code, run a case, plot it | [getting_started.md](getting_started.md) |
| know what every input file and key does | [input_reference.md](input_reference.md) |
| know what the solver and the plotters write | [outputs.md](outputs.md) |
| see the equations solved in each scenario | [theory/equations.md](theory/equations.md) |
| understand the staggered grid, zones and stretching | [theory/mesh.md](theory/mesh.md) |
| see every boundary condition and how it is implemented | [theory/boundary_conditions.md](theory/boundary_conditions.md) |
| understand one time step, class by class | [theory/solver.md](theory/solver.md) |
| run or understand one scenario | [scenarios/](#scenarios) (table below) |
| find which source file does what | [reference/code_map.md](reference/code_map.md) |
| follow the algorithm function by function | [reference/algorithm_walkthrough.md](reference/algorithm_walkthrough.md) |
| see the discrete operators and their MMS verification | [reference/operators_verification.md](reference/operators_verification.md) |
| follow the heated-cavity implementation line by line | [reference/dhc_implementation.md](reference/dhc_implementation.md) |
| know which bugs were fixed and which results they invalidate | [reference/development_log.md](reference/development_log.md) |

## Scenarios

| # | Scenario | Input folder | Equations | Boundary conditions | Reference | Page |
|---|---|---|---|---|---|---|
| 1 | Lid-driven cavity | `input/lid_driven_cavity/` | continuity + momentum | moving lid $u = 1$; no-slip walls | Ghia et al. (1982) | [lid_driven_cavity.md](scenarios/lid_driven_cavity.md) |
| 2 | Differentially heated cavity | `input/differentially_heated_cavity/` | + Boussinesq buoyancy + energy | no-slip walls; $\theta = 1$ west, $\theta = 0$ east, adiabatic top/bottom | de Vahl Davis (1983) | [differentially_heated_cavity.md](scenarios/differentially_heated_cavity.md) |
| 3 | Square cylinder, benchmark | `input/square_cylinder_breuer/` | continuity + momentum, immersed obstacle | parabolic inlet, zero-gradient outlet with $p = 0$, no-slip walls and cylinder | Breuer et al. (2000) | [square_cylinder.md](scenarios/square_cylinder.md) |
| 3b | Square cylinder, variant | `input/square_cylinder/` | 〃 | uniform inlet, otherwise as 3 | none | [square_cylinder.md §9](scenarios/square_cylinder.md#9-the-square_cylinder-variant) |

Literature data and configurations for the cylinder:
[scenarios/square_cylinder_benchmarks.md](scenarios/square_cylinder_benchmarks.md).

Every scenario page follows the same structure: **1** physical problem
(with a sketch) · **2** governing equations · **3** boundary and initial
conditions · **4** numerical set-up (shipped inputs and the settings of the
validated runs) · **5** code path (which functions run, with links) ·
**6** how to run · **7** outputs and post-processing · **8** validation
results · **9** known limitations · **10** references.

## Suggested reading order

1. [getting_started.md](getting_started.md) — get a case running.
2. [theory/equations.md](theory/equations.md) → [theory/mesh.md](theory/mesh.md)
   → [theory/boundary_conditions.md](theory/boundary_conditions.md)
   → [theory/solver.md](theory/solver.md) — how the method works.
3. The scenario pages, in order 1 → 3: each adds one ingredient (buoyancy
   and a transported scalar; then open boundaries and an immersed body).
4. [reference/](reference/) when you need the detail behind a specific
   function.

## Layout of `docs/`

```
docs/
├── README.md                      this index
├── getting_started.md             build, run, tools, plotting, troubleshooting
├── input_reference.md             every input file and key
├── outputs.md                     output folders, file formats, summary CSVs
├── theory/
│   ├── equations.md               governing equations, scaling, projection method, stability
│   ├── mesh.md                    staggered MAC grid, zones, stretching, VTK
│   ├── boundary_conditions.md     BC types, input format, implementation per variable
│   └── solver.md                  the time step, class by class
├── scenarios/
│   ├── lid_driven_cavity.md
│   ├── differentially_heated_cavity.md
│   ├── square_cylinder.md
│   └── square_cylinder_benchmarks.md   literature values, configurations, solver changes
├── reference/
│   ├── code_map.md                every source file, its role and its scenarios
│   ├── algorithm_walkthrough.md   function-by-function walkthrough (lid-driven cavity)
│   ├── operators_verification.md  operator stencils and the 10 MMS suites
│   ├── dhc_implementation.md      heated-cavity implementation in depth
│   └── development_log.md         bugs found and fixed
├── figures/                       curated figures used in these pages
└── legacy/Lid_driven_structured.cpp   original single-file prototype (not built)
```

## Conventions

- All quantities are non-dimensional; $\rho = 1$ and the reference length
  and velocity are 1 in every input.
- **Shipped inputs** are the files committed in `input/`; **validated runs**
  are the stored results in `output/`, whose exact inputs are kept in each
  `output/<scenario>/<case>/input/` when they were run with
  `tools/run_case.sh`. The two can differ — each scenario page lists both.
- Code links point to a file and line (`file.cpp#L120`); line numbers refer
  to the version of the code these pages were written against and may drift
  slightly as the code changes.

| Symbol | Meaning |
|---|---|
| $u, v$ | velocity components in $x, y$ |
| $p$ | pressure (kinematic, $\rho = 1$) |
| $\theta$ | non-dimensional temperature $(T - T_c)/(T_h - T_c)$ |
| $Re, Ra, Pr$ | Reynolds, Rayleigh, Prandtl numbers |
| $L, H$ | domain length and height |
| $D$ | cylinder side; $\beta = D/H$ blockage |
| $N, M$ | number of pressure cells in $x$ and $y$ |
| $h$ | local cell size; $h_{\min}$ the smallest |
| $C_D, C_L$ | drag and lift coefficients |
| $St$ | Strouhal number $fD/U$ |
| $L_r$ | recirculation length behind the cylinder |
| $t^*$ | convective time $tU/D$ |
| $\overline{Nu}$ | mean Nusselt number on the hot wall |
