# Getting Started — Build, Run, Plot

From a fresh clone to a validated result and its figures. Every command is
run from the **repository root**: the executables resolve `input/` and
`output/` as relative paths.

---

## Contents

1. [Requirements](#1-requirements)
2. [Build](#2-build)
3. [Running the active scenario](#3-running-the-active-scenario)
4. [Running individual cases (`tools/run_case.sh`)](#4-running-individual-cases-toolsrun_casesh)
5. [Batches and progress](#5-batches-and-progress)
6. [Plotting](#6-plotting)
7. [Reproducing the validated results](#7-reproducing-the-validated-results)
8. [Setting up a new case](#8-setting-up-a-new-case)
9. [Troubleshooting](#9-troubleshooting)

---

## 1. Requirements

- CMake ≥ 3.15 and a C++17 compiler (g++ or clang++). OpenMP is optional
  (used by the cylinder pressure solver when found).
- Python 3 with the packages in [requirements.txt](../requirements.txt)
  (`numpy`, `pandas`, `matplotlib`, `scipy`) for the plotters. A LaTeX
  installation is optional: without it the figures fall back to matplotlib's
  Computer Modern mathtext.

```bash
python3 -m venv .venv
.venv/bin/pip install -r requirements.txt
```

`.vscode/settings.json` points VS Code at `.venv/bin/python3`. If your
default `python3` is another environment, call
`.venv/bin/python3 plotters/<script>.py` explicitly.

## 2. Build

```bash
bash run.sh build        # cmake configure (Release) + build, runs nothing
```

or by hand:

```bash
cmake -S . -B build_cmake -DCMAKE_BUILD_TYPE=Release
cmake --build build_cmake
```

| Executable (in `build/`) | Purpose |
|---|---|
| `mesh_gen` | the solver: runs the scenario named in `input/scenario.txt` |
| `reynolds_sweep` | the same program under another name |
| `operator_verify` | MMS verification of the discrete operators → `output/verification/data/convergence.csv` |

Rebuild after any change in `src/`.

## 3. Running the active scenario

1. Choose the scenario in [input/scenario.txt](../input/scenario.txt):
   `lid_driven_cavity`, `differentially_heated_cavity`, `square_cylinder` or
   `square_cylinder_breuer`.
2. Edit the files in `input/<scenario>/` if needed (all keys:
   [input_reference.md](input_reference.md)); the Re or Ra values to sweep
   are in `reynolds_cases.txt` / `rayleigh_cases.txt`.
3. Run:

```bash
bash run.sh              # build + run mesh_gen
bash run.sh mesh_gen     # same, explicit
bash run.sh verify       # build + run operator_verify
```

Every case of the sweep starts from rest and writes to
`output/<scenario>/<Re_N | Ra_N>/` ([outputs.md](outputs.md)). The console
prints $\lVert\Delta u\rVert_\infty/\Delta t$ every 100 steps; the case ends
at `max_steps` or when that residual drops below `ss_tol` (never, if
`ss_tol = 0`). A case that produces NaN stops with `ERROR ... diverged` and
writes nothing.

## 4. Running individual cases (`tools/run_case.sh`)

Recommended for real runs: it copies `input/` and the executable to a
temporary folder, applies your overrides there (the committed inputs are
never edited), runs one case, and stores the result together with the exact
inputs used and the log.

```bash
tools/run_case.sh SCENARIO VALUE [options]
```

| Option | Effect |
|---|---|
| `--set KEY=VALUE` | override a key of `sim_params.txt` (`dt`, `max_steps`, `ss_tol`, …); repeatable |
| `--stretch KEY=VALUE` | override a key of `stretch.txt` (cell counts, `beta_*`); repeatable |
| `--scale F` | multiply every cell count in `stretch.txt` by `F` (grid studies) |
| `--threads N` | OpenMP threads for the pressure solver (default 4; the cavity CG is serial, use 1) |
| `--store SUBDIR` | store under `output/SCENARIO/SUBDIR/` instead of `output/SCENARIO/` |

Result: `output/SCENARIO[/SUBDIR]/<Re_N|Ra_N>/{fields,data,input}/`, log in
`data/run.log` (PCG iteration-cap warnings filtered out). If the case
diverges, nothing is stored, `FAILED` is printed and the log stays in the
temporary folder. Several cases can run at once; keep the total of
`--threads` below the number of cores. If the executable is missing, the
script asks you to run `bash run.sh build` first.

If `tools/*.sh` are not executable in your clone, run them with
`bash tools/run_case.sh ...`.

## 5. Batches and progress

```bash
nohup bash tools/run_pending_cases.sh > output/logs/pending.log 2>&1 &
bash tools/status.sh                 # snapshot of running and finished cases
watch -n 30 bash tools/status.sh     # refresh every 30 s
```

[tools/run_pending_cases.sh](../tools/run_pending_cases.sh) is the batch that
produced the stored validation runs (cylinder $Re = 70$–$200$ and the $D/20$
check, cavity $Re = 100$–$7500$, heated cavity $Ra = 10^6$); it ends by
regenerating every figure. It is also the record of the exact commands used.
`output/logs/` is git-ignored.

## 6. Plotting

```bash
python3 plotters/plot_results.py                                  # lid-driven cavity
python3 plotters/plot_dhc.py                                      # heated cavity
python3 plotters/plot_square_cylinder.py square_cylinder_breuer   # cylinder, benchmark
python3 plotters/plot_square_cylinder.py                          # cylinder, variant
python3 plotters/plot_convergence.py                              # operator verification
```

The plotters find every case present in `output/<scenario>/`, write per-case
figures to `<case>/figures/`, cross-case figures to
`output/<scenario>/figures/` and a summary CSV with the benchmark errors.
`figures/` folders are git-ignored; curated copies used in the docs are in
[docs/figures/](figures/). Details: [outputs.md](outputs.md).

## 7. Reproducing the validated results

| Scenario | Command(s) | Page |
|---|---|---|
| Lid-driven cavity, $Re = 100$ | `tools/run_case.sh lid_driven_cavity 100 --set dt=3e-4 --set max_steps=100000 --set ss_tol=1e-5 --threads 1` | [scenarios/lid_driven_cavity.md](scenarios/lid_driven_cavity.md) |
| Lid-driven cavity, $Re = 400$–$5000$ | `tools/run_case.sh lid_driven_cavity <Re> --set dt=1e-3 --set max_steps=200000 --set ss_tol=1e-5 --threads 1` | 〃 |
| Heated cavity, $Ra = 10^3$–$10^5$ | `tools/run_case.sh differentially_heated_cavity <Ra> --stretch N_cells_west=40 --stretch M_cells_west=40 --stretch beta_x_west=0.0 --stretch beta_y_south=0.0 --set dt=2e-5 --threads 1` | [scenarios/differentially_heated_cavity.md](scenarios/differentially_heated_cavity.md) |
| Heated cavity, $Ra = 10^6$ | as above with `N/M_cells_west=100`, `--set dt=2.5e-6 --set max_steps=160000 --set ss_tol=0.0` | 〃 |
| Cylinder, steady ($Re \leq 40$) | `tools/run_case.sh square_cylinder_breuer <Re> --threads 4` | [scenarios/square_cylinder.md](scenarios/square_cylinder.md) |
| Cylinder, periodic ($Re = 100$) | `tools/run_case.sh square_cylinder_breuer 100 --set max_steps=50000 --set ss_tol=0.0 --threads 4` | 〃 |
| Operators (MMS) | `bash run.sh verify && python3 plotters/plot_convergence.py` | [reference/operators_verification.md](reference/operators_verification.md) |

For orientation, [tools/run_pending_cases.sh](../tools/run_pending_cases.sh)
notes about 2.5 h for its cylinder batch and 3 h for its cavity batch
(cases run in parallel) on the 16-core development machine.

## 8. Setting up a new case

1. Copy an existing folder: `cp -r input/square_cylinder_breuer input/my_case`.
   The scenario type follows from the name: `differentially_heated_cavity`
   enables the energy equation, any name starting with `square_cylinder`
   enables the obstacle, anything else runs as a lid-driven-cavity-type
   problem (no obstacle, no energy, Re sweep, Ghia comparison printed).
2. Edit geometry (`data.txt`), grid (`stretch.txt`), parameters
   (`sim_params.txt`), boundary conditions (`boundaries.txt`) — see
   [input_reference.md](input_reference.md) and
   [theory/boundary_conditions.md](theory/boundary_conditions.md).
3. Check the explicit time-step limit for the new grid
   ([theory/equations.md §6](theory/equations.md#6-stability-of-the-explicit-scheme)).
4. Run a short test (`--set max_steps=500`), open
   `output/my_case/mesh/*.vtk` in ParaView to check the grid (and
   `obstacle_mask.vtk` for an obstacle), then run the full case.

## 9. Troubleshooting

| Symptom | Likely cause | What to do |
|---|---|---|
| `ERROR: ... diverged (non-finite velocity)` | `dt` above the explicit stability limit | reduce `dt` below $h_{\min}^2/(8\nu_{\max})$; see the scenario page |
| `CG WARNING: did not converge ... = -nan` | the run is already diverging | same as above |
| repeated `CG WARNING` / `PCG WARNING` with a small residual (≈ 1e-5) | iteration cap reached on a large or stretched grid; the step continues with the current iterate | usually harmless if small; more threads speed up PCG |
| `Error: could not open input/...` | executable not run from the repository root, or wrong scenario name | `cd` to the root; check `input/scenario.txt` |
| case never stops | `ss_tol = 0.0` disables the steady-state exit | set `ss_tol` (e.g. `1e-5`) or reduce `max_steps` |
| `WARNING mark_solid_block: requested coordinate ... snapped` | obstacle bounds not on a cell face | make `obstacle.txt` match the zone corners in `data.txt` |
| plotters fail to import `numpy`/`pandas` | wrong Python interpreter | use `.venv/bin/python3` |
