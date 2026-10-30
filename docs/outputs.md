# Outputs — What the Solver and the Plotters Write

Results are organised by scenario, then by case. The **solver** writes
`mesh/`, `fields/` and `data/`; the **plotters** write `figures/` and the
summary CSVs.

```
output/
├── lid_driven_cavity/
│   ├── mesh/                       mesh_P.vtk, mesh_U.vtk, mesh_V.vtk
│   ├── Re_<N>/
│   │   ├── fields/                 pressure.vtk, u_velocity.vtk, v_velocity.vtk
│   │   ├── data/                   centerline_u.csv, centerline_v.csv, run.log*
│   │   ├── input/                  exact inputs of the run*
│   │   └── figures/                flow_fields, centerlines                       (plotter)
│   ├── figures/                    centerline_comparison, centerline_trend,
│   │                               streamfunction_overview                        (plotter)
│   └── validation_summary.csv      errors vs. Ghia et al. (1982)                  (plotter)
├── differentially_heated_cavity/
│   ├── mesh/
│   ├── Ra_<N>/{fields, data, input*, figures}
│   │                               + temperature.vtk, dhc_metrics.csv
│   ├── figures/                    fields_overview, benchmark                     (plotter)
│   └── validation_summary.csv      errors vs. de Vahl Davis (1983)                (plotter)
├── square_cylinder_breuer/         (and square_cylinder/)
│   ├── mesh/                       + obstacle_mask.vtk
│   ├── Re_<N>/{fields, data, input*, figures}
│   │                               + forces.csv, strouhal.csv;
│   │                               figures: vorticity, force_coefficients,
│   │                               lift_spectrum, recirculation                   (plotter)
│   ├── figures/                    strouhal, recirculation_length                 (plotter)
│   ├── force_summary.csv           C_D, C_L, St, L_r and references per case      (plotter)
│   ├── grid_study_Re40.csv         grid study at Re = 40 (Breuer only)
│   └── grid_study/D20/             runs on the D/20 grid (Breuer only)
├── verification/
│   ├── data/convergence.csv        (operator_verify)
│   └── figures/convergence_rates   (plotter)
└── logs/                           console logs of tools/run_pending_cases.sh (git-ignored)
```

\* only for cases run with `tools/run_case.sh`.

`figures/` folders and `output/logs/` are git-ignored (regenerated on
demand). Everything else in `output/` is committed, so the validation data
is available without re-running. The figures used in the documentation are
copied to [docs/figures/](figures/).

---

## 1. VTK files (ParaView)

ASCII legacy VTK, `DATASET RECTILINEAR_GRID`, one `CELL_DATA` scalar per
file. Each variable is written **on its own staggered grid**:

| File | Grid | Cells | Scalar |
|---|---|---|---|
| `pressure.vtk`, `temperature.vtk` | P grid | $N \times M$ | `pressure`, `temperature` |
| `u_velocity.vtk` | U grid (cells centred on vertical faces) | $(N{+}1) \times M$ | `u_velocity` |
| `v_velocity.vtk` | V grid (cells centred on horizontal faces) | $N \times (M{+}1)$ | `v_velocity` |
| `mesh/mesh_P/U/V.vtk` | the three grids | | `cell_ID` |
| `mesh/obstacle_mask.vtk` | P grid | $N \times M$ | `obstacle_mask` (1 = solid) |

To combine $u$ and $v$ in ParaView, resample them onto a common grid (e.g.
*Resample With Dataset* onto `pressure.vtk`). The Python plotters instead
read the three files with `nsplot.load_staggered` and keep each value at its
true location.

The fields are those of the **last time step**: the steady solution for
steady cases, an instantaneous snapshot for vortex-shedding cases.

## 2. CSV files written by the solver

| File | Columns | Notes |
|---|---|---|
| `centerline_u.csv` (LDC) | `y, y_norm, u, u_norm` | $u$ on the U-grid column nearest $x = L/2$; `_norm` = divided by $H$ / by $U_{lid}$ (= 1) |
| `centerline_v.csv` (LDC) | `x, x_norm, v, v_norm` | $v$ on the V-grid row nearest $y = H/2$ |
| `dhc_metrics.csv` (DHC) | `Ra, Pr, u_max, v_max, Nu_avg` | definitions in [scenarios/differentially_heated_cavity.md §7](scenarios/differentially_heated_cavity.md#7-outputs-and-post-processing) |
| `forces.csv` (cylinder) | `step, time, Cd, Cl` | every step; `time` = step·`dt` in solver units (convective time $t^* = tU/D = 8\,t$ for $D = 0.125$) |
| `strouhal.csv` (cylinder) | `Re, St` | solver estimate from zero crossings of $C_L$ over the second half of the record; the plotter's FFT/saturated-window values in `force_summary.csv` are the ones quoted |
| `verification/data/convergence.csv` | `suite, N, error_linf` | written by `operator_verify` |

## 3. Summary CSVs written by the plotters

| File | Columns |
|---|---|
| `lid_driven_cavity/validation_summary.csv` | `Re, grid, psi_min, x_c, y_c` (primary vortex), `psi_min_ghia, x_c_ghia, y_c_ghia, psi_min_rel_err_pct`, `u_rms_err, u_max_err, v_rms_err, v_max_err` (centreline errors vs. Ghia; empty if $Re$ is not tabulated) |
| `differentially_heated_cavity/validation_summary.csv` | `Ra, grid, u_max, v_max, Nu_avg`, the de Vahl Davis references `*_ref` and errors `*_err_pct` |
| `square_cylinder*/force_summary.csv` | `Re, t_star_end, regime` (steady/developing/periodic), `t_star_stats_from` (start of the statistics window), `Cd_mean, Cl_mean, Cl_rms, Lr_D`, `St_fft, St_zc, St_solver`, and the reference values `*_ref` where available |
| `square_cylinder_breuer/grid_study_Re40.csv` | `grid, cells_across_D, Nx, Ny, dt, Cd, Cd_pressure, Cd_friction, Lr_D, Cd_ref, Lr_D_ref` |

## 4. Figures

Every figure is saved twice: `<name>.pdf` (vector, drawn at a 6.3 in text
width with 9 pt Computer Modern, to include in LaTeX at `width=\textwidth`
without rescaling) and `<name>.png` (300 dpi). What each plotter draws is
listed in the scenario pages (§7 of each) and in
[reference/code_map.md §3](reference/code_map.md#3-python-post-processing--plotters).
