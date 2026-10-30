# Input File Reference

Every input file and every key the solver reads. All files are plain text:
`key = value` lines (spaces around `=` required), or whitespace-separated
columns for the boundary and case files. Lines starting with `#` and empty
lines are ignored. An unknown key in `data.txt`, `stretch.txt` or
`sim_params.txt` prints `Warning: unknown key in configuration file: <key>`
and is otherwise ignored, so a misspelt key leaves the default in place —
check the console output after changing an input.

```
input/
├── scenario.txt                  which scenario to run
└── <scenario>/
    ├── data.txt                  domain size and zone corners
    ├── stretch.txt               cells per zone and stretching
    ├── sim_params.txt            physical and numerical parameters
    ├── boundaries.txt            velocity boundary conditions
    ├── boundaries_T.txt          temperature boundary conditions   (DHC only)
    ├── obstacle.txt              solid obstacle bounds             (cylinder only)
    ├── reynolds_cases.txt        Re values to sweep                (LDC, cylinder)
    ├── rayleigh_cases.txt        Ra values to sweep                (DHC)
    └── ghia_reference.txt        reference data, not read          (LDC only)
```

Per-scenario values and their meaning are discussed in each scenario page
([scenarios/](scenarios/)).

---

## `scenario.txt`

```
scenario = square_cylinder_breuer
```

| Value | Runs |
|---|---|
| `lid_driven_cavity` | momentum + pressure, Re sweep, Ghia centrelines |
| `differentially_heated_cavity` | + energy and buoyancy (requires `solve_energy = 1`), Ra sweep |
| `square_cylinder`, `square_cylinder_breuer`, any `square_cylinder*` | + immersed obstacle, forces, Strouhal number |
| anything else (folder must exist) | treated as `lid_driven_cavity` |

Missing file → `lid_driven_cavity`. `tools/run_case.sh` writes its own copy
of this file, so it does not need editing for those runs.

## `data.txt` — domain and zones

| Key | Meaning |
|---|---|
| `L_domain`, `H_domain` | domain length ($x$) and height ($y$) |
| `W_depth` | out-of-plane depth used for face areas and volumes (cancels in every non-dimensional result; keep 0.1) |
| `node_NW_x`, `node_NW_y` | top-left corner of the centre zone: $x_{NW}$ is the west/centre limit, $y_{NW}$ the centre/north limit |
| `node_NE_x`, `node_NE_y` | top-right corner: $x_{NE}$ is the centre/east limit |
| `node_SW_x`, `node_SW_y` | bottom-left corner: $y_{SW}$ is the south/centre limit |
| `node_SE_x`, `node_SE_y` | read but not used |

The domain is split into 3×3 zones by these corners
([theory/mesh.md §5](theory/mesh.md#5-zones-the-domain-is-split-into-a-33-grid-of-regions)):

```
 y = H ┌──────────┬──────────┬────────────┐
       │          │  north   │            │
 y_NW  ├──────────NW────────NE────────────┤
       │   west   │  centre  │    east    │     cylinder: centre zone = obstacle
 y_SW  ├──────────SW─────────SE───────────┤     cavities: every node at (L, H)
       │          │  south   │            │               → one single zone
 y = 0 └──────────┴──────────┴────────────┘
      x = 0     x_NW       x_NE          x = L
```

- **Cavities:** put all four nodes at $(L, H)$ and give the centre/east
  zones 0 cells; the whole domain is then the west/south zone. A node inside
  the domain with 0 cells in its zone creates zero-width cells (NaN).
- **Cylinder:** put the nodes on the obstacle corners so that the centre
  zone is exactly the obstacle, and repeat the same bounds in `obstacle.txt`.

## `stretch.txt` — cells and stretching

| Key | Zone |
|---|---|
| `N_cells_west`, `N_cells_center`, `N_cells_east` | cells in the west, centre, east zones ($x$) |
| `M_cells_west`, `M_cells_center`, `M_cells_east` | cells in the **south**, centre, **north** zones ($y$) — note the suffixes |
| `beta_x_west`, `beta_x_center`, `beta_x_east` | stretching of each $x$ zone |
| `beta_y_south`, `beta_y_center`, `beta_y_north` | stretching of each $y$ zone |

Stretching ([theory/mesh.md §4](theory/mesh.md#4-mesh-stretching-refinement)):
$\beta = 0$ uniform; $\beta > 0$ refines towards the **end** of the zone
(larger $x$ or $y$); $\beta < 0$ towards the **start**; larger $|\beta|$
means stronger refinement. Each zone is refined towards one end only.

## `sim_params.txt` — physical and numerical parameters

| Key | Type | Meaning | Used by |
|---|---|---|---|
| `U_lid` | real | lid velocity; also the velocity scale in $\nu = U_{lid}L/Re$ | LDC |
| `U_inf` | real | inflow velocity scale (uniform $U_\infty$ or parabolic $U_{\max}$); used in $\nu = U_\infty D/Re$, $C_D$, $C_L$, $St$ | cylinder |
| `rho` | real | density (keep 1) | all |
| `dt` | real | time step — limited by the explicit scheme, see [theory/equations.md §6](theory/equations.md#6-stability-of-the-explicit-scheme) | all |
| `max_steps` | int | maximum number of time steps per case | all |
| `ss_tol` | real | stop when $\lVert u^{n+1}-u^n\rVert_\infty/\Delta t <$ `ss_tol` (checked every 100 steps); `0.0` disables the check | all |
| `solve_energy` | 0/1 | energy equation on/off; only takes effect in the `differentially_heated_cavity` scenario | DHC |
| `Pr` | real | Prandtl number; momentum diffusivity is set to `Pr` | DHC |
| `Re` | real | read but **not used**: $Re$ comes from `reynolds_cases.txt` | — |

The viscosity is never read: it is derived per case from the swept number
([theory/equations.md §2](theory/equations.md#2-non-dimensional-form-used-by-each-scenario)).

## `boundaries.txt` — velocity boundary conditions

```
# BOCO_NAME  FACE   IDX_START IDX_END  TYPE               VALUE
lid          NORTH  0         9999     WALL_FIXED_VALUE   1.0
wall_south   SOUTH  0         9999     WALL               0.0
```

`FACE` is `WEST`, `EAST`, `SOUTH` or `NORTH`; `IDX_START IDX_END` is a
half-open range of faces along the side (`0 9999` = the whole side). Types:
`WALL`, `WALL_FIXED_VALUE`, `INLET`, `INLET_PARABOLIC`, `OUTLET`
(`SYMMETRY` is parsed but acts as a no-slip wall). The pressure condition
follows from the type. Meaning, allowed faces and implementation of each:
[theory/boundary_conditions.md](theory/boundary_conditions.md).

## `boundaries_T.txt` — temperature boundary conditions (DHC)

Same format. Types `WALL_ISOTHERMAL` ($\theta =$ `VALUE`) and
`WALL_ADIABATIC` ($\partial\theta/\partial n = 0$). Sides not listed are
adiabatic. Non-dimensional: hot wall 1, cold wall 0.

## `obstacle.txt` — solid obstacle (cylinder)

| Key | Meaning |
|---|---|
| `obstacle_x0`, `obstacle_x1` | $x$ extent of the solid block |
| `obstacle_y0`, `obstacle_y1` | $y$ extent |

The bounds are snapped to the nearest cell faces (a warning is printed if
they move by more than $10^{-6}$). The obstacle side $D = x_1 - x_0$ after
snapping is the length scale of $Re$, $C_D$, $C_L$ and $St$.

## `reynolds_cases.txt` / `rayleigh_cases.txt` — the sweep

One or more values separated by spaces, commas or new lines; each value is a
separate case, started from rest. Reynolds numbers are integers; Rayleigh
numbers may use exponents (`1e5`). Empty or missing file → default sweep:
$Re = 100, 200, \dots, 1000$ or $Ra = 10^3, 10^4, 10^5, 10^6$.

## `ghia_reference.txt` (LDC)

A readable copy of Tables I and II of Ghia et al. (1982). The code does
**not** read it: the same tables are embedded in
[src/main.cpp](../src/main.cpp) (console printout) and
[plotters/plot_results.py](../plotters/plot_results.py) (figures and errors).
