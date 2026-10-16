#!/bin/bash
# tools/run_pending_cases.sh — computes the cases that are still missing from
# output/, stores them there (via tools/run_case.sh) and regenerates every
# figure. Run from anywhere:
#
#   bash run.sh build                                   # once, after any code change
#   nohup bash tools/run_pending_cases.sh > output/logs/pending.log 2>&1 &
#
# Progress:   tail -f output/logs/*.log
# Summary:    grep -h "done\|FAILED" output/logs/*.log
#
# Thread budget: 16 cores. Batch 1 uses 16 threads (cylinder, OpenMP pressure
# solver); batch 2 runs the cavity cases, whose CG solver is serial (1 thread
# each). Approximate wall times on the 16-core development machine.
set -uo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
cd "$ROOT"
mkdir -p output/logs
RUN=tools/run_case.sh
log() { echo "output/logs/$1.log"; }

# ---------------------------------------------------------------------------
# Batch 1 — Square cylinder, Breuer et al. (2000) set-up          (~2.5 h)
#   dt = 1e-3 -> t U/D = 0.008 * steps
# ---------------------------------------------------------------------------
# Re = 100 to t U/D = 400 (was 256): sharper St, more saturated periods
$RUN square_cylinder_breuer 100 --set max_steps=50000 --set ss_tol=0.0 --threads 4 \
    > "$(log breuer_Re100)" 2>&1 &
# Re = 70 to t U/D = 650: the instability grows slowly near onset (Re_cr ~ 60)
$RUN square_cylinder_breuer 70  --set max_steps=81000 --set ss_tol=0.0 --threads 4 \
    > "$(log breuer_Re70)" 2>&1 &
# Re = 150 and 200 to t U/D = 300: Strouhal-number trend above Re = 100
$RUN square_cylinder_breuer 150 --set max_steps=37500 --set ss_tol=0.0 --threads 3 \
    > "$(log breuer_Re150)" 2>&1 &
$RUN square_cylinder_breuer 200 --set max_steps=37500 --set ss_tol=0.0 --threads 3 \
    > "$(log breuer_Re200)" 2>&1 &
# Grid check of St at Re = 100: D/20 cells across the cylinder
$RUN square_cylinder_breuer 100 --scale 1.25 --set max_steps=50000 --set ss_tol=0.0 \
    --threads 2 --store grid_study/D20 > "$(log breuer_Re100_D20)" 2>&1 &
wait

# ---------------------------------------------------------------------------
# Batch 2 — cavities (serial CG solver, 1 thread each)            (~3 h)
# ---------------------------------------------------------------------------
# Lid-driven cavity: every Ghia et al. (1982) Reynolds number on the shipped
# 120 x 120 stretched grid (the stored 50 x 50 runs are under-resolved at
# Re >= 1000). Explicit diffusion limit on this grid: dt < h^2/(8 nu), with
# h_min = 5.6e-3 -> Re = 100 needs dt = 3e-4.
$RUN lid_driven_cavity 100  --set dt=3e-4 --set max_steps=100000 --set ss_tol=1e-5 --threads 1 \
    > "$(log lid_Re100)" 2>&1 &
for re in 400 1000 3200 5000; do
    $RUN lid_driven_cavity $re --set dt=1e-3 --set max_steps=200000 --set ss_tol=1e-5 --threads 1 \
        > "$(log lid_Re$re)" 2>&1 &
done
$RUN lid_driven_cavity 7500 --set dt=1e-3 --set max_steps=250000 --set ss_tol=1e-5 --threads 1 \
    > "$(log lid_Re7500)" 2>&1 &

# Differentially heated cavity, Ra = 1e6 (de Vahl Davis: Nu = 8.800). The
# shipped stretched grid is impractical with the explicit scheme (dt ~ 1e-7),
# so a uniform 100 x 100 grid; dt = 2.5e-6 is needed for the start-up
# transient. 160 000 steps -> t = 0.4 (diffusive units).
$RUN differentially_heated_cavity 1e6 \
    --stretch N_cells_west=100 --stretch M_cells_west=100 \
    --stretch beta_x_west=0.0 --stretch beta_x_east=0.0 \
    --stretch beta_y_south=0.0 --stretch beta_y_north=0.0 \
    --set dt=2.5e-6 --set max_steps=160000 --set ss_tol=0.0 --threads 1 \
    > "$(log dhc_Ra1e6)" 2>&1 &
wait

# ---------------------------------------------------------------------------
# Plot everything (figures go to output/<scenario>/.../figures/)
# ---------------------------------------------------------------------------
PY="$ROOT/.venv/bin/python3"; [ -x "$PY" ] || PY=python3
"$PY" plotters/plot_square_cylinder.py square_cylinder_breuer
"$PY" plotters/plot_square_cylinder.py square_cylinder_breuer/grid_study/D20
"$PY" plotters/plot_results.py
"$PY" plotters/plot_dhc.py
echo "[$(date '+%F %T')] all pending cases finished"
grep -h "done\|FAILED" output/logs/*.log
