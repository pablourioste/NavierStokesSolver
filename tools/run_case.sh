#!/bin/bash
# tools/run_case.sh — run ONE case in an isolated working copy and store its
# results in output/<scenario>/, without touching the committed input/ files.
#
# Usage:
#   tools/run_case.sh SCENARIO VALUE [options]
#
#   SCENARIO   lid_driven_cavity | differentially_heated_cavity |
#              square_cylinder | square_cylinder_breuer
#   VALUE      Reynolds number (Rayleigh number for the heated cavity), e.g. 100, 1e6
#
# Options:
#   --set KEY=VALUE   override an entry of input/SCENARIO/sim_params.txt
#                     (dt, max_steps, ss_tol, ...); repeatable
#   --stretch KEY=VALUE  override an entry of stretch.txt (cell counts, betas);
#                     repeatable
#   --scale F         multiply every cell count in stretch.txt by F (grid studies)
#   --threads N       OpenMP threads for the pressure solver (default 4)
#   --store SUBDIR    store under output/SCENARIO/SUBDIR/ instead of
#                     output/SCENARIO/ (e.g. grid_study/D20)
#
# Stored per case, in output/SCENARIO[/SUBDIR]/<Re_N|Ra_N>/:
#   fields/ data/          solver output (VTK, CSV)
#   data/run.log           solver log (pressure-solver warnings filtered out)
#   input/                 the exact input files used for the run
#
# Several cases can run at the same time (each has its own working copy);
# keep the sum of --threads below the number of cores.
set -euo pipefail

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
[ $# -ge 2 ] || { sed -n '2,29p' "$0"; exit 1; }
SCEN="$1"; VAL="$2"; shift 2

SETS=(); STRETCH=(); SCALE=""; THREADS=4; STORE=""
while [ $# -gt 0 ]; do
    case "$1" in
        --set)     SETS+=("$2"); shift 2 ;;
        --stretch) STRETCH+=("$2"); shift 2 ;;
        --scale)   SCALE="$2"; shift 2 ;;
        --threads) THREADS="$2"; shift 2 ;;
        --store)   STORE="$2"; shift 2 ;;
        *) echo "Unknown option: $1" >&2; exit 1 ;;
    esac
done

[ -d "$ROOT/input/$SCEN" ] || { echo "No input/$SCEN/" >&2; exit 1; }
[ -x "$ROOT/build/mesh_gen" ] || { echo "Build first: bash run.sh build" >&2; exit 1; }

if [ "$SCEN" = differentially_heated_cavity ]; then
    CASES_FILE=rayleigh_cases.txt; PREFIX=Ra
else
    CASES_FILE=reynolds_cases.txt; PREFIX=Re
fi

WORK="$(mktemp -d "${TMPDIR:-/tmp}/ns_${SCEN}_${VAL}_XXXX")"
cp -r "$ROOT/input" "$WORK/"
cp "$ROOT/build/mesh_gen" "$WORK/"
IN="$WORK/input/$SCEN"
echo "scenario = $SCEN" > "$WORK/input/scenario.txt"
echo "$VAL" > "$IN/$CASES_FILE"

override() {   # override FILE KEY=VALUE
    local file="$1" key="${2%%=*}" val="${2#*=}"
    if grep -qE "^${key}[[:space:]]*=" "$file"; then
        sed -i -E "s|^${key}[[:space:]]*=.*|${key} = ${val}|" "$file"
    else
        echo "${key} = ${val}" >> "$file"
    fi
}
for kv in "${SETS[@]+"${SETS[@]}"}";       do override "$IN/sim_params.txt" "$kv"; done
for kv in "${STRETCH[@]+"${STRETCH[@]}"}"; do override "$IN/stretch.txt"    "$kv"; done

if [ -n "$SCALE" ]; then
    awk -v f="$SCALE" -F'=' '/^(N|M)_cells_/ { printf "%s= %d\n", $1, int($2 * f + 0.5); next } { print }' \
        "$IN/stretch.txt" > "$IN/stretch.tmp" && mv "$IN/stretch.tmp" "$IN/stretch.txt"
fi

echo "[$(date '+%F %T')] $SCEN $PREFIX=$VAL  (threads=$THREADS${SCALE:+, grid x$SCALE}${STORE:+, store=$STORE})  in $WORK"
set +e
( cd "$WORK" && OMP_NUM_THREADS="$THREADS" ./mesh_gen 2>&1 | grep -av --line-buffered "PCG WARNING" > "$WORK/run.log" )
set -e

CASE="$(ls "$WORK/output/$SCEN" 2>/dev/null | grep -E "^${PREFIX}_" | head -1 || true)"
if [ -z "$CASE" ] || [ ! -d "$WORK/output/$SCEN/$CASE/fields" ]; then
    echo "[$(date '+%F %T')] FAILED: $SCEN $PREFIX=$VAL produced no results. Log kept in $WORK/run.log" >&2
    grep -a "ERROR" "$WORK/run.log" >&2 || true
    exit 1
fi

DEST="$ROOT/output/$SCEN${STORE:+/$STORE}"
mkdir -p "$DEST"
rm -rf "${DEST:?}/$CASE"
cp -r "$WORK/output/$SCEN/$CASE" "$DEST/"
cp "$WORK/run.log" "$DEST/$CASE/data/run.log"
cp -r "$IN" "$DEST/$CASE/input"
[ -d "$DEST/mesh" ] || cp -r "$WORK/output/$SCEN/mesh" "$DEST/"
rm -rf "$WORK"

echo "[$(date '+%F %T')] done: ${DEST#$ROOT/}/$CASE  ($(grep -ao 'Steady state reached.*' "$DEST/$CASE/data/run.log" || echo "ran to max_steps"))"
