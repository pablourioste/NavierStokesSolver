#!/bin/bash
# tools/status.sh — progress of the cases started by tools/run_case.sh
#   bash tools/status.sh            # one snapshot
#   watch -n 30 bash tools/status.sh   # refresh every 30 s (Ctrl-C to quit)
ROOT="$(cd "$(dirname "$0")/.." && pwd)"

printf "%-30s %-10s %8s %10s %7s   %s\n" "CASE" "VALUE" "STEP" "MAX" "%" "LAST RESIDUAL ||du||/dt"
found=0
for d in "${TMPDIR:-/tmp}"/ns_*; do
    [ -f "$d/run.log" ] || continue
    found=1
    scen=$(sed -n 's/^scenario = //p' "$d/input/scenario.txt")
    cases=reynolds_cases.txt; [ "$scen" = differentially_heated_cavity ] && cases=rayleigh_cases.txt
    val=$(head -1 "$d/input/$scen/$cases")
    max=$(sed -n 's/^max_steps *= *//p' "$d/input/$scen/sim_params.txt")
    line=$(grep -a "Step" "$d/run.log" | tail -1)
    step=$(echo "$line" | sed -n 's/.*Step \([0-9]*\).*/\1/p')
    res=$(echo "$line" | sed -n 's/.*= *//p')
    # cylinder cases: the force history is flushed more often than the log
    f=$(ls "$d"/output/"$scen"/*/data/forces.csv 2>/dev/null | head -1)
    if [ -n "$f" ]; then
        fs=$(tail -1 "$f" | cut -d, -f1)
        [[ "$fs" =~ ^[0-9]+$ ]] && [ "$fs" -gt "${step:-0}" ] && step=$fs
    fi
    pct=$(awk -v s="${step:-0}" -v m="${max:-1}" 'BEGIN { printf "%.0f", 100 * s / m }')
    state=""
    grep -aq "Steady state reached" "$d/run.log" && state="  (steady)"
    grep -aq "ERROR" "$d/run.log" && state="  (DIVERGED)"
    printf "%-30s %-10s %8s %10s %6s%%   %s%s\n" "$scen" "$val" "${step:--}" "$max" "$pct" "${res:--}" "$state"
done
[ $found = 1 ] || echo "(no case running)"

echo
echo "Finished / failed:"
grep -ah "done:\|FAILED" "$ROOT"/output/logs/*.log 2>/dev/null | tr -d '\000' |
while IFS= read -r l; do
    # a FAILED line whose kept log no longer exists belongs to a run that was
    # deliberately stopped and cleaned up -> not a real failure
    if [[ "$l" == *FAILED* ]]; then
        k=$(echo "$l" | sed -n 's/.*Log kept in \(.*\)\/run.log.*/\1/p')
        [ -n "$k" ] && [ ! -d "$k" ] && continue
    fi
    echo "$l"
done
pgrep -f run_pending_cases.sh > /dev/null && echo "run_pending_cases.sh: running" \
                                          || echo "run_pending_cases.sh: not running"
