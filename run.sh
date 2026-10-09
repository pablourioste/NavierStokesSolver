#!/bin/bash
# run.sh — Configure, build, and run the NavierStokesSolver.
#
# Usage:
#   bash run.sh              # build + run the main solver (mesh_gen)
#   bash run.sh mesh_gen     # same as above, explicit
#   bash run.sh verify       # build + run the operator convergence check
#   bash run.sh sweep        # build + run the Reynolds/Rayleigh batch sweep
#   bash run.sh build        # configure + build only, run nothing
#
# Always run from the repository root so the relative input/ and output/
# paths resolve correctly.

set -e  # Exit immediately on any error

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
cd "$SCRIPT_DIR"

TARGET="${1:-mesh_gen}"

# ---------------------------------------------------------------------------
# 1. Configure + build with CMake (builds all three executables: mesh_gen,
#    operator_verify, reynolds_sweep — see CMakeLists.txt)
# ---------------------------------------------------------------------------
echo "=== Configuring CMake ==="
cmake -S . -B build_cmake -DCMAKE_BUILD_TYPE=Release

echo "=== Building ==="
cmake --build build_cmake

# ---------------------------------------------------------------------------
# 2. Run the requested executable
# ---------------------------------------------------------------------------
case "$TARGET" in
  build)
    echo "=== Build only, nothing to run ==="
    ;;
  mesh_gen)
    echo "=== Running mesh_gen (main solver: lid-driven cavity / DHC, per input/sim_params.txt) ==="
    ./build/mesh_gen
    ;;
  verify)
    echo "=== Running operator_verify (MMS convergence check for the discrete operators) ==="
    ./build/operator_verify
    ;;
  sweep)
    echo "=== Running reynolds_sweep (batch sweep over input/reynolds_cases.txt) ==="
    ./build/reynolds_sweep
    ;;
  *)
    echo "Unknown target: '$TARGET'"
    echo "Usage: bash run.sh [mesh_gen|verify|sweep|build]"
    exit 1
    ;;
esac
