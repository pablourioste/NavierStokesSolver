#!/bin/bash
# run.sh — Build and execute the NavierStokesSolver mesh generator.
# Usage: bash run.sh
#
# Must be run from the Mesh_creation/ directory so that relative paths
# (input/, output/) resolve correctly at runtime.

set -e  # Exit immediately on any error

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
cd "$SCRIPT_DIR"

# ---------------------------------------------------------------------------
# 1. Configure + Build with CMake
# ---------------------------------------------------------------------------
echo "=== Configuring CMake ==="
cmake -S . -B build_cmake -DCMAKE_BUILD_TYPE=Release

echo "=== Building ==="
cmake --build build_cmake

# ---------------------------------------------------------------------------
# 2. Run
# ---------------------------------------------------------------------------
echo "=== Running mesh_gen ==="
./build/mesh_gen
