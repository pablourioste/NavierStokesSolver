#!/bin/bash
# run.sh — Build and execute the NavierStokesSolver mesh generator.
# Usage: bash run.sh  (from the NavierStokesSolver/ root)

set -e

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
PROJECT_DIR="$SCRIPT_DIR/Mesh_creation"
cd "$PROJECT_DIR"

echo "=== Configuring CMake ==="
cmake -S . -B build_cmake -DCMAKE_BUILD_TYPE=Release

echo "=== Building ==="
cmake --build build_cmake

echo "=== Running mesh_gen ==="
./build/mesh_gen
