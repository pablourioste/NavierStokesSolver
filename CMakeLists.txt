cmake_minimum_required(VERSION 3.15)
project(NavierStokesSolver CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# ---------------------------------------------------------------------------
# Build output goes to build/ subfolder, keeping root clean
# ---------------------------------------------------------------------------
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/build)

# ---------------------------------------------------------------------------
# Sources
# ---------------------------------------------------------------------------
set(MESH_SOURCES
    src/mesh/MeshGenerator.cpp
)

set(SOLVER_SOURCES
    src/PressurePoissonSystem.cpp
    src/solver/sparse/solver_sparse.cpp
    src/MomentumSolver.cpp
    src/EnergySolver.cpp
)

set(MAIN_SOURCE
    src/main.cpp
)

set(OPERATOR_SOURCES
    src/operators/LaplaceOperator.cpp
    src/operators/GradientOperator.cpp
    src/operators/DivergenceOperator.cpp
    src/operators/DiffusionOperator.cpp
)

# ---------------------------------------------------------------------------
# Executable: mesh_gen
# ---------------------------------------------------------------------------
add_executable(mesh_gen
    ${MAIN_SOURCE}
    ${MESH_SOURCES}
    ${SOLVER_SOURCES}
    ${OPERATOR_SOURCES}
)

target_include_directories(mesh_gen PRIVATE
    ${CMAKE_SOURCE_DIR}/src
    ${CMAKE_SOURCE_DIR}/src/mesh
    ${CMAKE_SOURCE_DIR}/src/solver
    ${CMAKE_SOURCE_DIR}/src/operators
)

target_compile_options(mesh_gen PRIVATE -Wall -Wextra)

# ---------------------------------------------------------------------------
# Executable: reynolds_sweep  (batería de casos Re)
# ---------------------------------------------------------------------------
add_executable(reynolds_sweep
    src/reynolds_sweep.cpp
    ${MESH_SOURCES}
    ${SOLVER_SOURCES}
    ${OPERATOR_SOURCES}
)

target_include_directories(reynolds_sweep PRIVATE
    ${CMAKE_SOURCE_DIR}/src
    ${CMAKE_SOURCE_DIR}/src/mesh
    ${CMAKE_SOURCE_DIR}/src/solver
    ${CMAKE_SOURCE_DIR}/src/operators
)

target_compile_options(reynolds_sweep PRIVATE -Wall -Wextra)

# ---------------------------------------------------------------------------
# Executable: operator_verify  (verificación MMS independiente)
# ---------------------------------------------------------------------------
add_executable(operator_verify
    src/verification/verify_operators.cpp
    ${MESH_SOURCES}
    ${SOLVER_SOURCES}
    ${OPERATOR_SOURCES}
)

target_include_directories(operator_verify PRIVATE
    ${CMAKE_SOURCE_DIR}/src
    ${CMAKE_SOURCE_DIR}/src/mesh
    ${CMAKE_SOURCE_DIR}/src/solver
    ${CMAKE_SOURCE_DIR}/src/operators
)

target_compile_options(operator_verify PRIVATE -Wall -Wextra -O2)

# ---------------------------------------------------------------------------
# Convenience message
# ---------------------------------------------------------------------------
message(STATUS "Source dir : ${CMAKE_SOURCE_DIR}")
message(STATUS "Build dir  : ${CMAKE_BINARY_DIR}")
message(STATUS "Executable : ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/mesh_gen")
