#pragma once
#include "mesh/Mesh_config.h"
#include "solver/sparse/solver_class_sparse.hpp"  // Vector alias
#include "operators/DiffusionOperator.hpp"

// Convection discretisation scheme
enum class ConvectionScheme { CDS, UPWIND };

// =============================================================================
// MomentumSolver
//
// Implements the predictor step of Chorin's fractional-step method:
//
//   u* = u^n + dt * [ ν∇²u  −  (u·∇)u ]
//   v* = v^n + dt * [ ν∇²v  −  (u·∇)v ]
//
// Discretisation:
//   - Diffusion : 2nd-order finite differences on non-uniform grid
//   - Convection: CDS (central) or Upwind (1st-order), selectable at construction
//   - Time      : explicit Euler
//
// Grid conventions (1-based logical indices throughout):
//   U-nodes  (k,j)  k=1..N+1, j=1..M   — uIdx(k,j) = (j-1)*(N+1)+(k-1)
//   V-nodes  (i,l)  i=1..N,   l=1..M+1 — vIdx(i,l) = (l-1)*N    +(i-1)
//
// Boundary nodes in the U vector  : k=1 (west), k=N+1 (east)
// Boundary nodes in the V vector  : l=1 (south), l=M+1 (north)
// Ghost values (not in vectors)   : u at top/bottom walls, v at left/right walls
//
// BCs are read from pMesh.bound_* (populated by MeshConfig::define_boundaries).
//   WALL        → ghost/boundary value = 0
//   WALL_FIXED  → ghost/boundary value = BoundaryFace::value
//   INLET       → boundary value = BoundaryFace::value
//   OUTLET      → boundary value = 0 (zero-gradient approx.)
// =============================================================================
class MomentumSolver {
public:
    // nu    : kinematic viscosity
    // scheme: convection scheme (default CDS)
    // Ghost values for the stencil are derived automatically from pMesh.bound_*
    MomentumSolver(const GridData& pMesh,
                   const GridData& uMesh,
                   const GridData& vMesh,
                   double nu,
                   ConvectionScheme scheme = ConvectionScheme::CDS);

    // Apply boundary conditions to u and v vectors from pMesh.bound_* bocos.
    // Must be called at least once before advance(). Re-call after projection.
    void applyBoundaryConditions(Vector& u, Vector& v) const;

    // Advance one time step using Adams-Bashforth 2 (AB2):
    //   u* = u^n + dt * [ (3/2) R(u^n)  −  (1/2) R(u^{n-1}) ]
    // The first call falls back to explicit Euler (no u^{n-1} yet).
    // Interior nodes updated; boundary nodes left unchanged.
    void advance(Vector& u, Vector& v, double dt);

private:
    const GridData& pMesh_;
    const GridData& uMesh_;
    const GridData& vMesh_;
    double nu_;
    ConvectionScheme scheme_;

    // Ghost cell values derived from pMesh_.bound_* at construction time.
    // Used inside uRHS/vRHS stencils for wall-adjacent interior nodes.
    double u_ghost_north_;  // u above j=M row  (north wall)
    double u_ghost_south_;  // u below j=1 row  (south wall)
    double v_ghost_west_;   // v left  of i=1 col (west wall)
    double v_ghost_east_;   // v right of i=N col (east wall)

    // Diffusion operator: delegates ν∇²u and ν∇²v to DiffusionOperator.
    // Declared after ghost doubles so member init order is correct.
    DiffusionOperator diffOp_;

    // Adams-Bashforth 2 state: RHS from the previous time step.
    // Empty until the first call to advance() (first step uses Euler).
    Vector rhs_u_old_;
    Vector rhs_v_old_;
    bool   first_step_ = true;

    // Extract the ghost/boundary value implied by a BCType + stored face value.
    static double ghostValue(BCType type, double faceVal);
    // Extracts the ghost value from a boundary face list (returns 0 if empty).
    static double extractGhost(const std::vector<BoundaryFace>& bound);

    // Index helpers (consistent with PressurePoissonSystem)
    inline int uIdx(int k, int j) const { return (j-1)*uMesh_.N_cells_x + (k-1); }
    inline int vIdx(int i, int l) const { return (l-1)*vMesh_.N_cells_x + (i-1); }

    // RHS of u-momentum (ν∇²u − conv_u) at interior U-node (k,j), k=2..N, j=1..M
    double uRHS(const Vector& u, const Vector& v, int k, int j) const;

    // RHS of v-momentum (ν∇²v − conv_v) at interior V-node (i,l), i=1..N, l=2..M
    double vRHS(const Vector& u, const Vector& v, int i, int l) const;
};
