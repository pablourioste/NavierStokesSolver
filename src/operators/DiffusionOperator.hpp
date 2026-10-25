#pragma once
#include "mesh/Mesh_config.h"
#include "mesh/ObstacleMask.hpp"
#include "solver/sparse/solver_class_sparse.hpp"  // Vector

// =============================================================================
// DiffusionOperator
//
// Discretises the diffusion terms  ν∇²u  and  ν∇²v  of the momentum
// equations on the staggered U and V meshes, using 2nd-order finite
// differences on a non-uniform mesh.
//
// For an interior U node (k=2..N, j=1..M):
//
//   d²u/dx² = 2/(dx_e + dx_w) * [(u_E - u_P)/dx_e − (u_P - u_W)/dx_w]
//   d²u/dy² = 2/(dy_n + dy_s) * [(u_N - u_P)/dy_n − (u_P - u_S)/dy_s]
//   result = ν * (d²u/dx² + d²u/dy²)
//
// At the walls (j=1 or j=M) the image method (ghost cell) is applied:
//   u_ghost = 2·u_wall − u_P    →  (u_ghost + u_P)/2 = u_wall  ✓
//   dist    = 2·(y_P − y_wall)   →  correct distance to the boundary face
//
// The ghost-cell wall values must be passed to the constructor.
// For no-slip walls (velocity = 0): ghost = −u_P → pass u_wall = 0.
// For the lid (NORTH, u_wall = U_lid): ghost = 2·U_lid − u_P.
//
// NOTE: this operator contains only the DIFFUSIVE term. Convection stays in
// MomentumSolver.
// =============================================================================
class DiffusionOperator {
public:
    // nu                 : kinematic viscosity
    // u_ghost_north/south: prescribed u value on the N/S wall
    // v_ghost_west/east  : prescribed v value on the W/E wall
    DiffusionOperator(const GridData& pMesh,
                      const GridData& uMesh,
                      const GridData& vMesh,
                      double nu,
                      double u_ghost_north,
                      double u_ghost_south,
                      double v_ghost_west,
                      double v_ghost_east);

    // ν∇²u at every U node. Only interior nodes (k=2..N) are computed.
    // Boundary nodes (k=1, k=N+1) are set to 0.0 in the result vector.
    void applyU(const Vector& u, Vector& result) const;

    // ν∇²v at every V node. Only interior nodes (l=2..M) are computed.
    // Boundary nodes (l=1, l=M+1) are set to 0.0 in the result vector.
    void applyV(const Vector& v, Vector& result) const;

    // Per-node kernels — public so that MomentumSolver can delegate the
    // diffusive part of uRHS/vRHS without duplicating the formula.
    double diffuseU(const Vector& u, int k, int j) const;
    double diffuseV(const Vector& v, int i, int l) const;

private:
    const GridData& pMesh_;
    const GridData& uMesh_;
    const GridData& vMesh_;
    double nu_;
    double u_ghost_north_, u_ghost_south_;
    double v_ghost_west_,  v_ghost_east_;

    inline int uIdx(int k, int j) const {
        return (j - 1) * uMesh_.N_cells_x + (k - 1);
    }
    inline int vIdx(int i, int l) const {
        return (l - 1) * vMesh_.N_cells_x + (i - 1);
    }

};
