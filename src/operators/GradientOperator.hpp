#pragma once
#include "mesh/Mesh_config.h"
#include "solver/sparse/solver_class_sparse.hpp"  // Vector

// =============================================================================
// GradientOperator
//
// Discretises the pressure gradient  ∇p  from the P mesh onto the U/V meshes,
// using 2nd-order central finite differences on the staggered grid.
//
// For an interior U node (k=2..N, j=1..M):
//   (∂p/∂x)|_{k,j} = (p[pIdx(k,j)] − p[pIdx(k-1,j)]) / (x_P[k] − x_P[k-1])
//   where x_P[k] is the centre of the P cell with logical index k.
//
// For an interior V node (i=1..N, l=2..M):
//   (∂p/∂y)|_{i,l} = (p[pIdx(i,l)] − p[pIdx(i,l-1)]) / (y_P[l] − y_P[l-1])
//
// IMPORTANT: this operator is purely geometric, with no dt/ρ factors.
// The velocity correction in PressurePoissonSystem::projectVelocity()
// applies the additional −(dt/ρ) factor.
//
// Boundary nodes (k=1, k=N+1, l=1, l=M+1) are set to 0.0 in the output.
// =============================================================================
class GradientOperator {
public:
    GradientOperator(const GridData& pMesh,
                     const GridData& uMesh,
                     const GridData& vMesh);

    // Computes the discrete gradient of p.
    // Input:    p      — pressure field at P nodes, size pMesh.num_active_nodes
    // Outputs:  grad_u — ∂p/∂x at interior U nodes, size uMesh.num_active_nodes
    //           grad_v — ∂p/∂y at interior V nodes, size vMesh.num_active_nodes
    void apply(const Vector& p, Vector& grad_u, Vector& grad_v) const;

private:
    const GridData& pMesh_;
    const GridData& uMesh_;
    const GridData& vMesh_;

    inline int pIdx(int i, int j) const {
        return (j - 1) * pMesh_.N_cells_x + (i - 1);
    }
    inline int uIdx(int k, int j) const {
        return (j - 1) * uMesh_.N_cells_x + (k - 1);
    }
    inline int vIdx(int i, int l) const {
        return (l - 1) * vMesh_.N_cells_x + (i - 1);
    }
};
