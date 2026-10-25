#pragma once
#include "mesh/Mesh_config.h"
#include "solver/sparse/solver_class_sparse.hpp"  // Vector

// =============================================================================
// DivergenceOperator
//
// Discretises the divergence  ∇·u  by integrating the net flux over each P
// cell with the finite-volume method (FVM).
//
// For the P cell with logical indices (i,j):
//
//   flux = u_E * Se  −  u_W * Sw  +  v_N * Sn  −  v_S * Ss
//
//   Outward-normal convention:
//     East face  (n̂=+x): +u_E * Se
//     West face  (n̂=−x): −u_W * Sw
//     North face (n̂=+y): +v_N * Sn
//     South face (n̂=−y): −v_S * Ss
//
//   Staggered index mapping for P(i,j):
//     u_E = u[uIdx(i+1, j)],  u_W = u[uIdx(i,   j)]
//     v_N = v[vIdx(i,   j+1)], v_S = v[vIdx(i,   j)]
//
// The result is divided by the cell volume V[i][j], giving the divergence
// per unit volume (units of [1/time] if u has units of velocity). This
// allows a direct comparison with the analytical ∇·u in the MMS tests.
//
// IMPORTANT: no ρ/dt factors — purely geometric operator.
// =============================================================================
class DivergenceOperator {
public:
    DivergenceOperator(const GridData& pMesh,
                       const GridData& uMesh,
                       const GridData& vMesh);

    // Computes the discrete ∇·u per unit volume in every P cell.
    // Input:  u (size uMesh.num_active_nodes), v (size vMesh.num_active_nodes)
    // Output: vector of size pMesh.num_active_nodes
    Vector apply(const Vector& u, const Vector& v) const;

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
