#pragma once
#include "mesh/Mesh_config.h"
#include "PressurePoissonSystem.hpp"           // Triplet
#include "solver/sparse/solver_class_sparse.hpp"  // CSRMatrix, Vector

// =============================================================================
// LaplaceOperator
//
// Discretises the Laplace operator  ∇²p  on the pressure (P) mesh with the
// finite-volume method (FVM).
//
// For each P cell with logical indices (i,j):
//   a_E = Se[i][j] / (x[i+1] - x[i])      ← East face contribution
//   a_W = Sw[i][j] / (x[i]   - x[i-1])    ← West face contribution
//   a_N = Sn[i][j] / (y[j+1] - y[j])      ← North face contribution
//   a_S = Ss[i][j] / (y[j]   - y[j-1])    ← South face contribution
//   a_P = +(a_E + a_W + a_N + a_S)         ← diagonal (positive)
//
// Homogeneous Neumann boundary condition (∂p/∂n = 0) by default: boundary
// faces are skipped (coefficient = 0).
//
// The resulting matrix is SPD → suitable for CG without a preconditioner.
// If every face is Neumann the system is singular: use pinReference().
// =============================================================================
class LaplaceOperator {
public:
    // pMesh must outlive this object.
    explicit LaplaceOperator(const GridData& pMesh);

    // Builds the sparse Laplacian matrix A in triplet format.
    // Equivalent to PressurePoissonSystem::assembleMatrix().
    std::vector<Triplet> assembleToTriplets() const;

    // Builds the matrix directly in CSR format.
    CSRMatrix assemble() const;

    // Matrix-vector product  y = A·p  (without solving the system).
    // Useful for MMS tests: compare (A·φ)/V with the analytical ∇²φ.
    Vector apply(const Vector& p) const;

    // Regularisation for all-Neumann systems: fixes p[0] = 0.
    // Same as PressurePoissonSystem::pinReferencePressure().
    void pinReference(std::vector<Triplet>& triplets, Vector& b) const;

private:
    const GridData& pMesh_;

    inline int pIdx(int i, int j) const {
        return (j - 1) * pMesh_.N_cells_x + (i - 1);
    }
};
