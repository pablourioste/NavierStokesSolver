#pragma once
#include "mesh/Mesh_config.h"
#include "solver/sparse/solver_class_sparse.hpp"  // Vector alias

// =============================================================================
// addBuoyancy
//
// Boussinesq buoyancy coupling for the Differentially Heated Cavity, applied as
// an external explicit increment to the momentum predictor v* — this keeps
// MomentumSolver and PressurePoissonSystem completely untouched (modular).
//
// In the Ra-Pr (de Vahl Davis) scaling the vertical momentum source is
//   f_v = Ra * Pr * theta
// added to the predictor before the pressure projection, which then enforces
// incompressibility including the buoyant forcing.
//
// theta lives on P cells; the V-node (i,l) sits on the face between P cells
// (i,l-1) and (i,l), so theta is linearly interpolated there:
//   theta_v = 0.5 * (T[pIdx(i,l-1)] + T[pIdx(i,l)])
//
// Only interior V-nodes (l = 2..M) are updated; wall V-nodes (no penetration)
// are left untouched.
// =============================================================================
inline void addBuoyancy(const Vector& T, Vector& v, double dt,
                        double Ra, double Pr,
                        const GridData& pMesh, const GridData& vMesh)
{
    const int N = pMesh.N_cells_x;
    const int M = pMesh.M_cells_y;

    auto pIdx = [&](int i, int j) { return (j - 1) * pMesh.N_cells_x + (i - 1); };
    auto vIdx = [&](int i, int l) { return (l - 1) * vMesh.N_cells_x + (i - 1); };

    const double coeff = dt * Ra * Pr;
    for (int l = 2; l <= M; ++l)
        for (int i = 1; i <= N; ++i) {
            const double theta_v = 0.5 * (T[pIdx(i, l - 1)] + T[pIdx(i, l)]);
            v[vIdx(i, l)] += coeff * theta_v;
        }
}
