#pragma once
#include "mesh/Mesh_config.h"
#include "solver/sparse/solver_class_sparse.hpp"  // Vector
#include <utility>
#include <vector>

// =============================================================================
// ObstacleForces
//
// Surface-integrates pressure and viscous shear stress over the 4 flat faces
// of a rectangular immersed obstacle (Square Cylinder scenario), producing
// the net force (Fx, Fy) on the obstacle, and from it drag/lift coefficients.
//
// Pressure contribution: for each face, sum p_fluid * Area * n (outward
// normal) over the fluid-side P-cells immediately adjacent to that face —
// NEVER the solid cell's own (pinned-to-0, physically meaningless) pressure.
//
// Viscous contribution: tau = rho * nu * (d(u_tangential)/dn) evaluated at
// the face using a one-sided finite difference between the fluid-side
// velocity node and the (zero, no-slip) obstacle-surface node — the same
// wall-normal distance already used by the ghost-image stencils in
// MomentumSolver/DiffusionOperator.
// =============================================================================
class ObstacleForces {
public:
    ObstacleForces(const GridData& pMesh, const GridData& uMesh, const GridData& vMesh,
                   double nu, double rho);

    // Returns (Fx, Fy): net force on the obstacle from the current p, u, v
    // fields (pressure + viscous, integrated over the 4 obstacle faces).
    std::pair<double, double> computeForces(const Vector& p, const Vector& u, const Vector& v) const;

private:
    const GridData& pMesh_;
    const GridData& uMesh_;
    const GridData& vMesh_;
    double nu_;
    double rho_;

    inline int pIdx(int i, int j) const { return (j - 1) * pMesh_.N_cells_x + (i - 1); }
    inline int uIdx(int k, int j) const { return (j - 1) * uMesh_.N_cells_x + (k - 1); }
    inline int vIdx(int i, int l) const { return (l - 1) * vMesh_.N_cells_x + (i - 1); }
};

// Peak-detection Strouhal estimate on a (mostly) periodic Cl(t) time series.
// Discards the first half of the series as the startup transient, then
// averages the period between consecutive positive-going zero-crossings.
// Returns the dominant shedding frequency (1 / time unit), or 0.0 if fewer
// than 2 crossings are found after the transient (insufficient data).
double estimate_shedding_frequency(const std::vector<double>& time,
                                    const std::vector<double>& Cl);
