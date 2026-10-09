#pragma once
#include "mesh/Mesh_config.h"
#include "solver/sparse/solver_class_sparse.hpp"  // Vector alias
#include "MomentumSolver.hpp"                      // reuse ConvectionScheme enum

// =============================================================================
// EnergySolver
//
// Explicit advection-diffusion solver for the non-dimensional temperature
// field theta, used by the Differentially Heated Cavity (Boussinesq). Mirrors
// the structure of MomentumSolver, but the scalar lives on the PRESSURE grid
// (cell centred), exactly like the pressure field:
//
//   d(theta)/dt + (u . grad) theta = alpha * laplacian(theta)
//
// In the Ra-Pr (de Vahl Davis) scaling the thermal diffusivity alpha = 1.
//
// Why the P-grid is the natural home for theta:
//   - The staggered velocities ARE the face-normal velocities of each P cell,
//     so convective fluxes need no interpolation:
//       east face  u_e = u[uIdx(i+1,j)]   west face  u_w = u[uIdx(i,j)]
//       north face v_n = v[vIdx(i,j+1)]   south face v_s = v[vIdx(i,j)]
//     (same mapping as PressurePoissonSystem::computeRHS).
//   - Diffusion reuses pMesh face areas (Se/Sw/Sn/Ss) and cell centres, with
//     the same ghost-cell image method as DiffusionOperator.
//
// Boundary conditions are read from pMesh.bound_*_T:
//   WALL_ISOTHERMAL -> Dirichlet (theta = value), ghost = 2*theta_wall - theta_P
//   WALL_ADIABATIC  -> zero flux (d theta/dn = 0), ghost = theta_P
//
// All P cells are interior unknowns (walls are ghosts), so there are no
// boundary nodes to set inside the theta vector.
// =============================================================================
class EnergySolver {
public:
    // alpha : thermal diffusivity (= 1 in Ra-Pr scaling)
    // scheme: convection scheme (UPWIND recommended for high-Ra stability)
    EnergySolver(const GridData& pMesh,
                 const GridData& uMesh,
                 const GridData& vMesh,
                 double alpha,
                 ConvectionScheme scheme = ConvectionScheme::UPWIND);

    // Advance theta one step with Adams-Bashforth 2 (Euler on the first call),
    // using the (already corrected, divergence-free) velocity fields u, v.
    void advance(Vector& T, const Vector& u, const Vector& v, double dt);

private:
    const GridData& pMesh_;
    const GridData& uMesh_;
    const GridData& vMesh_;
    double alpha_;
    ConvectionScheme scheme_;

    // Per-wall temperature BC (uniform along each side, like MomentumSolver ghosts).
    BCType  bc_west_,  bc_east_,  bc_north_,  bc_south_;
    double  Tw_west_,  Tw_east_,  Tw_north_,  Tw_south_;

    // AB2 state
    Vector rhs_old_;
    bool   first_step_ = true;

    // Index helpers (1-based logical indices, 0-based storage)
    inline int pIdx(int i, int j) const { return (j-1)*pMesh_.N_cells_x + (i-1); }
    inline int uIdx(int k, int j) const { return (j-1)*uMesh_.N_cells_x + (k-1); }
    inline int vIdx(int i, int l) const { return (l-1)*vMesh_.N_cells_x + (i-1); }

    // RHS at interior P cell (i,j): alpha*lap(theta) - (u.grad)theta
    double tRHS(const Vector& T, const Vector& u, const Vector& v, int i, int j) const;

    // Extract (type,value) of a side's temperature BC from a bound_*_T list.
    static void extractWallT(const std::vector<BoundaryFace>& bound,
                             BCType& type, double& value);
};
