#ifndef PRESSURE_POISSON_SYSTEM_HPP
#define PRESSURE_POISSON_SYSTEM_HPP

#include "mesh/Mesh_config.h"
#include <vector>
#include <stdexcept>

// =============================================================================
// Triplet — one non-zero entry (row, col, value) of the sparse matrix A.
//
// The layout is intentionally compatible with Eigen::Triplet<double> so that,
// when a linear-algebra library is added, you can pass this vector directly:
//
//   std::vector<Eigen::Triplet<double>> eigen_trips(trips.begin(), trips.end());
//   A.setFromTriplets(eigen_trips.begin(), eigen_trips.end());
// =============================================================================
struct Triplet {
    int    row, col;
    double value;
};


// =============================================================================
// PressurePoissonSystem
//
// Implements the Pressure-Correction step of the fractional-step / projection
// method for incompressible 2D Navier-Stokes on a staggered MAC grid.
//
// Workflow per time step:
//   1.  (Momentum solver, external) compute predictor fields u*, v*.
//   2.  system.setVelocityFields(u_star, v_star);
//   3.  auto trips = system.assembleMatrix();
//       auto b = system.computeRHS(dt, rho);
//       system.pinReferencePressure(trips, b);  // only if all-Neumann BCs
//                                               // (also sets b[0] = 0 internally)
//   4.  <solve  A·p = b  with your sparse solver>
//   5.  std::vector<double> u_next, v_next;
//       system.projectVelocity(p_solution, dt, rho, u_next, v_next);
//   6.  u_next / v_next are now the divergence-free velocity fields.
//
// Grid layout (1-based logical indices):
//
//   pMesh  — P cells    (i=1..N , j=1..M )   centers at (x[i], y[j])
//   uMesh  — U cells    (k=1..N+1, j=1..M )  k=1 west wall, k=N+1 east wall
//   vMesh  — V cells    (i=1..N , l=1..M+1)  l=1 south wall, l=M+1 north wall
//
// U staggered mapping:
//   uMesh.xvc[k] = pMesh.x[k]  for k=1..N  (U cell boundaries = P centers)
//   → U node k lies between P cells (k-1) and (k) in x.
//   → East face of P cell (i,j) ↔ U node k_u = i+1
//   → West face of P cell (i,j) ↔ U node k_u = i
//
// V staggered mapping (symmetric argument in y):
//   vMesh.yvc[l] = pMesh.y[l]  for l=1..M
//   → North face of P cell (i,j) ↔ V node l_v = j+1
//   → South face of P cell (i,j) ↔ V node l_v = j
// =============================================================================
class PressurePoissonSystem {
public:

    // -------------------------------------------------------------------------
    // Constructor — stores const references to the three meshes.
    // The mesh objects must remain valid for the lifetime of this object.
    // -------------------------------------------------------------------------
    PressurePoissonSystem(const GridData& pMesh,
                          const GridData& uMesh,
                          const GridData& vMesh);

    // -------------------------------------------------------------------------
    // setVelocityFields
    //
    // Supply the predictor velocity fields u*, v* computed by the momentum
    // solver.  Sizes must match mesh.num_active_nodes exactly.
    //   u_star[global_id]  — for every u-cell (k,j), k=1..N+1, j=1..M
    //   v_star[global_id]  — for every v-cell (i,l), i=1..N,   l=1..M+1
    //
    // Boundary nodes (walls, inlets) must carry their prescribed values so
    // that computeRHS() picks up the correct boundary flux automatically.
    // -------------------------------------------------------------------------
    void setVelocityFields(const std::vector<double>& u_star,
                           const std::vector<double>& v_star);

    // -------------------------------------------------------------------------
    // assembleMatrix
    //
    // Builds the discrete Laplacian A for pressure using FVM.
    // For each P cell P:
    //   a_nb = Area_face / dist(x_P, x_nb)   (off-diagonal, positive)
    //   a_P  = -sum(a_nb over active neighbors)  (diagonal, negative)
    //
    // Boundary faces obey Homogeneous Neumann (dp/dn = 0) via the ghost-cell
    // method: the boundary coefficient is set to 0, and the central coefficient
    // sums only active (interior) neighbors.
    //
    // NOTE: if ALL boundary faces are Neumann the matrix is singular (null
    // space = constant pressure).  Call pinReferencePressure() to regularize.
    // If any boundary has a Dirichlet condition (OUTLET), modify the relevant
    // row here or in a separate post-processing step.
    // -------------------------------------------------------------------------
    std::vector<Triplet> assembleMatrix() const;

    // -------------------------------------------------------------------------
    // pinReferencePressure
    //
    // Modifies the triplet list in place to enforce p[0] = 0 (row replacement):
    //   - All off-diagonal entries in row 0 are zeroed.
    //   - The diagonal of row 0 is set to 1.
    // Also sets b[0] = 0.0 in the RHS vector to match the pinned row,
    // preventing mathematical corruption if the caller forgets to do so.
    //
    // Only needed when no Dirichlet pressure boundary is present.
    // -------------------------------------------------------------------------
    void pinReferencePressure(std::vector<Triplet>& triplets,
                              std::vector<double>& b) const;

    // -------------------------------------------------------------------------
    // computeRHS
    //
    // Computes b_P = (rho/dt) * div(u*)_P for every P cell, integrated with FVM:
    //
    //   b_P = (rho/dt) * [  u*_E * Se - u*_W * Sw
    //                      + v*_N * Sn - v*_S * Ss ]
    //
    // The outward-normal sign convention is:
    //   East/North faces: outward normal = +x / +y  →  +u*·A
    //   West/South faces: outward normal = -x / -y  →  -u*·A
    // -------------------------------------------------------------------------
    std::vector<double> computeRHS(double dt, double rho) const;

    // -------------------------------------------------------------------------
    // projectVelocity
    //
    // Corrects the predictor velocities using the solved pressure field and
    // writes the divergence-free result into u_next and v_next:
    //
    //   u_next = u* - (dt/rho) * (p_E - p_W) / dist_EW   [interior u-nodes]
    //   v_next = v* - (dt/rho) * (p_N - p_S) / dist_NS   [interior v-nodes]
    //
    // Boundary velocity nodes (walls, inlets) are NOT modified; they retain
    // their prescribed values from setVelocityFields().
    //
    // u_next and v_next are resized and overwritten on every call.
    // -------------------------------------------------------------------------
    void projectVelocity(const std::vector<double>& P_sol,
                         double dt, double rho,
                         std::vector<double>& u_next,
                         std::vector<double>& v_next);

private:
    // Mesh references
    const GridData& pMesh_;
    const GridData& uMesh_;
    const GridData& vMesh_;

    // Velocity fields (predictor → corrected after projectVelocity)
    std::vector<double> u_star_;
    std::vector<double> v_star_;

    // -------------------------------------------------------------------------
    // Index helpers — convert 1-based (i,j) / (k,j) / (i,l) logical indices
    // to the 0-based global_id used in cells[] and the solver vectors.
    //
    // All indices are 1-based on input, matching the mesh generation loops.
    // -------------------------------------------------------------------------

    // P cell (i=1..N, j=1..M)  →  k = (j-1)*N + (i-1)
    inline int pIdx(int i, int j) const {
        return (j - 1) * pMesh_.N_cells_x + (i - 1);
    }

    // U cell (k=1..N+1, j=1..M) → id = (j-1)*(N+1) + (k-1)
    inline int uIdx(int k, int j) const {
        return (j - 1) * uMesh_.N_cells_x + (k - 1);
    }

    // V cell (i=1..N, l=1..M+1) → id = (l-1)*N + (i-1)
    inline int vIdx(int i, int l) const {
        return (l - 1) * vMesh_.N_cells_x + (i - 1);
    }
};

#endif // PRESSURE_POISSON_SYSTEM_HPP
