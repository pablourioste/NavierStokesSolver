#include "PressurePoissonSystem.hpp"
#include "operators/LaplaceOperator.hpp"
#include "operators/DivergenceOperator.hpp"
#include "operators/GradientOperator.hpp"
#include <stdexcept>
#include <cmath>

// =============================================================================
//  STAGGERED GRID REFERENCE  (1-based logical indices throughout)
//
//  pMesh:  N × M cells     (i=1..N , j=1..M )
//    xvc[0..N+1]  — face x-positions   (xvc[0]=west wall, xvc[N]=east wall)
//    yvc[0..M+1]  — face y-positions
//    x[0..N+1]    — cell-center x      (x[i] = centre of P cell i)
//    y[0..M+1]    — cell-center y
//    Se[i][j]     — east face area  = (yvc[j]-yvc[j-1]) * W_depth
//    Sn[i][j]     — north face area = (xvc[i]-xvc[i-1]) * W_depth
//
//  uMesh:  (N+1) × M cells (k=1..N+1, j=1..M)
//    uMesh.xvc[0]   = pMesh.xvc[0]   (west domain wall)
//    uMesh.xvc[k]   = pMesh.x[k]     for k=1..N  (P cell centres)
//    uMesh.xvc[N+1] = pMesh.xvc[N]   (east domain wall)
//    → U node k=1    is the WEST boundary node  (k_u = 1)
//    → U node k=i+1  sits at x = pMesh.xvc[i]   (east face of P cell i)
//    → U node k=N+1  is the EAST boundary node
//
//    Face mapping for P cell (i,j):
//      East face  → U node (k_u = i+1, j)
//      West face  → U node (k_u = i,   j)
//
//  vMesh:  N × (M+1) cells (i=1..N , l=1..M+1)
//    (symmetric argument in y)
//    → V node l=1    is the SOUTH boundary node
//    → V node l=j+1  sits at y = pMesh.yvc[j]   (north face of P cell j)
//    → V node l=M+1  is the NORTH boundary node
//
//    Face mapping for P cell (i,j):
//      North face → V node (i, l_v = j+1)
//      South face → V node (i, l_v = j  )
// =============================================================================


// ---------------------------------------------------------------------------
// Constructor
// ---------------------------------------------------------------------------
PressurePoissonSystem::PressurePoissonSystem(const GridData& pMesh,
                                             const GridData& uMesh,
                                             const GridData& vMesh)
    : pMesh_(pMesh), uMesh_(uMesh), vMesh_(vMesh)
{
    if (pMesh_.Se.empty() || pMesh_.x.empty())
        throw std::runtime_error(
            "PressurePoissonSystem: pMesh face-area (Se/Sw/Sn/Ss) or cell-centre "
            "arrays (x/y) are empty. Run generate_mesh() and build_connectivity() first.");

    u_star_.assign(uMesh_.num_active_nodes, 0.0);
    v_star_.assign(vMesh_.num_active_nodes, 0.0);
}


// ---------------------------------------------------------------------------
// setVelocityFields
// ---------------------------------------------------------------------------
void PressurePoissonSystem::setVelocityFields(const std::vector<double>& u_star,
                                               const std::vector<double>& v_star)
{
    if (static_cast<int>(u_star.size()) != uMesh_.num_active_nodes)
        throw std::runtime_error(
            "PressurePoissonSystem::setVelocityFields: u_star size ("
            + std::to_string(u_star.size()) + ") != uMesh nodes ("
            + std::to_string(uMesh_.num_active_nodes) + ").");

    if (static_cast<int>(v_star.size()) != vMesh_.num_active_nodes)
        throw std::runtime_error(
            "PressurePoissonSystem::setVelocityFields: v_star size ("
            + std::to_string(v_star.size()) + ") != vMesh nodes ("
            + std::to_string(vMesh_.num_active_nodes) + ").");

    u_star_ = u_star;
    v_star_ = v_star;
}


// ---------------------------------------------------------------------------
// assembleMatrix
//
// Builds the sparse Laplacian A such that  A·p = b  (PPE).
//
// For each P cell k with logical index (i,j):
//
//   a_nb  =  Area_face / dist(x_P, x_nb)    [positive, off-diagonal]
//   a_P   = -sum(a_nb over active neighbours) [negative, diagonal]
//
// Homogeneous Neumann (dp/dn = 0) at every boundary face is handled via the
// ghost-cell method: setting p_ghost = p_P makes the face gradient zero, so
// the boundary face contributes 0 both to the off-diagonal AND the diagonal.
// Concretely: just skip boundary faces entirely.
//
// NOTE — Dirichlet pressure at an OUTLET face (e.g. p = 0):
//   Replace the boundary check with a Dirichlet treatment:
//     a_P  +=  a_face          (move known p to RHS instead of zeroing)
//     b[k] -= a_face * p_bc   (handled in computeRHS or a post-step)
//   and do NOT add an off-diagonal entry for that face.
// ---------------------------------------------------------------------------
std::vector<Triplet> PressurePoissonSystem::assembleMatrix() const
{
    return LaplaceOperator(pMesh_).assembleToTriplets();
}


// ---------------------------------------------------------------------------
// pinReferencePressure
//
// Regularises an all-Neumann system by fixing p[0] = 0 (row replacement).
// Zeros every entry in row 0 and column 0 (symmetric elimination), then sets
// the diagonal to 1.  Also sets b[0] = 0 to match the pinned row, preventing
// mathematical corruption if the caller forgets to do so.
// ---------------------------------------------------------------------------
void PressurePoissonSystem::pinReferencePressure(std::vector<Triplet>& triplets,
                                                  std::vector<double>& b) const
{
    LaplaceOperator(pMesh_).pinReference(triplets, b);
}


// ---------------------------------------------------------------------------
// computeRHS
//
// RHS of the PPE:  b_P = (rho/dt) * ∫ div(u*) dV
//
// In FVM integral form for a rectangular cell:
//
//   ∫ div(u*) dV = Σ_faces  (u*·n̂) * A_face
//                = u*_E * Se  −  u*_W * Sw
//                + v*_N * Sn  −  v*_S * Ss
//
// The outward-normal sign:
//   East  (n̂ = +x):  flux = +u*_E * Se
//   West  (n̂ = -x):  flux = -u*_W * Sw   (u* is +x, outward is -x → minus sign)
//   North (n̂ = +y):  flux = +v*_N * Sn
//   South (n̂ = -y):  flux = -v*_S * Ss
//
// Staggered index mapping for P cell (i,j):
//   East  U node → uIdx(i+1, j)     West  U node → uIdx(i,   j)
//   North V node → vIdx(i,   j+1)   South V node → vIdx(i,   j)
// ---------------------------------------------------------------------------
std::vector<double> PressurePoissonSystem::computeRHS(double dt, double rho) const
{
    // DivergenceOperator::apply() returns flux/V per cell (divergence per unit volume).
    // PPE RHS: b[k] = -(ρ/dt) * ∫div(u*)dV = -(ρ/dt) * div_vec[k] * V[i][j]
    const Vector div_vec =
        DivergenceOperator(pMesh_, uMesh_, vMesh_).apply(u_star_, v_star_);

    const int n_cells = pMesh_.num_active_nodes;
    std::vector<double> b(n_cells, 0.0);
    for (int k = 0; k < n_cells; ++k) {
        const int i = pMesh_.id_to_i[k];
        const int j = pMesh_.id_to_j[k];
        b[k] = -(rho / dt) * div_vec[k] * pMesh_.V[i][j];
    }
    return b;
}


// ---------------------------------------------------------------------------
// projectVelocity
//
// Velocity correction step:
//
//   u^{n+1} = u*  −  (dt/ρ) * ∂p/∂x
//   v^{n+1} = v*  −  (dt/ρ) * ∂p/∂y
//
// The discrete pressure gradient at each velocity node is evaluated between
// the two adjacent P cells that straddle that node.
//
// Interior U node (i, j),  i = 2..N  (i=1 and i=N+1 are wall/inlet nodes)
//   → straddles P cells (i-1, j) and (i, j)
//   → ∂p/∂x ≈ (P[pIdx(i,j)] − P[pIdx(i−1,j)]) / (pMesh.x[i] − pMesh.x[i−1])
//
// Interior V node (i, l),  l = 2..M  (l=1 and l=M+1 are wall/inlet nodes)
//   → straddles P cells (i, l-1) and (i, l)
//   → ∂p/∂y ≈ (P[pIdx(i,l)] − P[pIdx(i,l−1)]) / (pMesh.y[l] − pMesh.y[l−1])
//
// Boundary velocity nodes (i=1, i=N+1, l=1, l=M+1) are NOT corrected: they
// hold their prescribed values (wall no-slip or inlet conditions).
//
// Results are written into u_next and v_next (output parameters) so that the
// corrected fields are directly accessible to the calling solver.
// ---------------------------------------------------------------------------
void PressurePoissonSystem::projectVelocity(const std::vector<double>& P_sol,
                                             double dt, double rho,
                                             std::vector<double>& u_next,
                                             std::vector<double>& v_next)
{
    if (static_cast<int>(P_sol.size()) != pMesh_.num_active_nodes)
        throw std::runtime_error(
            "PressurePoissonSystem::projectVelocity: P_sol size ("
            + std::to_string(P_sol.size()) + ") != pMesh nodes ("
            + std::to_string(pMesh_.num_active_nodes) + ").");

    const int N = pMesh_.N_cells_x;
    const int M = pMesh_.M_cells_y;

    // Initialise from the predictor; boundary nodes are carried over unchanged.
    u_next = u_star_;
    v_next = v_star_;

    // GradientOperator::apply() computes ∂p/∂x at interior U nodes and
    // ∂p/∂y at interior V nodes (boundary nodes left at 0).
    Vector grad_u, grad_v;
    GradientOperator(pMesh_, uMesh_, vMesh_).apply(P_sol, grad_u, grad_v);

    // ---- Correct U field: interior nodes k=2..N ----
    for (int j = 1; j <= M; ++j)
        for (int i = 2; i <= N; ++i)
            u_next[ uIdx(i, j) ] -= (dt / rho) * grad_u[ uIdx(i, j) ];

    // ---- Correct V field: interior nodes l=2..M ----
    for (int l = 2; l <= M; ++l)
        for (int i = 1; i <= N; ++i)
            v_next[ vIdx(i, l) ] -= (dt / rho) * grad_v[ vIdx(i, l) ];
}
