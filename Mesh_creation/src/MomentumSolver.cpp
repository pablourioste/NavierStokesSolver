#include "MomentumSolver.hpp"
#include <stdexcept>

// =============================================================================
//  STAGGERED GRID LAYOUT RECAP
//
//  U-cell (k,j) occupies x ∈ [uMesh.xvc[k-1], uMesh.xvc[k]]
//                         y ∈ [uMesh.yvc[j-1], uMesh.yvc[j]]
//  U-node (k,j) sits at cell centroid: cells[uIdx(k,j)].x / .y
//
//  V-cell (i,l) occupies x ∈ [vMesh.xvc[i-1], vMesh.xvc[i]]
//                         y ∈ [vMesh.yvc[l-1], vMesh.yvc[l]]
//  V-node (i,l) sits at cell centroid: cells[vIdx(i,l)].x / .y
//
//  Wall positions (used as ghost references):
//    West  x = uMesh.xvc[0]  = vMesh.xvc[0]  = pMesh.xvc[0]
//    East  x = uMesh.xvc[N+1]= vMesh.xvc[N]  = pMesh.xvc[N]
//    South y = uMesh.yvc[0]  = vMesh.yvc[0]  = pMesh.yvc[0]
//    North y = uMesh.yvc[M]  = vMesh.yvc[M+1]= pMesh.yvc[M]
//
//  Ghost values for wall-adjacent stencils are derived from pMesh.bound_*
//  at construction time via ghostValue(BCType, faceVal).
// =============================================================================

// ---------------------------------------------------------------------------
// ghostValue — maps a BCType + stored face value to a velocity ghost value
// ---------------------------------------------------------------------------
double MomentumSolver::ghostValue(BCType type, double faceVal)
{
    switch (type) {
        case WALL_FIXED_VALUE: return faceVal;
        case INLET:      return faceVal;
        case WALL_ADIABATIC:
        case OUTLET:
        default:         return 0.0;
    }
}

// ---------------------------------------------------------------------------
// extractGhost — returns the ghost value from the first face of a boundary
//               list, or 0.0 if the list is empty (uniform wall assumed).
// ---------------------------------------------------------------------------
double MomentumSolver::extractGhost(const std::vector<BoundaryFace>& bound)
{
    return bound.empty() ? 0.0 : ghostValue(bound[0].type, bound[0].value);
}

// ---------------------------------------------------------------------------
// Constructor — derives ghost values from pMesh.bound_* bocos
// ---------------------------------------------------------------------------
MomentumSolver::MomentumSolver(const GridData& pMesh,
                                const GridData& uMesh,
                                const GridData& vMesh,
                                double nu,
                                ConvectionScheme scheme)
    : pMesh_(pMesh), uMesh_(uMesh), vMesh_(vMesh),
      nu_(nu), scheme_(scheme),
      u_ghost_north_(extractGhost(pMesh.bound_north)),
      u_ghost_south_(extractGhost(pMesh.bound_south)),
      v_ghost_west_ (extractGhost(pMesh.bound_west)),
      v_ghost_east_ (extractGhost(pMesh.bound_east)),
      diffOp_(pMesh, uMesh, vMesh, nu,
              extractGhost(pMesh.bound_north),
              extractGhost(pMesh.bound_south),
              extractGhost(pMesh.bound_west),
              extractGhost(pMesh.bound_east))
{
    if (uMesh_.cells.empty() || vMesh_.cells.empty())
        throw std::runtime_error(
            "MomentumSolver: uMesh/vMesh cells are empty. "
            "Call build_connectivity() before constructing MomentumSolver.");
}

// ---------------------------------------------------------------------------
// applyBoundaryConditions — sets boundary nodes in u/v from pMesh.bound_* bocos
// ---------------------------------------------------------------------------
void MomentumSolver::applyBoundaryConditions(Vector& u, Vector& v) const
{
    const int N = pMesh_.N_cells_x;
    const int M = pMesh_.M_cells_y;

    // -- U field boundary nodes (normal component at WEST/EAST faces) --
    // u-nodes at k=1 (west) and k=N+1 (east) sit ON the wall.
    // INLET  → prescribed normal velocity = face value
    // WALL*  → no-slip, u = 0
    // OUTLET → zero-gradient approx, u = 0
    // The tangential ghost (u above/below NORTH/SOUTH wall) lives in u_ghost_north_/south_
    // and is used inside uRHS — NOT set here.
    for (int j = 1; j <= M; ++j) {
        const BCType wt = pMesh_.bound_west.empty() ? WALL_ADIABATIC
                                                     : pMesh_.bound_west[j-1].type;
        const BCType et = pMesh_.bound_east.empty() ? WALL_ADIABATIC
                                                     : pMesh_.bound_east[j-1].type;
        u[uIdx(1,   j)] = (wt == INLET) ? pMesh_.bound_west[j-1].value : 0.0;
        u[uIdx(N+1, j)] = (et == INLET) ? pMesh_.bound_east[j-1].value : 0.0;
    }

    // -- V field boundary nodes (normal component at SOUTH/NORTH faces) --
    // v-nodes at l=1 (south) and l=M+1 (north) sit ON the wall.
    // Normal v must be 0 for any wall BC (no penetration).
    // WALL_FIXED_VALUE on NORTH encodes the tangential u of the lid — that goes
    // into u_ghost_north_ (set in the constructor), NOT into v here.
    // The tangential ghost (v left/right of WEST/EAST wall) lives in v_ghost_west_/east_
    // and is used inside vRHS — NOT set here.
    for (int i = 1; i <= N; ++i) {
        v[vIdx(i, 1  )] = 0.0;   // south wall: no penetration
        v[vIdx(i, M+1)] = 0.0;   // north wall: no penetration
    }
}


// ---------------------------------------------------------------------------
// uRHS — u-momentum right-hand side at interior U-node (k,j)
//
// Returns:  ν*(d²u/dx² + d²u/dy²)  −  (u*du/dx + v_interp*du/dy)
//
// Interior range: k=2..N, j=1..M
// ---------------------------------------------------------------------------
double MomentumSolver::uRHS(const Vector& u, const Vector& v,
                             int k, int j) const
{
    const int M = pMesh_.M_cells_y;

    const double x_P = uMesh_.cells[uIdx(k, j)].x;
    const double y_P = uMesh_.cells[uIdx(k, j)].y;

    // ---- Neighbours and distances in x ----
    // k=2..N guarantees k-1>=1 and k+1<=N+1, both in the vector.
    const double u_E  = u[uIdx(k+1, j)];
    const double u_W  = u[uIdx(k-1, j)];
    const double dx_e = uMesh_.cells[uIdx(k+1, j)].x - x_P;
    const double dx_w = x_P - uMesh_.cells[uIdx(k-1, j)].x;

    // ---- Neighbours and distances in y ----
    // South: j=1 hits the no-slip bottom wall (ghost, not in vector).
    double u_S, dy_s;
    if (j > 1) {
        u_S  = u[uIdx(k, j-1)];
        dy_s = y_P - uMesh_.cells[uIdx(k, j-1)].y;
    } else {
        // Image/ghost-cell method: ghost = 2*u_wall - u_P, distance = 2*(y_P - wall)
        // Ensures (u_ghost + u_P)/2 = u_wall at the south face, and gives
        // d²u/dy² = (2*u_wall - 3*u_P + u_N)/dy²  (matches reference)
        u_S  = 2.0 * u_ghost_south_ - u[uIdx(k, j)];
        dy_s = 2.0 * (y_P - uMesh_.yvc[0]);
    }

    // North: j=M hits the north wall ghost (not in vector).
    double u_N, dy_n;
    if (j < M) {
        u_N  = u[uIdx(k, j+1)];
        dy_n = uMesh_.cells[uIdx(k, j+1)].y - y_P;
    } else {
        // Image/ghost-cell method: ghost = 2*u_wall - u_P, distance = 2*(wall - y_P)
        // For lid: u_wall = U_ref → ghost = 2*U_ref - u_P
        u_N  = 2.0 * u_ghost_north_ - u[uIdx(k, j)];
        dy_n = 2.0 * (uMesh_.yvc[M] - y_P);
    }

    const double u_P = u[uIdx(k, j)];

    // ---- Diffusion: delegated to DiffusionOperator ----
    const double diff = diffOp_.diffuseU(u, k, j);

    // ---- v interpolated to U-node (k,j) ----
    // Four surrounding V-nodes: (k-1,j), (k,j), (k-1,j+1), (k,j+1)
    // k=2..N  →  i=k-1 ∈ [1,N] and i=k ∈ [2,N] (N is also valid), both in range.
    // l=j and l=j+1: j ∈ [1,M] → l ∈ [1,M] and l+1 ∈ [2,M+1], all valid.
    const double v_interp = 0.25 * (v[vIdx(k-1, j  )] + v[vIdx(k, j  )]
                                  + v[vIdx(k-1, j+1)] + v[vIdx(k, j+1)]);

    // ---- Convection ----
    double conv;
    if (scheme_ == ConvectionScheme::CDS) {
        conv = u_P * (u_E - u_W) / (dx_e + dx_w)
             + v_interp * (u_N - u_S) / (dy_n + dy_s);
    } else {
        // First-order upwind
        const double du_dx = (u_P >= 0.0) ? (u_P - u_W)/dx_w : (u_E - u_P)/dx_e;
        const double du_dy = (v_interp >= 0.0) ? (u_P - u_S)/dy_s : (u_N - u_P)/dy_n;
        conv = u_P * du_dx + v_interp * du_dy;
    }

    return diff - conv;
}


// ---------------------------------------------------------------------------
// vRHS — v-momentum right-hand side at interior V-node (i,l)
//
// Returns:  ν*(d²v/dx² + d²v/dy²)  −  (u_interp*dv/dx + v*dv/dy)
//
// Interior range: i=1..N, l=2..M
// ---------------------------------------------------------------------------
double MomentumSolver::vRHS(const Vector& u, const Vector& v,
                             int i, int l) const
{
    const int N = pMesh_.N_cells_x;

    const double x_P = vMesh_.cells[vIdx(i, l)].x;
    const double y_P = vMesh_.cells[vIdx(i, l)].y;

    // ---- Neighbours and distances in y ----
    // l=2..M guarantees l-1>=1 and l+1<=M+1, both in the vector.
    const double v_N  = v[vIdx(i, l+1)];
    const double v_S  = v[vIdx(i, l-1)];
    const double dy_n = vMesh_.cells[vIdx(i, l+1)].y - y_P;
    const double dy_s = y_P - vMesh_.cells[vIdx(i, l-1)].y;

    // ---- Neighbours and distances in x ----
    // West: i=1 hits the no-slip west wall (ghost, not in vector).
    double v_W, dx_w;
    if (i > 1) {
        v_W  = v[vIdx(i-1, l)];
        dx_w = x_P - vMesh_.cells[vIdx(i-1, l)].x;
    } else {
        // Image/ghost-cell method: ghost = 2*v_wall - v_P, distance = 2*(x_P - wall)
        v_W  = 2.0 * v_ghost_west_ - v[vIdx(i, l)];
        dx_w = 2.0 * (x_P - vMesh_.xvc[0]);
    }

    // East: i=N hits the east wall ghost (not in vector).
    double v_E, dx_e;
    if (i < N) {
        v_E  = v[vIdx(i+1, l)];
        dx_e = vMesh_.cells[vIdx(i+1, l)].x - x_P;
    } else {
        // Image/ghost-cell method: ghost = 2*v_wall - v_P, distance = 2*(wall - x_P)
        v_E  = 2.0 * v_ghost_east_ - v[vIdx(i, l)];
        dx_e = 2.0 * (vMesh_.xvc[N] - x_P);
    }

    const double v_P = v[vIdx(i, l)];

    // ---- Diffusion: delegated to DiffusionOperator ----
    const double diff = diffOp_.diffuseV(v, i, l);

    // ---- u interpolated to V-node (i,l) ----
    // Four surrounding U-nodes: (i,l-1), (i+1,l-1), (i,l), (i+1,l)
    // l=2..M  →  j=l-1 ∈ [1,M-1] and j=l ∈ [2,M], both valid U-rows.
    // i+1 for i=N: uIdx(N+1,j) is the east boundary node (set to 0 by applyBC) ✓
    const double u_interp = 0.25 * (u[uIdx(i,   l-1)] + u[uIdx(i+1, l-1)]
                                  + u[uIdx(i,   l  )] + u[uIdx(i+1, l  )]);

    // ---- Convection ----
    double conv;
    if (scheme_ == ConvectionScheme::CDS) {
        conv = u_interp * (v_E - v_W) / (dx_e + dx_w)
             + v_P * (v_N - v_S) / (dy_n + dy_s);
    } else {
        // First-order upwind
        const double dv_dx = (u_interp >= 0.0) ? (v_P - v_W)/dx_w : (v_E - v_P)/dx_e;
        const double dv_dy = (v_P >= 0.0) ? (v_P - v_S)/dy_s : (v_N - v_P)/dy_n;
        conv = u_interp * dv_dx + v_P * dv_dy;
    }

    return diff - conv;
}


// ---------------------------------------------------------------------------
// advance — Adams-Bashforth 2 predictor step
//
// AB2:  u* = u^n + dt * [ (3/2) R(u^n)  −  (1/2) R(u^{n-1}) ]
//
// The RHS  R = ν∇²u − (u·∇)u  is computed at every interior node into a
// temporary buffer (forward-Euler style, no in-place update) so that u/v
// changes in one node do not pollute neighbours within the same step.
//
// First call: R(u^{n-1}) is not yet available → falls back to Euler.
// Subsequent calls: uses the stored previous-step RHS for AB2.
// ---------------------------------------------------------------------------
void MomentumSolver::advance(Vector& u, Vector& v, double dt)
{
    const int N = pMesh_.N_cells_x;
    const int M = pMesh_.M_cells_y;

    // --- Compute current RHS R^n ---
    Vector rhs_u(u.size(), 0.0);
    for (int j = 1; j <= M; ++j)
        for (int k = 2; k <= N; ++k)
            rhs_u[uIdx(k, j)] = uRHS(u, v, k, j);

    Vector rhs_v(v.size(), 0.0);
    for (int l = 2; l <= M; ++l)
        for (int i = 1; i <= N; ++i)
            rhs_v[vIdx(i, l)] = vRHS(u, v, i, l);

    if (first_step_) {
        // Step 0: no R^{n-1} available → plain Euler
        for (std::size_t idx = 0; idx < u.size(); ++idx) u[idx] += dt * rhs_u[idx];
        for (std::size_t idx = 0; idx < v.size(); ++idx) v[idx] += dt * rhs_v[idx];
        first_step_ = false;
    } else {
        // AB2: u* = u^n + dt * (3/2 * R^n  -  1/2 * R^{n-1})
        for (std::size_t idx = 0; idx < u.size(); ++idx)
            u[idx] += dt * (1.5 * rhs_u[idx] - 0.5 * rhs_u_old_[idx]);
        for (std::size_t idx = 0; idx < v.size(); ++idx)
            v[idx] += dt * (1.5 * rhs_v[idx] - 0.5 * rhs_v_old_[idx]);
    }

    // Store R^n for next step
    rhs_u_old_ = std::move(rhs_u);
    rhs_v_old_ = std::move(rhs_v);
}
