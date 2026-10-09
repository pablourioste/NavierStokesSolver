#include "EnergySolver.hpp"
#include <stdexcept>

// ---------------------------------------------------------------------------
// extractWallT — read a side's temperature BC (type + wall value).
// Uniform along the side (uses the first face), mirroring MomentumSolver's
// per-wall ghost extraction. Defaults to adiabatic if the list is empty.
// ---------------------------------------------------------------------------
void EnergySolver::extractWallT(const std::vector<BoundaryFace>& bound,
                                BCType& type, double& value)
{
    if (bound.empty()) { type = WALL_ADIABATIC; value = 0.0; return; }
    type  = bound[0].type;
    value = bound[0].value;
}

EnergySolver::EnergySolver(const GridData& pMesh,
                           const GridData& uMesh,
                           const GridData& vMesh,
                           double alpha,
                           ConvectionScheme scheme)
    : pMesh_(pMesh), uMesh_(uMesh), vMesh_(vMesh),
      alpha_(alpha), scheme_(scheme)
{
    if (pMesh_.cells.empty())
        throw std::runtime_error(
            "EnergySolver: pMesh cells are empty. "
            "Call build_connectivity() before constructing EnergySolver.");

    extractWallT(pMesh_.bound_west_T,  bc_west_,  Tw_west_);
    extractWallT(pMesh_.bound_east_T,  bc_east_,  Tw_east_);
    extractWallT(pMesh_.bound_north_T, bc_north_, Tw_north_);
    extractWallT(pMesh_.bound_south_T, bc_south_, Tw_south_);
}

// ---------------------------------------------------------------------------
// tRHS — alpha*laplacian(theta) - (u.grad)theta at interior P cell (i,j)
//
// Diffusion: 2nd-order FD on the non-uniform P-grid with the ghost-cell image
//   method at walls (same scheme as DiffusionOperator):
//     isothermal: theta_ghost = 2*theta_wall - theta_P
//     adiabatic : theta_ghost = theta_P  (zero normal gradient)
//
// Convection: conservative FVM divergence of (u*theta) per unit volume. The
//   wall-normal velocity is zero (no penetration), so boundary-face fluxes
//   vanish and no temperature ghost is needed for convection.
// ---------------------------------------------------------------------------
double EnergySolver::tRHS(const Vector& T, const Vector& u, const Vector& v,
                          int i, int j) const
{
    const int N = pMesh_.N_cells_x;
    const int M = pMesh_.M_cells_y;

    const double T_P = T[pIdx(i, j)];
    const double x_P = pMesh_.cells[pIdx(i, j)].x;
    const double y_P = pMesh_.cells[pIdx(i, j)].y;

    // ---- ghost helper for a Dirichlet/Neumann wall ----
    auto wallGhost = [&](BCType type, double Twall) -> double {
        return (type == WALL_ISOTHERMAL) ? (2.0 * Twall - T_P) : T_P;
    };

    // ---- x-neighbours / distances ----
    double T_E, dx_e, T_W, dx_w;
    if (i < N) { T_E = T[pIdx(i+1, j)]; dx_e = pMesh_.cells[pIdx(i+1, j)].x - x_P; }
    else       { T_E = wallGhost(bc_east_, Tw_east_); dx_e = 2.0 * (pMesh_.xvc[N] - x_P); }
    if (i > 1) { T_W = T[pIdx(i-1, j)]; dx_w = x_P - pMesh_.cells[pIdx(i-1, j)].x; }
    else       { T_W = wallGhost(bc_west_, Tw_west_); dx_w = 2.0 * (x_P - pMesh_.xvc[0]); }

    // ---- y-neighbours / distances ----
    double T_N, dy_n, T_S, dy_s;
    if (j < M) { T_N = T[pIdx(i, j+1)]; dy_n = pMesh_.cells[pIdx(i, j+1)].y - y_P; }
    else       { T_N = wallGhost(bc_north_, Tw_north_); dy_n = 2.0 * (pMesh_.yvc[M] - y_P); }
    if (j > 1) { T_S = T[pIdx(i, j-1)]; dy_s = y_P - pMesh_.cells[pIdx(i, j-1)].y; }
    else       { T_S = wallGhost(bc_south_, Tw_south_); dy_s = 2.0 * (y_P - pMesh_.yvc[0]); }

    const double d2T_dx2 = 2.0 / (dx_e + dx_w) * ((T_E - T_P) / dx_e - (T_P - T_W) / dx_w);
    const double d2T_dy2 = 2.0 / (dy_n + dy_s) * ((T_N - T_P) / dy_n - (T_P - T_S) / dy_s);
    const double diff = alpha_ * (d2T_dx2 + d2T_dy2);

    // ---- Convection: conservative flux through the four faces ----
    // Face-normal velocities are the staggered velocities directly.
    const double u_e = u[uIdx(i+1, j)];  // east  face
    const double u_w = u[uIdx(i,   j)];  // west  face
    const double v_n = v[vIdx(i, j+1)];  // north face
    const double v_s = v[vIdx(i, j  )];  // south face

    const double Se = pMesh_.Se[i][j];
    const double Sw = pMesh_.Sw[i][j];
    const double Sn = pMesh_.Sn[i][j];
    const double Ss = pMesh_.Ss[i][j];
    const double Vol = pMesh_.V[i][j];

    double Fe = 0.0, Fw = 0.0, Fn = 0.0, Fs = 0.0;
    if (scheme_ == ConvectionScheme::CDS) {
        Fe = (i < N) ? u_e * 0.5 * (T_P + T[pIdx(i+1, j)]) * Se : 0.0;
        Fw = (i > 1) ? u_w * 0.5 * (T[pIdx(i-1, j)] + T_P) * Sw : 0.0;
        Fn = (j < M) ? v_n * 0.5 * (T_P + T[pIdx(i, j+1)]) * Sn : 0.0;
        Fs = (j > 1) ? v_s * 0.5 * (T[pIdx(i, j-1)] + T_P) * Ss : 0.0;
    } else {
        // First-order upwind: take the value from the upstream cell.
        Fe = (i < N) ? u_e * ((u_e >= 0.0) ? T_P : T[pIdx(i+1, j)]) * Se : 0.0;
        Fw = (i > 1) ? u_w * ((u_w >= 0.0) ? T[pIdx(i-1, j)] : T_P) * Sw : 0.0;
        Fn = (j < M) ? v_n * ((v_n >= 0.0) ? T_P : T[pIdx(i, j+1)]) * Sn : 0.0;
        Fs = (j > 1) ? v_s * ((v_s >= 0.0) ? T[pIdx(i, j-1)] : T_P) * Ss : 0.0;
    }

    const double conv = (Fe - Fw + Fn - Fs) / Vol;

    return diff - conv;
}

// ---------------------------------------------------------------------------
// advance — Adams-Bashforth 2 (Euler on first call). Forward-Euler-style
// buffering: the whole RHS is computed from T^n before any update so node
// updates do not pollute neighbours within a step.
// ---------------------------------------------------------------------------
void EnergySolver::advance(Vector& T, const Vector& u, const Vector& v, double dt)
{
    const int N = pMesh_.N_cells_x;
    const int M = pMesh_.M_cells_y;

    Vector rhs(T.size(), 0.0);
    for (int j = 1; j <= M; ++j)
        for (int i = 1; i <= N; ++i)
            rhs[pIdx(i, j)] = tRHS(T, u, v, i, j);

    if (first_step_) {
        for (std::size_t k = 0; k < T.size(); ++k) T[k] += dt * rhs[k];
        first_step_ = false;
    } else {
        for (std::size_t k = 0; k < T.size(); ++k)
            T[k] += dt * (1.5 * rhs[k] - 0.5 * rhs_old_[k]);
    }
    rhs_old_ = std::move(rhs);
}
