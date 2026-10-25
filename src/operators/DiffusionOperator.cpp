#include "DiffusionOperator.hpp"

DiffusionOperator::DiffusionOperator(const GridData& pMesh,
                                     const GridData& uMesh,
                                     const GridData& vMesh,
                                     double nu,
                                     double u_ghost_north,
                                     double u_ghost_south,
                                     double v_ghost_west,
                                     double v_ghost_east)
    : pMesh_(pMesh), uMesh_(uMesh), vMesh_(vMesh),
      nu_(nu),
      u_ghost_north_(u_ghost_north), u_ghost_south_(u_ghost_south),
      v_ghost_west_(v_ghost_west),   v_ghost_east_(v_ghost_east)
{}

// ---------------------------------------------------------------------------
// diffuseU — ν∇²u at the interior U node (k,j), k=2..N, j=1..M
//
// Code extracted from MomentumSolver::uRHS() — diffusive part only.
// ---------------------------------------------------------------------------
double DiffusionOperator::diffuseU(const Vector& u, int k, int j) const
{
    const int M = pMesh_.M_cells_y;

    const double x_P = uMesh_.cells[uIdx(k, j)].x;
    const double y_P = uMesh_.cells[uIdx(k, j)].y;
    const double u_P0 = u[uIdx(k, j)];

    // Neighbours in x (k=2..N guarantees k-1>=1 and k+1<=N+1 inside the global
    // domain, but the neighbouring U node may fall inside an immersed obstacle).
    const bool west_blocked_x = uNodeBlocked(pMesh_, k - 1, j);
    const bool east_blocked_x = uNodeBlocked(pMesh_, k + 1, j);

    double u_W, dx_w;
    if (!west_blocked_x) {
        u_W  = u[uIdx(k - 1, j)];
        dx_w = x_P - uMesh_.cells[uIdx(k - 1, j)].x;
    } else {
        u_W  = -u_P0;
        dx_w = 2.0 * (x_P - pMesh_.obstacle_xmax);   // solid to the west: wall = obstacle EAST face
    }

    double u_E, dx_e;
    if (!east_blocked_x) {
        u_E  = u[uIdx(k + 1, j)];
        dx_e = uMesh_.cells[uIdx(k + 1, j)].x - x_P;
    } else {
        u_E  = -u_P0;
        dx_e = 2.0 * (pMesh_.obstacle_xmin - x_P);   // solid to the east: wall = obstacle WEST face
    }

    // Neighbours in y — ghost cells at the N/S walls, or at the surface of an
    // immersed obstacle if the neighbouring U node (k,j-1)/(k,j+1) falls inside
    // the solid block (uNodeBlocked). If node (k,j) itself lies ON the surface
    // (uNodeOnObstacleSurface), it was already set to 0 by
    // applyBoundaryConditions before this call — no special case needed.
    const bool south_blocked = (j > 1) && uNodeBlocked(pMesh_, k, j - 1);
    const bool north_blocked = (j < M) && uNodeBlocked(pMesh_, k, j + 1);

    double u_S, dy_s;
    if (j > 1 && !south_blocked) {
        u_S  = u[uIdx(k, j - 1)];
        dy_s = y_P - uMesh_.cells[uIdx(k, j - 1)].y;
    } else if (south_blocked) {
        // No-slip on the obstacle surface: ghost = -u_P (u_wall = 0),
        // distance to the obstacle's top face (exact, snapped coordinate).
        u_S  = -u_P0;
        dy_s = 2.0 * (y_P - pMesh_.obstacle_ymax);   // solid to the south: wall = obstacle TOP face
    } else {
        // South wall: image method → ghost = 2·u_wall − u_P
        u_S  = 2.0 * u_ghost_south_ - u_P0;
        dy_s = 2.0 * (y_P - uMesh_.yvc[0]);
    }

    double u_N, dy_n;
    if (j < M && !north_blocked) {
        u_N  = u[uIdx(k, j + 1)];
        dy_n = uMesh_.cells[uIdx(k, j + 1)].y - y_P;
    } else if (north_blocked) {
        u_N  = -u_P0;
        dy_n = 2.0 * (pMesh_.obstacle_ymin - y_P);   // solid to the north: wall = obstacle BOTTOM face
    } else {
        // North wall: image method → ghost = 2·u_wall − u_P
        u_N  = 2.0 * u_ghost_north_ - u_P0;
        dy_n = 2.0 * (uMesh_.yvc[M] - y_P);
    }

    const double u_P = u_P0;

    const double d2u_dx2 = 2.0 / (dx_e + dx_w) * ((u_E - u_P) / dx_e - (u_P - u_W) / dx_w);
    const double d2u_dy2 = 2.0 / (dy_n + dy_s) * ((u_N - u_P) / dy_n - (u_P - u_S) / dy_s);

    return nu_ * (d2u_dx2 + d2u_dy2);
}

// ---------------------------------------------------------------------------
// diffuseV — ν∇²v at the interior V node (i,l), i=1..N, l=2..M
//
// Code extracted from MomentumSolver::vRHS() — diffusive part only.
// ---------------------------------------------------------------------------
double DiffusionOperator::diffuseV(const Vector& v, int i, int l) const
{
    const int N = pMesh_.N_cells_x;

    const double x_P = vMesh_.cells[vIdx(i, l)].x;
    const double y_P = vMesh_.cells[vIdx(i, l)].y;
    const double v_P0 = v[vIdx(i, l)];

    // Neighbours in y (l=2..M guarantees l-1>=1 and l+1<=M+1 inside the global
    // domain, but the neighbouring V node may fall inside an immersed obstacle).
    const bool south_blocked_y = vNodeBlocked(pMesh_, i, l - 1);
    const bool north_blocked_y = vNodeBlocked(pMesh_, i, l + 1);

    double v_S, dy_s;
    if (!south_blocked_y) {
        v_S  = v[vIdx(i, l - 1)];
        dy_s = y_P - vMesh_.cells[vIdx(i, l - 1)].y;
    } else {
        v_S  = -v_P0;
        dy_s = 2.0 * (y_P - pMesh_.obstacle_ymax);   // solid to the south: wall = obstacle TOP face
    }

    double v_N, dy_n;
    if (!north_blocked_y) {
        v_N  = v[vIdx(i, l + 1)];
        dy_n = vMesh_.cells[vIdx(i, l + 1)].y - y_P;
    } else {
        v_N  = -v_P0;
        dy_n = 2.0 * (pMesh_.obstacle_ymin - y_P);   // solid to the north: wall = obstacle BOTTOM face
    }

    // Neighbours in x — ghost cells at the W/E walls, or at the surface of an
    // immersed obstacle if the neighbouring V node falls inside the solid block.
    const bool west_blocked = (i > 1) && vNodeBlocked(pMesh_, i - 1, l);
    const bool east_blocked = (i < N) && vNodeBlocked(pMesh_, i + 1, l);

    double v_W, dx_w;
    if (i > 1 && !west_blocked) {
        v_W  = v[vIdx(i - 1, l)];
        dx_w = x_P - vMesh_.cells[vIdx(i - 1, l)].x;
    } else if (west_blocked) {
        v_W  = -v_P0;
        dx_w = 2.0 * (x_P - pMesh_.obstacle_xmax);   // solid to the west: wall = obstacle EAST face
    } else {
        // West wall: image method → ghost = 2·v_wall − v_P
        v_W  = 2.0 * v_ghost_west_ - v_P0;
        dx_w = 2.0 * (x_P - vMesh_.xvc[0]);
    }

    double v_E, dx_e;
    if (i < N && !east_blocked) {
        v_E  = v[vIdx(i + 1, l)];
        dx_e = vMesh_.cells[vIdx(i + 1, l)].x - x_P;
    } else if (east_blocked) {
        v_E  = -v_P0;
        dx_e = 2.0 * (pMesh_.obstacle_xmin - x_P);   // solid to the east: wall = obstacle WEST face
    } else {
        // East wall: image method → ghost = 2·v_wall − v_P
        v_E  = 2.0 * v_ghost_east_ - v_P0;
        dx_e = 2.0 * (vMesh_.xvc[N] - x_P);
    }

    const double v_P = v_P0;

    const double d2v_dx2 = 2.0 / (dx_e + dx_w) * ((v_E - v_P) / dx_e - (v_P - v_W) / dx_w);
    const double d2v_dy2 = 2.0 / (dy_n + dy_s) * ((v_N - v_P) / dy_n - (v_P - v_S) / dy_s);

    return nu_ * (d2v_dx2 + d2v_dy2);
}

void DiffusionOperator::applyU(const Vector& u, Vector& result) const
{
    const int N = pMesh_.N_cells_x;
    const int M = pMesh_.M_cells_y;
    result.assign(uMesh_.num_active_nodes, 0.0);

    for (int j = 1; j <= M; ++j)
        for (int k = 2; k <= N; ++k)
            result[uIdx(k, j)] = diffuseU(u, k, j);
}

void DiffusionOperator::applyV(const Vector& v, Vector& result) const
{
    const int N = pMesh_.N_cells_x;
    const int M = pMesh_.M_cells_y;
    result.assign(vMesh_.num_active_nodes, 0.0);

    for (int l = 2; l <= M; ++l)
        for (int i = 1; i <= N; ++i)
            result[vIdx(i, l)] = diffuseV(v, i, l);
}
