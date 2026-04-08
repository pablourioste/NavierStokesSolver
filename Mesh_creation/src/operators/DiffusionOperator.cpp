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
// diffuseU — ν∇²u en el nodo U interior (k,j), k=2..N, j=1..M
//
// Código extraído de MomentumSolver::uRHS() — solo la parte difusiva.
// ---------------------------------------------------------------------------
double DiffusionOperator::diffuseU(const Vector& u, int k, int j) const
{
    const int M = pMesh_.M_cells_y;

    const double x_P = uMesh_.cells[uIdx(k, j)].x;
    const double y_P = uMesh_.cells[uIdx(k, j)].y;

    // Vecinos en x (k=2..N garantiza k-1>=1 y k+1<=N+1)
    const double u_E  = u[uIdx(k + 1, j)];
    const double u_W  = u[uIdx(k - 1, j)];
    const double dx_e = uMesh_.cells[uIdx(k + 1, j)].x - x_P;
    const double dx_w = x_P - uMesh_.cells[uIdx(k - 1, j)].x;

    // Vecinos en y — celdas fantasma en las paredes N/S
    double u_S, dy_s;
    if (j > 1) {
        u_S  = u[uIdx(k, j - 1)];
        dy_s = y_P - uMesh_.cells[uIdx(k, j - 1)].y;
    } else {
        // Pared Sur: método imagen → ghost = 2·u_wall − u_P
        u_S  = 2.0 * u_ghost_south_ - u[uIdx(k, j)];
        dy_s = 2.0 * (y_P - uMesh_.yvc[0]);
    }

    double u_N, dy_n;
    if (j < M) {
        u_N  = u[uIdx(k, j + 1)];
        dy_n = uMesh_.cells[uIdx(k, j + 1)].y - y_P;
    } else {
        // Pared Norte: método imagen → ghost = 2·u_wall − u_P
        u_N  = 2.0 * u_ghost_north_ - u[uIdx(k, j)];
        dy_n = 2.0 * (uMesh_.yvc[M] - y_P);
    }

    const double u_P = u[uIdx(k, j)];

    const double d2u_dx2 = 2.0 / (dx_e + dx_w) * ((u_E - u_P) / dx_e - (u_P - u_W) / dx_w);
    const double d2u_dy2 = 2.0 / (dy_n + dy_s) * ((u_N - u_P) / dy_n - (u_P - u_S) / dy_s);

    return nu_ * (d2u_dx2 + d2u_dy2);
}

// ---------------------------------------------------------------------------
// diffuseV — ν∇²v en el nodo V interior (i,l), i=1..N, l=2..M
//
// Código extraído de MomentumSolver::vRHS() — solo la parte difusiva.
// ---------------------------------------------------------------------------
double DiffusionOperator::diffuseV(const Vector& v, int i, int l) const
{
    const int N = pMesh_.N_cells_x;

    const double x_P = vMesh_.cells[vIdx(i, l)].x;
    const double y_P = vMesh_.cells[vIdx(i, l)].y;

    // Vecinos en y (l=2..M garantiza l-1>=1 y l+1<=M+1)
    const double v_N  = v[vIdx(i, l + 1)];
    const double v_S  = v[vIdx(i, l - 1)];
    const double dy_n = vMesh_.cells[vIdx(i, l + 1)].y - y_P;
    const double dy_s = y_P - vMesh_.cells[vIdx(i, l - 1)].y;

    // Vecinos en x — celdas fantasma en las paredes O/E
    double v_W, dx_w;
    if (i > 1) {
        v_W  = v[vIdx(i - 1, l)];
        dx_w = x_P - vMesh_.cells[vIdx(i - 1, l)].x;
    } else {
        // Pared Oeste: método imagen → ghost = 2·v_wall − v_P
        v_W  = 2.0 * v_ghost_west_ - v[vIdx(i, l)];
        dx_w = 2.0 * (x_P - vMesh_.xvc[0]);
    }

    double v_E, dx_e;
    if (i < N) {
        v_E  = v[vIdx(i + 1, l)];
        dx_e = vMesh_.cells[vIdx(i + 1, l)].x - x_P;
    } else {
        // Pared Este: método imagen → ghost = 2·v_wall − v_P
        v_E  = 2.0 * v_ghost_east_ - v[vIdx(i, l)];
        dx_e = 2.0 * (vMesh_.xvc[N] - x_P);
    }

    const double v_P = v[vIdx(i, l)];

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
