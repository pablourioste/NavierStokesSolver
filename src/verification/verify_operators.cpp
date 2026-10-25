// =============================================================================
// verify_operators.cpp
//
// Verification of the discrete operators with the Method of Manufactured
// Solutions (MMS).
//
// For each operator a smooth analytical function is chosen, sampled on the
// mesh, the discrete operator is applied and the result is compared with the
// analytical one. This is repeated on uniform meshes N = 8, 16, 32, 64 and the
// convergence rate is estimated (expected ≈ 2.0 for 2nd-order operators).
//
// ---------------------------------------------------------------------------
// SIGN CONVENTION OF THE LAPLACIAN MATRIX
// ---------------------------------------------------------------------------
// The FVM Laplacian assembles the positive-definite form:
//   diagonal      = +(a_E + a_W + a_N + a_S)
//   off-diagonal  = −a_nb
//
// Hence:      (A·φ)[k] = Σ a_nb*(φ_P − φ_nb) = −∫∇²φ dV ≈ −V·∇²φ
// Therefore:  (A·φ)[k] / V  ≈  −∇²φ
//
// For φ = cos(πx)cos(πy): ∇²φ = −2π²cos(πx)cos(πy)
//   ⟹  (A·φ)/V  ≈  +2π²cos(πx)cos(πy)
//
// ---------------------------------------------------------------------------
// DISCRETE IDENTITY D·G = −L/V
// ---------------------------------------------------------------------------
// GradientOperator  G : P → (U,V)  eval   ∂p/∂x, ∂p/∂y
// DivergenceOperator D : (U,V) → P  eval  ∫∇·u dV / V
// LaplaceOperator    L : P → P       eval  ∫∫∇²p dV = −A·p
//
// For interior and boundary cells:
//   D·(G·φ) = ∇²φ   and   L·φ/V = −∇²φ
//   ⟹  D·G·φ + L·φ/V = 0   (machine precision)
//
// ---------------------------------------------------------------------------
// SAMPLING U/V AT FACE POSITIONS
// ---------------------------------------------------------------------------
// The boundary U nodes (k=1, k=N+1) have their centroid at dx/4 (not at the
// wall x=0). For the DivergenceOperator test to work correctly the field must
// be sampled at the FACE position: pMesh.xvc[k-1].
// For interior U nodes (k=2..N), centroid == face (they coincide).
// =============================================================================

#include <cmath>
#include <cstdio>
#include <fstream>
#include <filesystem>
#include <iostream>
#include <iomanip>
#include <algorithm>
#include <numeric>
#include <random>
#include <vector>
#include <string>

#include "mesh/Mesh_config.h"
#include "operators/LaplaceOperator.hpp"
#include "operators/GradientOperator.hpp"
#include "operators/DivergenceOperator.hpp"
#include "operators/DiffusionOperator.hpp"

static const double PI = std::acos(-1.0);

// =============================================================================
// Uniform mesh builders — no input files needed
// =============================================================================

// Builds the uniform N×M P mesh on [0,L]×[0,H].
static GridData makeUniformPMesh(int N, int M,
                                 double L = 1.0, double H = 1.0,
                                 double W = 0.1)
{
    GridData g;
    g.N_cells_x       = N;
    g.M_cells_y       = M;
    g.num_active_nodes = N * M;

    const double dx = L / N;
    const double dy = H / M;

    // xvc[0]=0, xvc[1..N]=i*dx, xvc[N+1]=L (for sizing; the last one is the East wall)
    g.xvc.resize(N + 2);
    for (int i = 0; i <= N; ++i) g.xvc[i] = i * dx;
    g.xvc[N + 1] = L;

    g.yvc.resize(M + 2);
    for (int j = 0; j <= M; ++j) g.yvc[j] = j * dy;
    g.yvc[M + 1] = H;

    // Cell centres, with ghost cells at 0 and N+1/M+1
    g.x.resize(N + 2);
    for (int i = 1; i <= N; ++i) g.x[i] = (i - 0.5) * dx;
    g.x[0]     = -0.5 * dx;
    g.x[N + 1] =  L + 0.5 * dx;

    g.y.resize(M + 2);
    for (int j = 1; j <= M; ++j) g.y[j] = (j - 0.5) * dy;
    g.y[0]     = -0.5 * dy;
    g.y[M + 1] =  H + 0.5 * dy;

    // Face areas (2D arrays [0..N+1][0..M+1])
    g.Se.assign(N + 2, std::vector<double>(M + 2, 0.0));
    g.Sw.assign(N + 2, std::vector<double>(M + 2, 0.0));
    g.Sn.assign(N + 2, std::vector<double>(M + 2, 0.0));
    g.Ss.assign(N + 2, std::vector<double>(M + 2, 0.0));
    g.V .assign(N + 2, std::vector<double>(M + 2, 0.0));

    for (int i = 1; i <= N; ++i) {
        for (int j = 1; j <= M; ++j) {
            g.Se[i][j] = g.Sw[i][j] = dy * W;
            g.Sn[i][j] = g.Ss[i][j] = dx * W;
            g.V [i][j] = dx * dy * W;
        }
    }

    // Cells and connectivity
    g.cells.resize(N * M);
    g.id_to_i.resize(N * M);
    g.id_to_j.resize(N * M);

    for (int j = 1; j <= M; ++j) {
        for (int i = 1; i <= N; ++i) {
            const int k = (j - 1) * N + (i - 1);
            InternalCell& c = g.cells[k];
            c.global_id = k;
            c.x   = g.x[i];
            c.y   = g.y[j];
            c.Vol = g.V[i][j];

            c.bc_east  = (i == N);
            c.bc_west  = (i == 1);
            c.bc_north = (j == M);
            c.bc_south = (j == 1);

            c.n_east  = c.bc_east  ? -1 : (j - 1) * N + i;
            c.n_west  = c.bc_west  ? -1 : (j - 1) * N + (i - 2);
            c.n_north = c.bc_north ? -1 : j * N + (i - 1);
            c.n_south = c.bc_south ? -1 : (j - 2) * N + (i - 1);

            g.id_to_i[k] = i;
            g.id_to_j[k] = j;
        }
    }
    return g;
}

// Staggered (N+1)×M U mesh — same logic as MeshConfig::generate_u_mesh.
// xvc_u[0]   = pMesh.xvc[0]  (West wall)
// xvc_u[k]   = pMesh.x[k]    for k=1..N  (P cell centres)
// xvc_u[N+1] = pMesh.xvc[N]  (East wall)
static GridData makeUMesh(const GridData& p)
{
    const int N = p.N_cells_x;
    const int M = p.M_cells_y;

    GridData u;
    u.N_cells_x       = N + 1;
    u.M_cells_y       = M;
    u.num_active_nodes = (N + 1) * M;

    u.xvc.resize(N + 2);
    u.xvc[0] = p.xvc[0];
    for (int k = 1; k <= N; ++k) u.xvc[k] = p.x[k];
    u.xvc[N + 1] = p.xvc[N];

    u.yvc = p.yvc;

    u.cells.resize((N + 1) * M);
    for (int j = 1; j <= M; ++j) {
        for (int k = 1; k <= N + 1; ++k) {
            const int id = (j - 1) * (N + 1) + (k - 1);
            u.cells[id].global_id = id;
            u.cells[id].x = (u.xvc[k] + u.xvc[k - 1]) / 2.0;
            u.cells[id].y = (u.yvc[j] + u.yvc[j - 1]) / 2.0;
        }
    }
    return u;
}

// Staggered N×(M+1) V mesh.
// yvc_v[0]   = pMesh.yvc[0]  (South wall)
// yvc_v[l]   = pMesh.y[l]    for l=1..M  (P cell centres)
// yvc_v[M+1] = pMesh.yvc[M]  (North wall)
static GridData makeVMesh(const GridData& p)
{
    const int N = p.N_cells_x;
    const int M = p.M_cells_y;

    GridData v;
    v.N_cells_x       = N;
    v.M_cells_y       = M + 1;
    v.num_active_nodes = N * (M + 1);

    v.xvc = p.xvc;

    v.yvc.resize(M + 2);
    v.yvc[0] = p.yvc[0];
    for (int l = 1; l <= M; ++l) v.yvc[l] = p.y[l];
    v.yvc[M + 1] = p.yvc[M];

    v.cells.resize(N * (M + 1));
    for (int l = 1; l <= M + 1; ++l) {
        for (int i = 1; i <= N; ++i) {
            const int id = (l - 1) * N + (i - 1);
            v.cells[id].global_id = id;
            v.cells[id].x = (v.xvc[i] + v.xvc[i - 1]) / 2.0;
            v.cells[id].y = (v.yvc[l] + v.yvc[l - 1]) / 2.0;
        }
    }
    return v;
}

// =============================================================================
// L∞ error metric over two vectors of the same size
// =============================================================================
static double errorLinf(const Vector& computed, const Vector& exact)
{
    double err = 0.0;
    for (std::size_t i = 0; i < computed.size(); ++i)
        err = std::max(err, std::abs(computed[i] - exact[i]));
    return err;
}

// =============================================================================
// Suite 1 — LaplaceOperator
//
// φ(x,y) = cos(πx)cos(πy)
//   ∂φ/∂n = 0 on every wall → homogeneous Neumann satisfied ✓
//   ∇²φ = −2π²cos(πx)cos(πy)
//   (A·φ)/V ≈ −∇²φ = +2π²cos(πx)cos(πy)   (positive-definite convention)
// =============================================================================
static double runLaplace(int N)
{
    const GridData pMesh = makeUniformPMesh(N, N);
    LaplaceOperator L(pMesh);

    Vector phi(pMesh.num_active_nodes);
    for (int k = 0; k < pMesh.num_active_nodes; ++k) {
        phi[k] = std::cos(PI * pMesh.cells[k].x)
               * std::cos(PI * pMesh.cells[k].y);
    }

    const Vector Aphi = L.apply(phi);

    // (A·φ)/V  vs  +2π²cos(πx)cos(πy)
    Vector disc(pMesh.num_active_nodes), exact(pMesh.num_active_nodes);
    for (int k = 0; k < pMesh.num_active_nodes; ++k) {
        const int i = pMesh.id_to_i[k];
        const int j = pMesh.id_to_j[k];
        disc[k]  = Aphi[k] / pMesh.V[i][j];
        exact[k] = +2.0 * PI * PI
                 * std::cos(PI * pMesh.cells[k].x)
                 * std::cos(PI * pMesh.cells[k].y);
    }
    return errorLinf(disc, exact);
}

// =============================================================================
// Suite 2 — GradientOperator
//
// φ(x,y) = cos(πx)cos(πy)
//   ∂φ/∂x = −πsin(πx)cos(πy)  at interior U nodes (face between P cells)
//   ∂φ/∂y = −πcos(πx)sin(πy)  at interior V nodes
//
// Interior U nodes have centroid == P face, so the comparison is direct.
// =============================================================================
static std::pair<double,double> runGradient(int N)
{
    const GridData pMesh = makeUniformPMesh(N, N);
    const GridData uMesh = makeUMesh(pMesh);
    const GridData vMesh = makeVMesh(pMesh);
    GradientOperator G(pMesh, uMesh, vMesh);

    // Sample φ at P centres
    Vector phi(pMesh.num_active_nodes);
    for (int k = 0; k < pMesh.num_active_nodes; ++k) {
        phi[k] = std::cos(PI * pMesh.cells[k].x)
               * std::cos(PI * pMesh.cells[k].y);
    }

    Vector gu, gv;
    G.apply(phi, gu, gv);

    // Extract the interior U nodes (k=2..N) and compare
    Vector gu_int, gu_ex, gv_int, gv_ex;
    for (int j = 1; j <= N; ++j) {
        for (int k = 2; k <= N; ++k) {
            const int id = (j - 1) * (N + 1) + (k - 1);
            const double xu = uMesh.cells[id].x;  // centroid = face
            const double yu = uMesh.cells[id].y;
            gu_int.push_back(gu[id]);
            gu_ex .push_back(-PI * std::sin(PI * xu) * std::cos(PI * yu));
        }
    }
    for (int l = 2; l <= N; ++l) {
        for (int i = 1; i <= N; ++i) {
            const int id = (l - 1) * N + (i - 1);
            const double xv = vMesh.cells[id].x;
            const double yv = vMesh.cells[id].y;  // centroid = face
            gv_int.push_back(gv[id]);
            gv_ex .push_back(-PI * std::cos(PI * xv) * std::sin(PI * yv));
        }
    }
    return { errorLinf(gu_int, gu_ex), errorLinf(gv_int, gv_ex) };
}

// =============================================================================
// Suite 3 — DivergenceOperator
//
// u(x,y) = sin(2πx)cos(πy)   →  ∂u/∂x = 2πcos(2πx)cos(πy)
// v(x,y) = cos(πx)sin(2πy)   →  ∂v/∂y = 2πcos(πx)cos(2πy)
// ∇·(u,v) = 2πcos(2πx)cos(πy) + 2πcos(πx)cos(2πy)
//
// IMPORTANT: u is sampled at the FACE POSITION x = pMesh.xvc[k-1] (not at the
// centroid of the boundary U node, which sits at dx/4). For interior nodes
// centroid == face, so there is no difference. For k=1 (West wall): x=0.
// Likewise v is sampled at y = pMesh.yvc[l-1].
// =============================================================================
static double runDivergence(int N)
{
    const GridData pMesh = makeUniformPMesh(N, N);
    const GridData uMesh = makeUMesh(pMesh);
    const GridData vMesh = makeVMesh(pMesh);
    DivergenceOperator D(pMesh, uMesh, vMesh);

    // Sample u at the FACE position: x_face = pMesh.xvc[k-1]
    Vector u(uMesh.num_active_nodes), v(vMesh.num_active_nodes);
    for (int j = 1; j <= N; ++j) {
        for (int k = 1; k <= N + 1; ++k) {
            const int id  = (j - 1) * (N + 1) + (k - 1);
            const double xf = pMesh.xvc[k - 1];   // face = pMesh.xvc[k-1]
            const double yc = pMesh.y[j];           // P-centre y = U-centre y
            u[id] = std::sin(2.0 * PI * xf) * std::cos(PI * yc);
        }
    }
    for (int l = 1; l <= N + 1; ++l) {
        for (int i = 1; i <= N; ++i) {
            const int id  = (l - 1) * N + (i - 1);
            const double xc = pMesh.x[i];            // P-centre x = V-centre x
            const double yf = pMesh.yvc[l - 1];      // face = pMesh.yvc[l-1]
            v[id] = std::cos(PI * xc) * std::sin(2.0 * PI * yf);
        }
    }

    const Vector div_h = D.apply(u, v);

    // Compare with ∇·(u,v) at P cell centres
    Vector exact(pMesh.num_active_nodes);
    for (int k = 0; k < pMesh.num_active_nodes; ++k) {
        const double xc = pMesh.cells[k].x;
        const double yc = pMesh.cells[k].y;
        exact[k] = 2.0 * PI * std::cos(2.0 * PI * xc) * std::cos(PI * yc)
                 + 2.0 * PI * std::cos(PI * xc)        * std::cos(2.0 * PI * yc);
    }
    return errorLinf(div_h, exact);
}

// =============================================================================
// Suite 4 — DiffusionOperator
//
// u(x,y) = sin(πy)   (independent of x  →  ∂²u/∂x²=0 exactly)
//   ν∇²u = ν·d²u/dy² = −νπ²sin(πy),   ν=1
//
// Since the field is constant in x, the non-uniform x spacing (caused by the
// boundary node k=1 sitting at dx/4 instead of x=0) does not affect the result.
// The y term uses the image method (ghost cell) at j=1 and j=M, which is
// 2nd order.
//
// Note: sin(πy)=0 at y=0 and y=1 → u_ghost = 2·0 − u_P = −u_P, which is the
// exact field value in the mirror cell (sin(π·(−dy/2)) = −sin(π·dy/2) = −u_P). ✓
// =============================================================================
static double runDiffusion(int N)
{
    const GridData pMesh = makeUniformPMesh(N, N);
    const GridData uMesh = makeUMesh(pMesh);
    const GridData vMesh = makeVMesh(pMesh);

    const double nu = 1.0;
    // sin(πy)=0 on the South and North walls → ghost=0
    DiffusionOperator Diff(pMesh, uMesh, vMesh, nu,
                           /*u_ghost_N=*/0.0, /*u_ghost_S=*/0.0,
                           /*v_ghost_W=*/0.0, /*v_ghost_E=*/0.0);

    // Sample u = sin(πy) at U node centres
    Vector u(uMesh.num_active_nodes);
    for (int id = 0; id < uMesh.num_active_nodes; ++id) {
        u[id] = std::sin(PI * uMesh.cells[id].y);
    }

    Vector diff_u;
    Diff.applyU(u, diff_u);

    // Compare at interior U nodes k=2..N, j=1..M
    // (j=1 and j=M use the 2nd-order ghost cell ✓)
    Vector computed, exact;
    for (int j = 1; j <= N; ++j) {
        for (int k = 2; k <= N; ++k) {
            const int id = (j - 1) * (N + 1) + (k - 1);
            computed.push_back(diff_u[id]);
            const double y = uMesh.cells[id].y;
            exact.push_back(-nu * PI * PI * std::sin(PI * y));
        }
    }
    return errorLinf(computed, exact);
}

// =============================================================================
// Suite 5 — Discrete identity D·G + L/V = 0
//
// For any field φ and for ALL cells (boundary cells included):
//   D·(G·φ)  +  (L·φ)/V  =  0   to machine precision
//
// Proof for the interior East face of cell (i,j):
//   Contribution of D·G:  Se/V * (φ[i+1,j]−φ[i,j]) / (x[i+1]−x[i]) = +a_E/V*(φ_E−φ_P)
//   Contribution of L/V:  [−a_E*(φ_E−φ_P) + ...] / V = −a_E/V*(φ_E−φ_P)
//   Sum = 0 ✓
//
// For a boundary face (bc=true): G returns 0 and L omits that term → both = 0 ✓
// =============================================================================
static double runDGL(int N)
{
    const GridData pMesh = makeUniformPMesh(N, N);
    const GridData uMesh = makeUMesh(pMesh);
    const GridData vMesh = makeVMesh(pMesh);

    LaplaceOperator    L(pMesh);
    GradientOperator   G(pMesh, uMesh, vMesh);
    DivergenceOperator D(pMesh, uMesh, vMesh);

    // Field φ = cos(πx)cos(πy) at P centres
    Vector phi(pMesh.num_active_nodes);
    for (int k = 0; k < pMesh.num_active_nodes; ++k) {
        phi[k] = std::cos(PI * pMesh.cells[k].x)
               * std::cos(PI * pMesh.cells[k].y);
    }

    // D·(G·φ)
    Vector gu, gv;
    G.apply(phi, gu, gv);
    const Vector DGphi = D.apply(gu, gv);

    // (L·φ) / V
    const Vector Lphi = L.apply(phi);
    Vector LphiV(pMesh.num_active_nodes);
    for (int k = 0; k < pMesh.num_active_nodes; ++k) {
        const int i = pMesh.id_to_i[k];
        const int j = pMesh.id_to_j[k];
        LphiV[k] = Lphi[k] / pMesh.V[i][j];
    }

    // D·G·φ + L·φ/V must be 0 (not D·G·φ - L·φ/V)
    Vector sum(pMesh.num_active_nodes);
    for (int k = 0; k < pMesh.num_active_nodes; ++k)
        sum[k] = DGphi[k] + LphiV[k];

    Vector zero(pMesh.num_active_nodes, 0.0);
    return errorLinf(sum, zero);
}

// =============================================================================
// CSV helper — appends one row to output/verification/data/convergence.csv
// =============================================================================
static void csvRow(std::ofstream& f,
                   const std::string& suite, int N, double error,
                   const std::string& error2_label = "", double error2 = -1.0)
{
    f << suite << "," << N << "," << std::scientific << std::setprecision(6) << error;
    if (error2 >= 0.0)
        f << "," << error2_label << "," << error2;
    f << "\n";
}

// =============================================================================
// Suite 6 — DiffusionOperator applyV
//
// v(x,y) = sin(πx)   (independent of y  →  ∂²v/∂y²=0 exactly)
//   ν∇²v = ν·d²v/dx² = −νπ²sin(πx),   ν=1
//
// The West (x=0) and East (x=1) walls have v=0: sin(π·0)=sin(π·1)=0,
// so v_ghost_W=0 and v_ghost_E=0 are consistent with the function.
// The boundary V nodes (l=1, l=M+1) are set to 0 by applyV and are not
// included in the comparison.
// =============================================================================
static double runDiffusionV(int N)
{
    const GridData pMesh = makeUniformPMesh(N, N);
    const GridData uMesh = makeUMesh(pMesh);
    const GridData vMesh = makeVMesh(pMesh);

    const double nu = 1.0;
    DiffusionOperator Diff(pMesh, uMesh, vMesh, nu,
                           /*u_ghost_N=*/0.0, /*u_ghost_S=*/0.0,
                           /*v_ghost_W=*/0.0, /*v_ghost_E=*/0.0);

    // Sample v = sin(πx) at V node centres
    Vector v(vMesh.num_active_nodes);
    for (int id = 0; id < vMesh.num_active_nodes; ++id)
        v[id] = std::sin(PI * vMesh.cells[id].x);

    Vector diff_v;
    Diff.applyV(v, diff_v);

    // Compare at interior V nodes: l=2..M (all i=1..N)
    // (i=1 and i=N use the ghost cell — they check the W/E logic)
    Vector computed, exact;
    for (int l = 2; l <= N; ++l) {
        for (int i = 1; i <= N; ++i) {
            const int id = (l - 1) * N + (i - 1);
            computed.push_back(diff_v[id]);
            const double x = vMesh.cells[id].x;
            exact.push_back(-nu * PI * PI * std::sin(PI * x));
        }
    }
    return errorLinf(computed, exact);
}

// =============================================================================
// Suite 7 — DiffusionOperator with a non-zero boundary condition (moving lid)
//
// u(x,y) = y   (linear in y  →  ∂²u/∂y²=0  →  ν∇²u=0 exactly)
//
// South wall: u_wall=0  →  ghost = 2·0 − u_P = −u_P
// North wall: u_wall=1  →  ghost = 2·1 − u_P
//
// For a linear field the second-order difference cancels exactly
// (analytical proof in docs/OPERATORS_REPORT.md, Section 4, Suite 7).
// The expected error is < 1e-12 for every N (algebraic cancellation, not O(h²)).
// =============================================================================
static double runDiffusionLinear(int N)
{
    const GridData pMesh = makeUniformPMesh(N, N);
    const GridData uMesh = makeUMesh(pMesh);
    const GridData vMesh = makeVMesh(pMesh);

    const double nu  = 1.0;
    const double U_S = 0.0;  // no-slip south
    const double U_N = 1.0;  // lid north

    DiffusionOperator Diff(pMesh, uMesh, vMesh, nu,
                           /*u_ghost_N=*/U_N, /*u_ghost_S=*/U_S,
                           /*v_ghost_W=*/0.0, /*v_ghost_E=*/0.0);

    // Sample u = y at U node centres
    Vector u(uMesh.num_active_nodes);
    for (int id = 0; id < uMesh.num_active_nodes; ++id)
        u[id] = uMesh.cells[id].y;

    Vector diff_u;
    Diff.applyU(u, diff_u);

    // Compare at interior U nodes k=2..N, all j
    // The exact result is 0 for a linear field
    Vector computed, exact_zero;
    for (int j = 1; j <= N; ++j) {
        for (int k = 2; k <= N; ++k) {
            const int id = (j - 1) * (N + 1) + (k - 1);
            computed.push_back(diff_u[id]);
            exact_zero.push_back(0.0);
        }
    }
    return errorLinf(computed, exact_zero);
}

// =============================================================================
// Suite 8 — Discrete adjoint: Gᵀ = −D  (machine precision)
//
// Symmetry-preserving property:
//   φᵀ·(D·u)  +  (G·φ)ᵀ·u  =  0   ∀ φ∈P, u∈U (interior nodes)
//
// This is the discrete form of the integration-by-parts identity:
//   ∫(∇·u)φ dV + ∫u·∇φ dV = ∮φu·n dS = 0   (with u·n=0 on the walls)
//
// Repeated for the v component.  Expected error: ~1e-13.
// =============================================================================
static double runAdjoint(int N)
{
    const GridData pMesh = makeUniformPMesh(N, N);
    const GridData uMesh = makeUMesh(pMesh);
    const GridData vMesh = makeVMesh(pMesh);

    GradientOperator   G(pMesh, uMesh, vMesh);
    DivergenceOperator D(pMesh, uMesh, vMesh);

    // Reproducible random vectors (fixed seed)
    std::mt19937 rng(42);
    std::uniform_real_distribution<double> dist(-1.0, 1.0);

    // Field φ on P (any values)
    Vector phi(pMesh.num_active_nodes);
    for (auto& v : phi) v = dist(rng);

    // Field u on U — boundary nodes (k=1, k=N+1) set to zero so that the
    // boundary term of the integration-by-parts identity vanishes
    Vector u_vec(uMesh.num_active_nodes, 0.0);
    const int Np1 = N + 1;
    for (int j = 1; j <= N; ++j)
        for (int k = 2; k <= N; ++k)
            u_vec[(j - 1) * Np1 + (k - 1)] = dist(rng);

    // Field v on V — boundary nodes (l=1, l=M+1) set to zero
    Vector v_vec(vMesh.num_active_nodes, 0.0);
    for (int l = 2; l <= N; ++l)
        for (int i = 1; i <= N; ++i)
            v_vec[(l - 1) * N + (i - 1)] = dist(rng);

    // φᵀ·(D·u) — D.apply() returns flux/V, so the product without volume
    // weighting is equivalent to the product in the pure P space.
    // The identity D^T = -G is checked in the unweighted norm:
    //   φᵀ·(D·u) + (G·φ)ᵀ·u = 0
    const Vector Du = D.apply(u_vec, v_vec);
    double phiDu = 0.0;
    for (int k = 0; k < pMesh.num_active_nodes; ++k) {
        phiDu += phi[k] * Du[k];  // unweighted: D already divides by V internally
    }

    // (G·φ)ᵀ·u — gradient of φ on interior U/V nodes
    Vector gu, gv;
    G.apply(phi, gu, gv);
    double Gphi_u = 0.0;
    for (int j = 1; j <= N; ++j)
        for (int k = 2; k <= N; ++k)
            Gphi_u += gu[(j - 1) * Np1 + (k - 1)] * u_vec[(j - 1) * Np1 + (k - 1)];
    for (int l = 2; l <= N; ++l)
        for (int i = 1; i <= N; ++i)
            Gphi_u += gv[(l - 1) * N + (i - 1)] * v_vec[(l - 1) * N + (i - 1)];

    return std::abs(phiDu + Gphi_u);
}

// =============================================================================
// Suite 9 — LaplaceOperator: symmetry and positive definiteness
//
// Test A (symmetry):  |φᵀ·(A·ψ) − ψᵀ·(A·φ)| < ε_machine
// Test B (pos. def.): φᵀ·A·φ > 0  for a random non-zero φ
//
// Returns the symmetry error (Test A) and the quadratic form (Test B).
// =============================================================================
static std::pair<double, double> runSPD(int N)
{
    const GridData pMesh = makeUniformPMesh(N, N);
    LaplaceOperator L(pMesh);

    std::mt19937 rng(123);
    std::uniform_real_distribution<double> dist(-1.0, 1.0);

    Vector phi(pMesh.num_active_nodes), psi(pMesh.num_active_nodes);
    for (auto& v : phi) v = dist(rng);
    for (auto& v : psi) v = dist(rng);

    const Vector Aphi = L.apply(phi);
    const Vector Apsi = L.apply(psi);

    // Test A: φᵀ·Aψ vs ψᵀ·Aφ
    double phiApsi = 0.0, psiAphi = 0.0;
    for (int k = 0; k < pMesh.num_active_nodes; ++k) {
        phiApsi += phi[k] * Apsi[k];
        psiAphi += psi[k] * Aphi[k];
    }
    const double sym_err = std::abs(phiApsi - psiAphi);

    // Test B: φᵀ·Aφ  (pin a reference value to remove the constant null mode)
    // Set phi[0]=0 so that the constant does not contribute
    Vector phi2 = phi;
    phi2[0] = 0.0;
    const Vector Aphi2 = L.apply(phi2);
    double quad = 0.0;
    for (int k = 0; k < pMesh.num_active_nodes; ++k)
        quad += phi2[k] * Aphi2[k];

    return { sym_err, quad };
}

// =============================================================================
// Suite 10 — LaplaceOperator on a stretched mesh (tanh, γ=2)
//
// Same manufactured function as Suite 1: φ=cos(πx)cos(πy)
// Mesh generated with the hyperbolic-tangent mapping:
//   x_map = 1/2 · (1 + tanh(γ·(2ξ−1)) / tanh(γ)),  γ=2
//
// NOTE on convergence:
// The FVM operator uses a_E = Se / (x[i+1] − x[i]) where x[i] are cell centres.
// On non-uniform meshes (h_e ≠ h_w), the cross term φ'·(h_e − h_w)/dx_i
// introduces a LOCAL error of order O(1/N) (first order), since:
//   h_e − h_w = f''(ξ)/N² + O(N⁻³)   and   dx_i = f'(ξ)/N + O(N⁻²)
//   → error = φ'·f''(ξ)/(f'(ξ)·N) = O(1/N)
//
// With the correct factor 2/(h_e + h_w) the scheme would recover 2nd order.
// This test exposes the limitation of the current stencil on stretched meshes.
// Expected rate: ~1.0
// =============================================================================
static GridData makeStretchedPMesh(int N, int M,
                                   double L = 1.0, double H = 1.0,
                                   double W = 0.1, double gamma = 2.0)
{
    GridData g;
    g.N_cells_x        = N;
    g.M_cells_y        = M;
    g.num_active_nodes = N * M;

    const double tg = std::tanh(gamma);

    // Face positions in x (0..N) with tanh stretching
    auto mapX = [&](double xi) {
        return 0.5 * L * (1.0 + std::tanh(gamma * (2.0 * xi - 1.0)) / tg);
    };
    auto mapY = [&](double eta) {
        return 0.5 * H * (1.0 + std::tanh(gamma * (2.0 * eta - 1.0)) / tg);
    };

    g.xvc.resize(N + 2);
    for (int i = 0; i <= N; ++i) g.xvc[i] = mapX(static_cast<double>(i) / N);
    g.xvc[N + 1] = L;

    g.yvc.resize(M + 2);
    for (int j = 0; j <= M; ++j) g.yvc[j] = mapY(static_cast<double>(j) / M);
    g.yvc[M + 1] = H;

    // Cell centres (average of the faces)
    g.x.resize(N + 2);
    for (int i = 1; i <= N; ++i) g.x[i] = 0.5 * (g.xvc[i - 1] + g.xvc[i]);
    g.x[0]     = 2.0 * g.xvc[0] - g.x[1];
    g.x[N + 1] = 2.0 * g.xvc[N] - g.x[N];

    g.y.resize(M + 2);
    for (int j = 1; j <= M; ++j) g.y[j] = 0.5 * (g.yvc[j - 1] + g.yvc[j]);
    g.y[0]     = 2.0 * g.yvc[0] - g.y[1];
    g.y[M + 1] = 2.0 * g.yvc[M] - g.y[M];

    // Areas and volumes
    g.Se.assign(N + 2, std::vector<double>(M + 2, 0.0));
    g.Sw.assign(N + 2, std::vector<double>(M + 2, 0.0));
    g.Sn.assign(N + 2, std::vector<double>(M + 2, 0.0));
    g.Ss.assign(N + 2, std::vector<double>(M + 2, 0.0));
    g.V .assign(N + 2, std::vector<double>(M + 2, 0.0));

    for (int i = 1; i <= N; ++i) {
        const double dx = g.xvc[i] - g.xvc[i - 1];
        for (int j = 1; j <= M; ++j) {
            const double dy = g.yvc[j] - g.yvc[j - 1];
            g.Se[i][j] = g.Sw[i][j] = dy * W;
            g.Sn[i][j] = g.Ss[i][j] = dx * W;
            g.V [i][j] = dx * dy * W;
        }
    }

    // Cells and connectivity
    g.cells.resize(N * M);
    g.id_to_i.resize(N * M);
    g.id_to_j.resize(N * M);
    for (int j = 1; j <= M; ++j) {
        for (int i = 1; i <= N; ++i) {
            const int k = (j - 1) * N + (i - 1);
            InternalCell& c = g.cells[k];
            c.global_id = k;
            c.x   = g.x[i];
            c.y   = g.y[j];
            c.Vol = g.V[i][j];
            c.bc_east  = (i == N);
            c.bc_west  = (i == 1);
            c.bc_north = (j == M);
            c.bc_south = (j == 1);
            c.n_east  = c.bc_east  ? -1 : (j - 1) * N + i;
            c.n_west  = c.bc_west  ? -1 : (j - 1) * N + (i - 2);
            c.n_north = c.bc_north ? -1 : j * N + (i - 1);
            c.n_south = c.bc_south ? -1 : (j - 2) * N + (i - 1);
            g.id_to_i[k] = i;
            g.id_to_j[k] = j;
        }
    }
    return g;
}

static double runLaplaceStretched(int N)
{
    const GridData pMesh = makeStretchedPMesh(N, N);
    LaplaceOperator L(pMesh);

    Vector phi(pMesh.num_active_nodes);
    for (int k = 0; k < pMesh.num_active_nodes; ++k) {
        phi[k] = std::cos(PI * pMesh.cells[k].x)
               * std::cos(PI * pMesh.cells[k].y);
    }

    const Vector Aphi = L.apply(phi);

    Vector disc(pMesh.num_active_nodes), exact(pMesh.num_active_nodes);
    for (int k = 0; k < pMesh.num_active_nodes; ++k) {
        const int i = pMesh.id_to_i[k];
        const int j = pMesh.id_to_j[k];
        disc[k]  = Aphi[k] / pMesh.V[i][j];
        exact[k] = +2.0 * PI * PI
                 * std::cos(PI * pMesh.cells[k].x)
                 * std::cos(PI * pMesh.cells[k].y);
    }
    return errorLinf(disc, exact);
}

// =============================================================================
// main — runs the 10 suites and prints the convergence tables
// =============================================================================
int main()
{
    const std::vector<int> sizes = { 8, 16, 32, 64 };
    const int W = 14;

    // ---- Open the output CSV ----
    // Create the output directory if it does not exist
    std::system("mkdir -p output");
    std::filesystem::create_directories("output/verification/data");
    std::ofstream csv("output/verification/data/convergence.csv");
    csv << "suite,N,error_linf\n";

    auto rate = [](double e_coarse, double e_fine) -> std::string {
        if (e_coarse < 1e-300 || e_fine < 1e-300) return "  —  ";
        double r = std::log2(e_coarse / e_fine);
        char buf[16];
        std::snprintf(buf, sizeof(buf), "%.3f", r);
        return std::string(buf);
    };

    std::cout << std::fixed << std::setprecision(4);

    // ---- Suite 1: Laplacian ----
    std::cout << "\n=== Suite 1: LaplaceOperator  φ=cos(πx)cos(πy)  (A·φ)/V ≈ +2π²cos(πx)cos(πy) ===\n";
    std::cout << std::setw(6)  << "N"
              << std::setw(W)  << "Error L∞"
              << std::setw(W)  << "Rate\n";
    double prev1 = 0.0;
    for (int N : sizes) {
        const double e = runLaplace(N);
        std::cout << std::setw(6) << N
                  << std::setw(W) << std::scientific << e
                  << std::setw(W) << rate(prev1, e) << "\n";
        csvRow(csv, "Laplace_uniform", N, e);
        prev1 = e;
    }

    // ---- Suite 2: Gradient ----
    std::cout << "\n=== Suite 2: GradientOperator  φ=cos(πx)cos(πy) ===\n";
    std::cout << std::setw(6)    << "N"
              << std::setw(W+2)  << "Error ∂p/∂x"
              << std::setw(W)    << "Rate"
              << std::setw(W+2)  << "Error ∂p/∂y"
              << std::setw(W)    << "Rate\n";
    double prev_eu = 0.0, prev_ev = 0.0;
    for (int N : sizes) {
        auto [eu, ev] = runGradient(N);
        std::cout << std::setw(6)    << N
                  << std::setw(W+2)  << std::scientific << eu
                  << std::setw(W)    << rate(prev_eu, eu)
                  << std::setw(W+2)  << std::scientific << ev
                  << std::setw(W)    << rate(prev_ev, ev) << "\n";
        csvRow(csv, "Gradient_dpx", N, eu);
        csvRow(csv, "Gradient_dpy", N, ev);
        prev_eu = eu;
        prev_ev = ev;
    }

    // ---- Suite 3: Divergence ----
    std::cout << "\n=== Suite 3: DivergenceOperator  u=sin(2πx)cos(πy), v=cos(πx)sin(2πy) ===\n";
    std::cout << "    (sampled at face positions: x=pMesh.xvc[k-1], y=pMesh.yvc[l-1])\n";
    std::cout << std::setw(6) << "N"
              << std::setw(W) << "Error L∞"
              << std::setw(W) << "Rate\n";
    double prev3 = 0.0;
    for (int N : sizes) {
        const double e = runDivergence(N);
        std::cout << std::setw(6) << N
                  << std::setw(W) << std::scientific << e
                  << std::setw(W) << rate(prev3, e) << "\n";
        csvRow(csv, "Divergence", N, e);
        prev3 = e;
    }

    // ---- Suite 4: Diffusion U ----
    std::cout << "\n=== Suite 4: DiffusionOperator  u=sin(πy) [const. in x], ν=1 ===\n";
    std::cout << "    (independence of x removes the 1st-order error at nodes k=2, k=N)\n";
    std::cout << std::setw(6) << "N"
              << std::setw(W) << "Error L∞"
              << std::setw(W) << "Rate\n";
    double prev4 = 0.0;
    for (int N : sizes) {
        const double e = runDiffusion(N);
        std::cout << std::setw(6) << N
                  << std::setw(W) << std::scientific << e
                  << std::setw(W) << rate(prev4, e) << "\n";
        csvRow(csv, "Diffusion_U", N, e);
        prev4 = e;
    }

    // ---- Suite 5: Identity D·G + L/V = 0 ----
    std::cout << "\n=== Suite 5: Identity D·G·φ + (L·φ)/V = 0  (machine precision) ===\n";
    std::cout << std::setw(6) << "N"
              << std::setw(W+2) << "max|D·G·φ + L·φ/V|\n";
    for (int N : sizes) {
        const double e = runDGL(N);
        std::cout << std::setw(6) << N
                  << std::setw(W+2) << std::scientific << e << "\n";
        csvRow(csv, "DGL_identity", N, e);
    }

    std::cout << "\nExpected rates: ~2.0 for Suites 1-4, ~1e-13 for Suite 5.\n";

    // =========================================================================
    // ADDITIONAL SUITES (6-10)
    // =========================================================================

    // ---- Suite 6: Diffusion V ----
    std::cout << "\n=== Suite 6: DiffusionOperator applyV  v=sin(πx) [const. in y], ν=1 ===\n";
    std::cout << "    (checks the West/East ghost cells, independent v-component path)\n";
    std::cout << std::setw(6) << "N"
              << std::setw(W) << "Error L∞"
              << std::setw(W) << "Rate\n";
    double prev6 = 0.0;
    for (int N : sizes) {
        const double e = runDiffusionV(N);
        std::cout << std::setw(6) << N
                  << std::setw(W) << std::scientific << e
                  << std::setw(W) << rate(prev6, e) << "\n";
        csvRow(csv, "Diffusion_V", N, e);
        prev6 = e;
    }

    // ---- Suite 7: Diffusion with a non-zero BC (moving lid) ----
    std::cout << "\n=== Suite 7: DiffusionOperator u=y (linear), u_N=1, u_S=0  →  ν∇²u = 0 exactly ===\n";
    std::cout << "    (algebraic cancellation for a linear field — expected error < 1e-12)\n";
    std::cout << std::setw(6) << "N"
              << std::setw(W+2) << "max|applyU(y)|\n";
    for (int N : sizes) {
        const double e = runDiffusionLinear(N);
        std::cout << std::setw(6) << N
                  << std::setw(W+2) << std::scientific << e << "\n";
        csvRow(csv, "Diffusion_linear_BC", N, e);
    }

    // ---- Suite 8: Adjoint Gᵀ = −D ----
    std::cout << "\n=== Suite 8: Discrete adjoint  φᵀ·(D·u) + (G·φ)ᵀ·u = 0  (machine precision) ===\n";
    std::cout << std::setw(6) << "N"
              << std::setw(W+2) << "Adjoint residual\n";
    for (int N : sizes) {
        const double e = runAdjoint(N);
        std::cout << std::setw(6) << N
                  << std::setw(W+2) << std::scientific << e << "\n";
        csvRow(csv, "Adjoint_GD", N, e);
    }

    // ---- Suite 9: LaplaceOperator SPD ----
    std::cout << "\n=== Suite 9: LaplaceOperator SPD — symmetry and positive definiteness ===\n";
    std::cout << std::setw(6) << "N"
              << std::setw(W+2) << "Symmetry err"
              << std::setw(W+2) << "φᵀAφ (>0?)\n";
    bool all_pd = true;
    for (int N : sizes) {
        auto [sym, quad] = runSPD(N);
        std::cout << std::setw(6) << N
                  << std::setw(W+2) << std::scientific << sym
                  << std::setw(W+2) << std::scientific << quad;
        if (quad <= 0.0) { std::cout << "  *** PD FAILURE ***"; all_pd = false; }
        std::cout << "\n";
        csvRow(csv, "SPD_symmetry", N, sym);
    }
    if (all_pd)
        std::cout << "    Positive definiteness: OK for every N.\n";

    // ---- Suite 10: Laplacian on a stretched mesh ----
    std::cout << "\n=== Suite 10: LaplaceOperator on a tanh mesh (γ=2)  φ=cos(πx)cos(πy) ===\n";
    std::cout << "    (stencil a_E=Se/(x[i+1]-x[i]) is 1st order on a non-uniform mesh)\n";
    std::cout << "    (expected rate ~1.0; for 2nd order use the factor 2/(h_e+h_w))\n";
    std::cout << std::setw(6)  << "N"
              << std::setw(W)  << "Error L∞"
              << std::setw(W)  << "Rate\n";
    double prev10 = 0.0;
    for (int N : sizes) {
        const double e = runLaplaceStretched(N);
        std::cout << std::setw(6) << N
                  << std::setw(W) << std::scientific << e
                  << std::setw(W) << rate(prev10, e) << "\n";
        csvRow(csv, "Laplace_stretched", N, e);
        prev10 = e;
    }

    std::cout << "\nExpected rates: ~2.0 for Suites 1-4, 6; ~1.0 for Suite 10 (non-uniform mesh).\n";
    std::cout << "Expected errors < 1e-12: Suites 5, 7, 8, 9 (algebraic identities).\n";
    std::cout << "\nConvergence data written to: output/verification/data/convergence.csv\n";

    return 0;
}
