// =============================================================================
// verify_operators.cpp
//
// Verificación de los operadores discretos mediante el Método de Soluciones
// Manufacturadas (MMS).
//
// Para cada operador se elige una función analítica suave, se muestrea en la
// malla, se aplica el operador discreto y se compara con el resultado analítico.
// Se repite en mallas uniformes N = 8, 16, 32, 64 y se estima la tasa de
// convergencia (esperada ≈ 2.0 para operadores de 2.º orden).
//
// ---------------------------------------------------------------------------
// CONVENCIÓN DE SIGNO DE LA MATRIZ LAPLACIANA
// ---------------------------------------------------------------------------
// El Laplaciano FVM ensambla la forma positivo-definida:
//   diagonal    = +(a_E + a_W + a_N + a_S)
//   fuera-diag  = −a_nb
//
// Por ello:   (A·φ)[k] = Σ a_nb*(φ_P − φ_nb) = −∫∇²φ dV ≈ −V·∇²φ
// Luego:      (A·φ)[k] / V  ≈  −∇²φ
//
// Para φ = cos(πx)cos(πy): ∇²φ = −2π²cos(πx)cos(πy)
//   ⟹  (A·φ)/V  ≈  +2π²cos(πx)cos(πy)
//
// ---------------------------------------------------------------------------
// IDENTIDAD DISCRETA D·G = −L/V
// ---------------------------------------------------------------------------
// GradientOperator  G : P → (U,V)  eval   ∂p/∂x, ∂p/∂y
// DivergenceOperator D : (U,V) → P  eval  ∫∇·u dV / V
// LaplaceOperator    L : P → P       eval  ∫∫∇²p dV = −A·p
//
// Para celdas interiores y de frontera:
//   D·(G·φ) = ∇²φ   y   L·φ/V = −∇²φ
//   ⟹  D·G·φ + L·φ/V = 0   (precisión máquina)
//
// ---------------------------------------------------------------------------
// MUESTREO EN POSICIONES DE CARA PARA U/V
// ---------------------------------------------------------------------------
// Los nodos U de frontera (k=1, k=N+1) tienen su centroide en dx/4 (no en
// la pared x=0). Para que la DivergenceOperator funcione correctamente el
// campo debe ser muestreado en la posición de CARA: pMesh.xvc[k-1].
// Para los nodos U interiores (k=2..N), centroide == cara (coinciden).
// =============================================================================

#include <cmath>
#include <cstdio>
#include <fstream>
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
// Constructores de malla uniforme — sin archivos de entrada
// =============================================================================

// Construye la malla P uniforme N×M en [0,L]×[0,H].
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

    // xvc[0]=0, xvc[1..N]=i*dx, xvc[N+1]=L (por tamaño; último es la pared Este)
    g.xvc.resize(N + 2);
    for (int i = 0; i <= N; ++i) g.xvc[i] = i * dx;
    g.xvc[N + 1] = L;

    g.yvc.resize(M + 2);
    for (int j = 0; j <= M; ++j) g.yvc[j] = j * dy;
    g.yvc[M + 1] = H;

    // Centros de celda con celdas fantasma en 0 y N+1/M+1
    g.x.resize(N + 2);
    for (int i = 1; i <= N; ++i) g.x[i] = (i - 0.5) * dx;
    g.x[0]     = -0.5 * dx;
    g.x[N + 1] =  L + 0.5 * dx;

    g.y.resize(M + 2);
    for (int j = 1; j <= M; ++j) g.y[j] = (j - 0.5) * dy;
    g.y[0]     = -0.5 * dy;
    g.y[M + 1] =  H + 0.5 * dy;

    // Áreas de cara (arreglos 2D [0..N+1][0..M+1])
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

    // Celdas y conectividad
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

// Malla U escalonada (N+1)×M — misma lógica que MeshConfig::generate_u_mesh.
// xvc_u[0]   = pMesh.xvc[0]  (pared Oeste)
// xvc_u[k]   = pMesh.x[k]    para k=1..N  (centros de celda P)
// xvc_u[N+1] = pMesh.xvc[N]  (pared Este)
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

// Malla V escalonada N×(M+1).
// yvc_v[0]   = pMesh.yvc[0]  (pared Sur)
// yvc_v[l]   = pMesh.y[l]    para l=1..M  (centros de celda P)
// yvc_v[M+1] = pMesh.yvc[M]  (pared Norte)
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
// Métrica de error L∞ sobre vectores del mismo tamaño
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
//   ∂φ/∂n = 0 en todas las paredes → Neumann homogéneo satisfecho ✓
//   ∇²φ = −2π²cos(πx)cos(πy)
//   (A·φ)/V ≈ −∇²φ = +2π²cos(πx)cos(πy)   (convención positivo-definida)
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
//   ∂φ/∂x = −πsin(πx)cos(πy)  en nodos U interiores (cara entre P-celdas)
//   ∂φ/∂y = −πcos(πx)sin(πy)  en nodos V interiores
//
// Los nodos U interiores tienen centroide == cara de P, así que la comparación
// es directa.
// =============================================================================
static std::pair<double,double> runGradient(int N)
{
    const GridData pMesh = makeUniformPMesh(N, N);
    const GridData uMesh = makeUMesh(pMesh);
    const GridData vMesh = makeVMesh(pMesh);
    GradientOperator G(pMesh, uMesh, vMesh);

    // Muestrear φ en centros P
    Vector phi(pMesh.num_active_nodes);
    for (int k = 0; k < pMesh.num_active_nodes; ++k) {
        phi[k] = std::cos(PI * pMesh.cells[k].x)
               * std::cos(PI * pMesh.cells[k].y);
    }

    Vector gu, gv;
    G.apply(phi, gu, gv);

    // Extraer nodos U interiores (k=2..N) y comparar
    Vector gu_int, gu_ex, gv_int, gv_ex;
    for (int j = 1; j <= N; ++j) {
        for (int k = 2; k <= N; ++k) {
            const int id = (j - 1) * (N + 1) + (k - 1);
            const double xu = uMesh.cells[id].x;  // centroide = cara
            const double yu = uMesh.cells[id].y;
            gu_int.push_back(gu[id]);
            gu_ex .push_back(-PI * std::sin(PI * xu) * std::cos(PI * yu));
        }
    }
    for (int l = 2; l <= N; ++l) {
        for (int i = 1; i <= N; ++i) {
            const int id = (l - 1) * N + (i - 1);
            const double xv = vMesh.cells[id].x;
            const double yv = vMesh.cells[id].y;  // centroide = cara
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
// IMPORTANTE: se muestrea u en la POSICIÓN DE CARA x = pMesh.xvc[k-1] (no en
// el centroide del nodo U de frontera, que está en dx/4).  Para nodos interiores
// centroide == cara, así que no hay diferencia.  Para k=1 (pared Oeste): x=0.
// Igualmente v se muestrea en y = pMesh.yvc[l-1].
// =============================================================================
static double runDivergence(int N)
{
    const GridData pMesh = makeUniformPMesh(N, N);
    const GridData uMesh = makeUMesh(pMesh);
    const GridData vMesh = makeVMesh(pMesh);
    DivergenceOperator D(pMesh, uMesh, vMesh);

    // Muestrear u en posición de CARA: x_face = pMesh.xvc[k-1]
    Vector u(uMesh.num_active_nodes), v(vMesh.num_active_nodes);
    for (int j = 1; j <= N; ++j) {
        for (int k = 1; k <= N + 1; ++k) {
            const int id  = (j - 1) * (N + 1) + (k - 1);
            const double xf = pMesh.xvc[k - 1];   // cara = pMesh.xvc[k-1]
            const double yc = pMesh.y[j];           // y centro P = y centro U
            u[id] = std::sin(2.0 * PI * xf) * std::cos(PI * yc);
        }
    }
    for (int l = 1; l <= N + 1; ++l) {
        for (int i = 1; i <= N; ++i) {
            const int id  = (l - 1) * N + (i - 1);
            const double xc = pMesh.x[i];            // x centro P = x centro V
            const double yf = pMesh.yvc[l - 1];      // cara = pMesh.yvc[l-1]
            v[id] = std::cos(PI * xc) * std::sin(2.0 * PI * yf);
        }
    }

    const Vector div_h = D.apply(u, v);

    // Comparar con ∇·(u,v) en centros de celda P
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
// u(x,y) = sin(πy)   (independiente de x  →  ∂²u/∂x²=0 exactamente)
//   ν∇²u = ν·d²u/dy² = −νπ²sin(πy),   ν=1
//
// Al ser constante en x, la no-uniformidad del espaciado en x (debida al nodo
// de frontera k=1 en dx/4 en lugar de x=0) no afecta al resultado.
// El término en y usa el método imagen (celda fantasma) en j=1 y j=M, que
// es de 2.º orden.
//
// Nota: sin(πy)=0 en y=0 y y=1 → u_ghost = 2·0 − u_P = −u_P, que es el
// valor exacto del campo en la celda espejo (sin(π·(−dy/2)) = −sin(π·dy/2) = −u_P). ✓
// =============================================================================
static double runDiffusion(int N)
{
    const GridData pMesh = makeUniformPMesh(N, N);
    const GridData uMesh = makeUMesh(pMesh);
    const GridData vMesh = makeVMesh(pMesh);

    const double nu = 1.0;
    // sin(πy)=0 en paredes Sur y Norte → ghost=0
    DiffusionOperator Diff(pMesh, uMesh, vMesh, nu,
                           /*u_ghost_N=*/0.0, /*u_ghost_S=*/0.0,
                           /*v_ghost_W=*/0.0, /*v_ghost_E=*/0.0);

    // Muestrear u = sin(πy) en centros de nodos U
    Vector u(uMesh.num_active_nodes);
    for (int id = 0; id < uMesh.num_active_nodes; ++id) {
        u[id] = std::sin(PI * uMesh.cells[id].y);
    }

    Vector diff_u;
    Diff.applyU(u, diff_u);

    // Comparar en nodos U interiores k=2..N, j=1..M
    // (j=1 y j=M usan celda fantasma de 2.º orden ✓)
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
// Suite 5 — Identidad discreta D·G + L/V = 0
//
// Para cualquier campo φ y para TODAS las celdas (incluidas las de frontera):
//   D·(G·φ)  +  (L·φ)/V  =  0   a precisión máquina
//
// Demostración para cara interior E de celda (i,j):
//   Contribución de G·L al divergente: Se/V * (φ[i+1,j]−φ[i,j]) / (x[i+1]−x[i]) = +a_E/V*(φ_E−φ_P)
//   Contribución de L/V:  [−a_E*(φ_E−φ_P) + ...] / V = −a_E/V*(φ_E−φ_P)
//   Suma = 0 ✓
//
// Para cara de frontera (bc=true): G devuelve 0 y L omite ese término → ambos = 0 ✓
// =============================================================================
static double runDGL(int N)
{
    const GridData pMesh = makeUniformPMesh(N, N);
    const GridData uMesh = makeUMesh(pMesh);
    const GridData vMesh = makeVMesh(pMesh);

    LaplaceOperator    L(pMesh);
    GradientOperator   G(pMesh, uMesh, vMesh);
    DivergenceOperator D(pMesh, uMesh, vMesh);

    // Campo φ = cos(πx)cos(πy) en centros P
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

    // D·G·φ + L·φ/V debe ser 0 (no D·G·φ - L·φ/V)
    Vector sum(pMesh.num_active_nodes);
    for (int k = 0; k < pMesh.num_active_nodes; ++k)
        sum[k] = DGphi[k] + LphiV[k];

    Vector zero(pMesh.num_active_nodes, 0.0);
    return errorLinf(sum, zero);
}

// =============================================================================
// CSV helper — añade una fila a output/convergence.csv
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
// v(x,y) = sin(πx)   (independiente de y  →  ∂²v/∂y²=0 exactamente)
//   ν∇²v = ν·d²v/dx² = −νπ²sin(πx),   ν=1
//
// Las paredes Oeste (x=0) y Este (x=1) tienen v=0: sin(π·0)=sin(π·1)=0,
// por eso v_ghost_W=0 y v_ghost_E=0 son consistentes con la función.
// Los nodos V de contorno (l=1, l=M+1) son puestos a 0 por applyV y no
// se incluyen en la comparación.
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

    // Muestrear v = sin(πx) en centros de nodos V
    Vector v(vMesh.num_active_nodes);
    for (int id = 0; id < vMesh.num_active_nodes; ++id)
        v[id] = std::sin(PI * vMesh.cells[id].x);

    Vector diff_v;
    Diff.applyV(v, diff_v);

    // Comparar en nodos V interiores: l=2..M (todos los i=1..N)
    // (i=1 e i=N usan celda fantasma — verifican la lógica O/E)
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
// Suite 7 — DiffusionOperator con condición de contorno no nula (tapa deslizante)
//
// u(x,y) = y   (lineal en y  →  ∂²u/∂y²=0  →  ν∇²u=0 exactamente)
//
// Pared Sur:  u_wall=0  →  ghost = 2·0 − u_P = −u_P
// Pared Norte: u_wall=1 →  ghost = 2·1 − u_P
//
// Para un campo lineal, la diferencia de segundo orden cancela exactamente
// (demostración analítica en OPERATORS_REPORT.md, Sección 4, Suite 7).
// El error esperado es < 1e-12 a todas las N (cancelación algebraica, no O(h²)).
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

    // Muestrear u = y en centros de nodos U
    Vector u(uMesh.num_active_nodes);
    for (int id = 0; id < uMesh.num_active_nodes; ++id)
        u[id] = uMesh.cells[id].y;

    Vector diff_u;
    Diff.applyU(u, diff_u);

    // Comparar en nodos U interiores k=2..N, todos j
    // Resultado exacto es 0 para campo lineal
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
// Suite 8 — Adjunto discreto: Gᵀ = −D  (precisión máquina)
//
// Propiedad de simetría-preservación:
//   φᵀ·(D·u)  +  (G·φ)ᵀ·u  =  0   ∀ φ∈P, u∈U (nodos int.)
//
// Es la forma discreta de la identidad de integración por partes:
//   ∫(∇·u)φ dV + ∫u·∇φ dV = ∮φu·n dS = 0   (con u·n=0 en paredes)
//
// Se repite para el componente v.  Error esperado: ~1e-13.
// =============================================================================
static double runAdjoint(int N)
{
    const GridData pMesh = makeUniformPMesh(N, N);
    const GridData uMesh = makeUMesh(pMesh);
    const GridData vMesh = makeVMesh(pMesh);

    GradientOperator   G(pMesh, uMesh, vMesh);
    DivergenceOperator D(pMesh, uMesh, vMesh);

    // Vectores aleatorios reproducibles (semilla fija)
    std::mt19937 rng(42);
    std::uniform_real_distribution<double> dist(-1.0, 1.0);

    // Campo φ en P (cualquier valor)
    Vector phi(pMesh.num_active_nodes);
    for (auto& v : phi) v = dist(rng);

    // Campo u en U — nodos de frontera (k=1, k=N+1) a cero para anular el
    // término de contorno en la identidad de integración por partes
    Vector u_vec(uMesh.num_active_nodes, 0.0);
    const int Np1 = N + 1;
    for (int j = 1; j <= N; ++j)
        for (int k = 2; k <= N; ++k)
            u_vec[(j - 1) * Np1 + (k - 1)] = dist(rng);

    // Campo v en V — nodos de frontera (l=1, l=M+1) a cero
    Vector v_vec(vMesh.num_active_nodes, 0.0);
    for (int l = 2; l <= N; ++l)
        for (int i = 1; i <= N; ++i)
            v_vec[(l - 1) * N + (i - 1)] = dist(rng);

    // φᵀ·(D·u) — D.apply() devuelve flux/V, por lo que el producto
    // sin ponderación de volumen es equivalente al producto en el espacio P puro.
    // La identidad D^T = -G se verifica en la norma no ponderada:
    //   φᵀ·(D·u) + (G·φ)ᵀ·u = 0
    const Vector Du = D.apply(u_vec, v_vec);
    double phiDu = 0.0;
    for (int k = 0; k < pMesh.num_active_nodes; ++k) {
        phiDu += phi[k] * Du[k];  // sin ponderación: D ya divide por V internamente
    }

    // (G·φ)ᵀ·u — gradiente de φ sobre nodos U/V interiores
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
// Suite 9 — LaplaceOperator: simetría y definición positiva
//
// Test A (simetría):  |φᵀ·(A·ψ) − ψᵀ·(A·φ)| < ε_maquina
// Test B (def. pos.): φᵀ·A·φ > 0  para φ aleatorio no nulo
//
// Retorna el error de simetría (Test A); la definición positiva se imprime.
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

    // Test B: φᵀ·Aφ  (pinear referencia para eliminar modo nulo cte.)
    // Fijamos phi[0]=0 para que la constante no contribuya
    Vector phi2 = phi;
    phi2[0] = 0.0;
    const Vector Aphi2 = L.apply(phi2);
    double quad = 0.0;
    for (int k = 0; k < pMesh.num_active_nodes; ++k)
        quad += phi2[k] * Aphi2[k];

    return { sym_err, quad };
}

// =============================================================================
// Suite 10 — LaplaceOperator en malla estirada (tanh, γ=2)
//
// Misma función manufacturada que Suite 1: φ=cos(πx)cos(πy)
// Malla generada con la transformación hiperbólica-tangente del artículo:
//   x_map = 1/2 · (1 + tanh(γ·(2ξ−1)) / tanh(γ)),  γ=2
//
// NOTA sobre convergencia:
// El operador FVM usa a_E = Se / (x[i+1] − x[i]) donde x[i] son centros de celda.
// En mallas no uniformes (h_e ≠ h_w), el término cruzado φ'·(h_e − h_w)/dx_i
// introduce un error LOCAL de orden O(1/N) (primer orden), ya que:
//   h_e − h_w = f''(ξ)/N² + O(N⁻³)   y   dx_i = f'(ξ)/N + O(N⁻²)
//   → error = φ'·f''(ξ)/(f'(ξ)·N) = O(1/N)
//
// Con el factor correcto 2/(h_e + h_w) el esquema recuperaría orden 2.
// Este test expone la limitación del estencil actual sobre mallas estiradas.
// Tasa esperada: ~1.0
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

    // Posiciones de cara en x (0..N) con estiramiento tanh
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

    // Centros de celda (promedio de caras)
    g.x.resize(N + 2);
    for (int i = 1; i <= N; ++i) g.x[i] = 0.5 * (g.xvc[i - 1] + g.xvc[i]);
    g.x[0]     = 2.0 * g.xvc[0] - g.x[1];
    g.x[N + 1] = 2.0 * g.xvc[N] - g.x[N];

    g.y.resize(M + 2);
    for (int j = 1; j <= M; ++j) g.y[j] = 0.5 * (g.yvc[j - 1] + g.yvc[j]);
    g.y[0]     = 2.0 * g.yvc[0] - g.y[1];
    g.y[M + 1] = 2.0 * g.yvc[M] - g.y[M];

    // Áreas y volúmenes
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

    // Celdas y conectividad
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
// main — ejecuta las 5 suites e imprime tablas de convergencia
// =============================================================================
int main()
{
    const std::vector<int> sizes = { 8, 16, 32, 64 };
    const int W = 14;

    // ---- Abrir CSV de salida ----
    // Crear directorio output si no existe
    std::system("mkdir -p output");
    std::ofstream csv("output/convergence.csv");
    csv << "suite,N,error_linf\n";

    auto rate = [](double e_coarse, double e_fine) -> std::string {
        if (e_coarse < 1e-300 || e_fine < 1e-300) return "  —  ";
        double r = std::log2(e_coarse / e_fine);
        char buf[16];
        std::snprintf(buf, sizeof(buf), "%.3f", r);
        return std::string(buf);
    };

    std::cout << std::fixed << std::setprecision(4);

    // ---- Suite 1: Laplaciano ----
    std::cout << "\n=== Suite 1: LaplaceOperator  φ=cos(πx)cos(πy)  (A·φ)/V ≈ +2π²cos(πx)cos(πy) ===\n";
    std::cout << std::setw(6)  << "N"
              << std::setw(W)  << "Error L∞"
              << std::setw(W)  << "Tasa\n";
    double prev1 = 0.0;
    for (int N : sizes) {
        const double e = runLaplace(N);
        std::cout << std::setw(6) << N
                  << std::setw(W) << std::scientific << e
                  << std::setw(W) << rate(prev1, e) << "\n";
        csvRow(csv, "Laplace_uniform", N, e);
        prev1 = e;
    }

    // ---- Suite 2: Gradiente ----
    std::cout << "\n=== Suite 2: GradientOperator  φ=cos(πx)cos(πy) ===\n";
    std::cout << std::setw(6)    << "N"
              << std::setw(W+2)  << "Error ∂p/∂x"
              << std::setw(W)    << "Tasa"
              << std::setw(W+2)  << "Error ∂p/∂y"
              << std::setw(W)    << "Tasa\n";
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

    // ---- Suite 3: Divergencia ----
    std::cout << "\n=== Suite 3: DivergenceOperator  u=sin(2πx)cos(πy), v=cos(πx)sin(2πy) ===\n";
    std::cout << "    (muestreo en posicion de cara: x=pMesh.xvc[k-1], y=pMesh.yvc[l-1])\n";
    std::cout << std::setw(6) << "N"
              << std::setw(W) << "Error L∞"
              << std::setw(W) << "Tasa\n";
    double prev3 = 0.0;
    for (int N : sizes) {
        const double e = runDivergence(N);
        std::cout << std::setw(6) << N
                  << std::setw(W) << std::scientific << e
                  << std::setw(W) << rate(prev3, e) << "\n";
        csvRow(csv, "Divergence", N, e);
        prev3 = e;
    }

    // ---- Suite 4: Difusión U ----
    std::cout << "\n=== Suite 4: DiffusionOperator  u=sin(πy) [cte. en x], ν=1 ===\n";
    std::cout << "    (la independencia en x elimina el error de 1er orden en nodos k=2, k=N)\n";
    std::cout << std::setw(6) << "N"
              << std::setw(W) << "Error L∞"
              << std::setw(W) << "Tasa\n";
    double prev4 = 0.0;
    for (int N : sizes) {
        const double e = runDiffusion(N);
        std::cout << std::setw(6) << N
                  << std::setw(W) << std::scientific << e
                  << std::setw(W) << rate(prev4, e) << "\n";
        csvRow(csv, "Diffusion_U", N, e);
        prev4 = e;
    }

    // ---- Suite 5: Identidad D·G + L/V = 0 ----
    std::cout << "\n=== Suite 5: Identidad D·G·φ + (L·φ)/V = 0  (precision maquina) ===\n";
    std::cout << std::setw(6) << "N"
              << std::setw(W+2) << "max|D·G·φ + L·φ/V|\n";
    for (int N : sizes) {
        const double e = runDGL(N);
        std::cout << std::setw(6) << N
                  << std::setw(W+2) << std::scientific << e << "\n";
        csvRow(csv, "DGL_identity", N, e);
    }

    std::cout << "\nTasas esperadas: ~2.0 para Suites 1-4, ~1e-13 para Suite 5.\n";

    // =========================================================================
    // NUEVAS SUITES (6-10)
    // =========================================================================

    // ---- Suite 6: Difusión V ----
    std::cout << "\n=== Suite 6: DiffusionOperator applyV  v=sin(πx) [cte. en y], ν=1 ===\n";
    std::cout << "    (verifica celdas fantasma Oeste/Este, ruta v-component independiente)\n";
    std::cout << std::setw(6) << "N"
              << std::setw(W) << "Error L∞"
              << std::setw(W) << "Tasa\n";
    double prev6 = 0.0;
    for (int N : sizes) {
        const double e = runDiffusionV(N);
        std::cout << std::setw(6) << N
                  << std::setw(W) << std::scientific << e
                  << std::setw(W) << rate(prev6, e) << "\n";
        csvRow(csv, "Diffusion_V", N, e);
        prev6 = e;
    }

    // ---- Suite 7: Difusión con BC no nula (tapa deslizante) ----
    std::cout << "\n=== Suite 7: DiffusionOperator u=y (lineal), u_N=1, u_S=0  →  ν∇²u = 0 exacto ===\n";
    std::cout << "    (cancelación algebraica para campo lineal — error esperado < 1e-12)\n";
    std::cout << std::setw(6) << "N"
              << std::setw(W+2) << "max|applyU(y)|\n";
    for (int N : sizes) {
        const double e = runDiffusionLinear(N);
        std::cout << std::setw(6) << N
                  << std::setw(W+2) << std::scientific << e << "\n";
        csvRow(csv, "Diffusion_linear_BC", N, e);
    }

    // ---- Suite 8: Adjunto Gᵀ = −D ----
    std::cout << "\n=== Suite 8: Adjunto discreto  φᵀ·(D·u) + (G·φ)ᵀ·u = 0  (precision maquina) ===\n";
    std::cout << std::setw(6) << "N"
              << std::setw(W+2) << "Residuo adjunto\n";
    for (int N : sizes) {
        const double e = runAdjoint(N);
        std::cout << std::setw(6) << N
                  << std::setw(W+2) << std::scientific << e << "\n";
        csvRow(csv, "Adjoint_GD", N, e);
    }

    // ---- Suite 9: SPD de LaplaceOperator ----
    std::cout << "\n=== Suite 9: LaplaceOperator SPD — simetría y definición positiva ===\n";
    std::cout << std::setw(6) << "N"
              << std::setw(W+2) << "Err simetría"
              << std::setw(W+2) << "φᵀAφ (>0?)\n";
    bool all_pd = true;
    for (int N : sizes) {
        auto [sym, quad] = runSPD(N);
        std::cout << std::setw(6) << N
                  << std::setw(W+2) << std::scientific << sym
                  << std::setw(W+2) << std::scientific << quad;
        if (quad <= 0.0) { std::cout << "  *** FALLO PD ***"; all_pd = false; }
        std::cout << "\n";
        csvRow(csv, "SPD_symmetry", N, sym);
    }
    if (all_pd)
        std::cout << "    Definición positiva: OK en todos los N.\n";

    // ---- Suite 10: Laplaciano en malla estirada ----
    std::cout << "\n=== Suite 10: LaplaceOperator en malla tanh (γ=2)  φ=cos(πx)cos(πy) ===\n";
    std::cout << "    (estencil a_E=Se/(x[i+1]-x[i]) es 1.er orden en malla no uniforme)\n";
    std::cout << "    (tasa esperada ~1.0; para 2.º orden usar factor 2/(h_e+h_w))\n";
    std::cout << std::setw(6)  << "N"
              << std::setw(W)  << "Error L∞"
              << std::setw(W)  << "Tasa\n";
    double prev10 = 0.0;
    for (int N : sizes) {
        const double e = runLaplaceStretched(N);
        std::cout << std::setw(6) << N
                  << std::setw(W) << std::scientific << e
                  << std::setw(W) << rate(prev10, e) << "\n";
        csvRow(csv, "Laplace_stretched", N, e);
        prev10 = e;
    }

    std::cout << "\nTasas esperadas: ~2.0 para Suites 1-4, 6; ~1.0 para Suite 10 (malla no uniforme).\n";
    std::cout << "Errores esperados < 1e-12: Suites 5, 7, 8, 9 (identidades algebraicas).\n";
    std::cout << "\nDatos de convergencia escritos en: output/convergence.csv\n";

    return 0;
}
