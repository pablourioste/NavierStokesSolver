#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
#include <algorithm> // per a la funció std::max

using namespace std;

// Funció per inicialitzar una matriu 2D de mida (N+2) x (M+2) a zero
vector<vector<double>> crearMatriu(int N, int M) {
    return vector<vector<double>>(N + 2, vector<double>(M + 2, 0.0));
}

// =========================================================
// FUNCIÓ: CONDICIONS DE CONTORN (Lid-Driven Cavity)
// =========================================================
void aplicar_condicions_contorno(int N, int M, double U_ref,
                                 vector<vector<double>>& u,
                                 vector<vector<double>>& v,
                                 vector<vector<double>>& P) {

    // 1. Parets VERTICALS (Esquerra i Dreta)
    for (int j = 1; j <= M; ++j) {
        // Velocitat U (Perpendicular a la paret vertical): El node cau EXACTAMENT a la paret
        u[1][j]   = 0.0; // Paret Oest (Esquerra)
        u[N+1][j] = 0.0; // Paret Est (Dreta)

        // Velocitat V (Paral·lela a la paret vertical): El node no està a la paret.
        // Interpolem amb el "Ghost cell" perquè la mitjana a la paret doni 0.
        v[0][j]   = -v[1][j]; // Ghost cell Oest
        v[N+1][j] = -v[N][j]; // Ghost cell Est

        // Pressió (Neumann dp/dn = 0): El gradient és nul, el ghost cell copia l'interior
        P[0][j]   = P[1][j];
        P[N+1][j] = P[N][j];
    }

    // 2. Parets HORITZONTALS (Inferior i Superior)
    for (int i = 1; i <= N; ++i) {
        // Velocitat V (Perpendicular a la paret horitzontal): El node cau EXACTAMENT a la paret
        v[i][1]   = 0.0; // Paret Sud (Sòl)
        v[i][M+1] = 0.0; // Paret Nord (Tapa)

        // Velocitat U (Paral·lela a la paret horitzontal):
        // Paret Sud (Sòl fix): u = 0
        u[i][0]   = -u[i][1];

        // Paret Nord (Tapa mòbil Lid-Driven): u = U_ref
        // (u_ghost + u_interior) / 2 = U_ref  =>  u_ghost = 2*U_ref - u_interior
        u[i][M+1] = 2.0 * U_ref - u[i][M];

        // Pressió (Neumann dp/dn = 0)
        P[i][0]   = P[i][1];
        P[i][M+1] = P[i][M];
    }
}

// =========================================================
// FUNCIÓ: CÀLCUL DE LA VELOCITAT PREDICTORA (Passos 1 i 2)
// =========================================================
void calcular_velocitat_predictora(int N, int M, double dx, double dt, double rho, double mu,
                                   const vector<vector<double>>& u,
                                   const vector<vector<double>>& v,
                                   vector<vector<double>>& u_p,
                                   vector<vector<double>>& v_p) {

    // Assumim malla quadrada uniforme: dx = dy
    double dVol = dx * dx; // Volum de la cel·la 2D (Area = dx*dy = dx^2)

    // ---------------------------------------------------------
    // 1. Malla Staggered-X (Eix horitzontal, velocitat U)
    // ---------------------------------------------------------
    // Comencem en i=2 perquè u[1][j] és la paret esquerra (que ja és 0 per BCs)
    for (int i = 2; i <= N; ++i) {
        for (int j = 1; j <= M; ++j) {

            // A. Interpolació a les cares (Esquema CDS - Mitjana aritmètica)
            double u_e = 0.5 * (u[i+1][j] + u[i][j]);   // Cara Est
            double u_w = 0.5 * (u[i-1][j] + u[i][j]);   // Cara Oest
            double u_n = 0.5 * (u[i][j+1] + u[i][j]);   // Cara Nord
            double u_s = 0.5 * (u[i][j-1] + u[i][j]);   // Cara Sud

            // EL TRENCACLOSQUES: Velocitat vertical v a les cares Nord/Sud de U
            // Com no hi ha node V al centre de la cara Nord de U, interpolem els 2 nodes V adjacents.
            double v_n = 0.5 * (v[i-1][j+1] + v[i][j+1]);
            double v_s = 0.5 * (v[i-1][j]   + v[i][j]);

            // B. Terme Convectiu: -[m_e*u_e - m_w*u_w + m_n*u_n - m_s*u_s]
            // Nota: El flux màssic m_e és (rho * u_e * Area_cara). Area_cara = dx.
            double conv_x = -rho * dx * (u_e*u_e - u_w*u_w + v_n*u_n - v_s*u_s);

            // C. Terme Difusiu: mu * [ (u_E - u_P)/dx * A_e ... ]
            // Com que A_e = dx i la distància és dx, ens queda la clàssica creu de diferències finites:
            double diff_x = mu * (u[i+1][j] + u[i-1][j] + u[i][j+1] + u[i][j-1] - 4.0 * u[i][j]);

            // D. Operador Total R(u)
            double R_u = conv_x + diff_x;

            // E. Integració temporal explícita (Euler)
            u_p[i][j] = u[i][j] + (dt / (rho * dVol)) * R_u;
        }
    }

    // ---------------------------------------------------------
    // 2. Malla Staggered-Y (Eix vertical, velocitat V)
    // ---------------------------------------------------------
    // Comencem en j=2 perquè v[i][1] és la paret inferior (terra)
    for (int i = 1; i <= N; ++i) {
        for (int j = 2; j <= M; ++j) {

            // A. Interpolació a les cares (CDS)
            double v_e = 0.5 * (v[i+1][j] + v[i][j]);
            double v_w = 0.5 * (v[i-1][j] + v[i][j]);
            double v_n = 0.5 * (v[i][j+1] + v[i][j]);
            double v_s = 0.5 * (v[i][j-1] + v[i][j]);

            // EL TRENCACLOSQUES INVERTIT: Velocitat U a les cares Est/Oest de V
            double u_e = 0.5 * (u[i+1][j] + u[i+1][j-1]);
            double u_w = 0.5 * (u[i][j]   + u[i][j-1]);

            // B. Terme Convectiu
            double conv_y = -rho * dx * (u_e*v_e - u_w*v_w + v_n*v_n - v_s*v_s);

            // C. Terme Difusiu
            double diff_y = mu * (v[i+1][j] + v[i-1][j] + v[i][j+1] + v[i][j-1] - 4.0 * v[i][j]);

            // D. Operador Total R(v)
            double R_v = conv_y + diff_y;

            // E. Integració temporal
            v_p[i][j] = v[i][j] + (dt / (rho * dVol)) * R_v;
        }
    }
}

// =========================================================
// FUNCIÓ: RESOLDRE PRESSIÓ (Pas 3 - Eq. de Poisson amb Gauss-Seidel)
// =========================================================
void resoldre_pressio(int N, int M, double dx, double dt, double rho, double Tol,
                      const vector<vector<double>>& u_p,
                      const vector<vector<double>>& v_p,
                      vector<vector<double>>& P) {

    // Matriu per guardar el terme font b_P (la divergència)
    auto b_P = crearMatriu(N, M);

    // ---------------------------------------------------------
    // 1. Calcular els coeficients de la malla (Uniforme dx=dy)
    // ---------------------------------------------------------
    // a_E = Area_Est / distancia_E = dx / dx = 1.0; (el mateix per tots)
    double a_E = 1.0;
    double a_W = 1.0;
    double a_N = 1.0;
    double a_S = 1.0;

    // ---------------------------------------------------------
    // 2. Calcular el Terme Font b_P (Divergència de u_p, v_p)
    // ---------------------------------------------------------
    double factor = rho / dt; // Nota: el document divideix per dt i posa -1,
                              // ho agrupem tot per coincidir amb l'eq matricial.

    for(int i = 1; i <= N; ++i) {
        for(int j = 1; j <= M; ++j) {
            // Com estem en malla staggered, les velocitats a les cares són directament les de les matrius
            double u_est  = u_p[i+1][j];
            double u_oest = u_p[i][j];
            double v_nord = v_p[i][j+1];
            double v_sud  = v_p[i][j];

            // Flux màssic net que "crea" la cel·la (per això és negatiu a la teoria)
            // b_P = -(1/dt) * [ (rho*u_e)*Area ... ] -> Amb dVol simplificat:
            b_P[i][j] = -factor * dx * (u_est - u_oest + v_nord - v_sud);
        }
    }

    // ---------------------------------------------------------
    // 3. SOLVER ITERATIU DE GAUSS-SEIDEL
    // ---------------------------------------------------------
    double error = 1.0;
    int iter = 0;
    int max_iter = 10000;

    while(error > Tol && iter < max_iter) {
        error = 0.0;

        for(int i = 1; i <= N; ++i) {
            for(int j = 1; j <= M; ++j) {

                // --- Inicialitzar coeficients per aquest node ---
                double aE_node = a_E; double aW_node = a_W;
                double aN_node = a_N; double aS_node = a_S;

                // --- Aplicar Condicions de Contorn (dp/dn = 0)
                // Si estem tocant una paret, anul·lem el coeficient veí cap a la paret.
                if (i == 1) aW_node = 0.0; // Paret Esquerra
                if (i == N) aE_node = 0.0; // Paret Dreta
                if (j == 1) aS_node = 0.0; // Paret Inferior
                if (j == M) aN_node = 0.0; // Paret Superior

                // El coeficient central és sempre la suma dels veïns ACTIUS
                double aP_node = aE_node + aW_node + aN_node + aS_node;

                // --- Equació de Gauss-Seidel [cite: 242] ---
                // a_P * P_P = a_E * P_E + a_W * P_W + a_N * P_N + a_S * P_S + b_P
                double P_new = (aE_node * P[i+1][j] +
                                aW_node * P[i-1][j] +
                                aN_node * P[i][j+1] +
                                aS_node * P[i][j-1] +
                                b_P[i][j]) / aP_node;

                // Calcular l'error màxim de tota la malla per saber si hem acabat
                double diff = abs(P_new - P[i][j]);
                if (diff > error) {
                    error = diff;
                }

                // Actualitzar el valor de la pressió IMEDIATAMENT (això és Gauss-Seidel)
                P[i][j] = P_new;
            }
        }

        // --- Condició de Referència per a la Pressió ---
        // Com totes les parets tenen dp/dn=0, hi ha infinites solucions.
        // Fixem un node central a 0.0 per "clavar" la pressió.
        double P_ref = P[N/2][M/2];
        for(int i = 1; i <= N; ++i) {
            for(int j = 1; j <= M; ++j) {
                P[i][j] = P[i][j] - P_ref; // Desplacem tota la matriu perquè el centre valgui 0
            }
        }

        // --- Actualitzar les ghost cells per complir dp/dn=0 a fora del bucle iteratiu ---
        for(int j = 1; j <= M; ++j) {
            P[0][j]   = P[1][j];
            P[N+1][j] = P[N][j];
        }
        for(int i = 1; i <= N; ++i) {
            P[i][0]   = P[i][1];
            P[i][M+1] = P[i][M];
        }

        iter++;
    }

    cout << "Gauss-Seidel convergit en " << iter << " iteracions. Error: " << error << endl;
}

// =========================================================
// FUNCIÓ: CORREGIR VELOCITATS (Pas 4)
// =========================================================
void corregir_velocitats(int N, int M, double dx, double dt, double rho,
                         const vector<vector<double>>& u_p,
                         const vector<vector<double>>& v_p,
                         const vector<vector<double>>& P,
                         vector<vector<double>>& u,
                         vector<vector<double>>& v) {

    // Corregir Velocitat U (Malla Staggered-X)
    for (int i = 2; i <= N; ++i) {
        for (int j = 1; j <= M; ++j) {
            u[i][j] = u_p[i][j] - (dt / (rho * dx)) * (P[i][j] - P[i-1][j]);
        }
    }

    // Corregir Velocitat V (Malla Staggered-Y)
    for (int i = 1; i <= N; ++i) {
        for (int j = 2; j <= M; ++j) {
            v[i][j] = v_p[i][j] - (dt / (rho * dx)) * (P[i][j] - P[i][j-1]);
        }
    }
}

// =========================================================
// FUNCIÓ: EXPORTAR RESULTATS A CSV
// =========================================================
#include <string> // Afegeix això a dalt de tot

void exportar_resultats(int N, int M, double dx, double dy, double Re,
                        const vector<vector<double>>& u,
                        const vector<vector<double>>& v) {

    // Creem un nom d'arxiu dinàmic, ex: "resultats_LDC_Re100.csv"
    string nom_arxiu = "resultats_LDC_Re" + to_string((int)Re) + ".csv";
    ofstream arxiu(nom_arxiu);

    for (int i = 1; i <= N; ++i) {
        for (int j = 1; j <= M; ++j) {
            double x_pos = (i - 0.5) * dx;
            double y_pos = (j - 0.5) * dy;
            double u_centre = 0.5 * (u[i][j] + u[i+1][j]);
            double v_centre = 0.5 * (v[i][j] + v[i][j+1]);
            double magnitud = sqrt(u_centre*u_centre + v_centre*v_centre);

            arxiu << x_pos << "," << y_pos << ","
                  << u_centre << "," << v_centre << "," << magnitud << "\n";
        }
    }
    arxiu.close();
    cout << "Arxiu '" << nom_arxiu << "' generat!" << endl;
}

int main() {
    int N = 60; int M = 60; double L = 1.0; double dx = L / N;
    double U_ref = 1.0; double rho = 1.0;

    // --- LLISTA DE REYNOLDS A SIMULAR ---
    vector<double> llista_Re = {100.0, 400.0, 1000.0};

    for (double Re : llista_Re) {
        cout << "\n=============================================" << endl;
        cout << "INICIANT SIMULACIO PER Re = " << Re << endl;
        cout << "=============================================" << endl;

        double mu = (rho * U_ref * L) / Re;

        // RE-INICIALITZAR TOT A ZERO PER A CADA REYNOLDS
        auto u = crearMatriu(N, M); auto v = crearMatriu(N, M); auto P = crearMatriu(N, M);
        auto u_p = crearMatriu(N, M); auto v_p = crearMatriu(N, M);
        auto u_old = u; auto v_old = v;

        aplicar_condicions_contorno(N, M, U_ref, u, v, P);

        double t = 0.0; double t_final = 50.0; double dt = 0.001;
        double Tol_Presio = 1e-5; double Tol_Steady = 1e-6;
        bool steady_state = false; int iteracio_global = 0;

        while (t < t_final && !steady_state) {
            iteracio_global++;

            calcular_velocitat_predictora(N, M, dx, dt, rho, mu, u, v, u_p, v_p);
            resoldre_pressio(N, M, dx, dt, rho, Tol_Presio, u_p, v_p, P);
            corregir_velocitats(N, M, dx, dt, rho, u_p, v_p, P, u, v);
            aplicar_condicions_contorno(N, M, U_ref, u, v, P);

            double max_diff = 0.0; double max_vel = 0.0;
            for (int i = 1; i <= N; ++i) {
                for (int j = 1; j <= M; ++j) {
                    double diff_u = abs(u[i][j] - u_old[i][j]);
                    double diff_v = abs(v[i][j] - v_old[i][j]);
                    if (diff_u > max_diff) max_diff = diff_u;
                    if (diff_v > max_diff) max_diff = diff_v;

                    if (abs(u[i][j]) > max_vel) max_vel = abs(u[i][j]);
                    if (abs(v[i][j]) > max_vel) max_vel = abs(v[i][j]);

                    u_old[i][j] = u[i][j]; v_old[i][j] = v[i][j];
                }
            }

            if (max_diff < Tol_Steady) steady_state = true;

            double dt_c = 0.35 * dx / (max_vel + 1e-9);
            double dt_d = 0.20 * (dx * dx) / (mu / rho);
            dt = min(dt_c, dt_d);
            t += dt;
        }

        cout << ">>> ESTAT ESTACIONARI (Re=" << Re << ") ASSOLIT A t = " << t << "s" << endl;

        // Exportem l'arxiu amb el Re al nom
        exportar_resultats(N, M, dx, dx, Re, u, v);
    }

    return 0;
}
