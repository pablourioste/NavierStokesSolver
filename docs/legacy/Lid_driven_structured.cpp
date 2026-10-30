// Legacy single-file reference prototype of the lid-driven cavity solver
// (uniform staggered grid, explicit CDS momentum, Gauss-Seidel pressure).
// Kept for reference only: it is not part of the CMake build.

#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
#include <string>
#include <algorithm> // for std::max

using namespace std;

// Creates a 2D matrix of size (N+2) x (M+2) initialised to zero
vector<vector<double>> create_matrix(int N, int M) {
    return vector<vector<double>>(N + 2, vector<double>(M + 2, 0.0));
}

// =========================================================
// FUNCTION: BOUNDARY CONDITIONS (Lid-Driven Cavity)
// =========================================================
void apply_boundary_conditions(int N, int M, double U_ref,
                               vector<vector<double>>& u,
                               vector<vector<double>>& v,
                               vector<vector<double>>& P) {

    // 1. VERTICAL walls (left and right)
    for (int j = 1; j <= M; ++j) {
        // Velocity U (normal to the vertical wall): the node lies EXACTLY on the wall
        u[1][j]   = 0.0; // West wall (left)
        u[N+1][j] = 0.0; // East wall (right)

        // Velocity V (parallel to the vertical wall): the node is not on the wall.
        // Use a ghost cell so that the average at the wall is 0.
        v[0][j]   = -v[1][j]; // West ghost cell
        v[N+1][j] = -v[N][j]; // East ghost cell

        // Pressure (Neumann dp/dn = 0): zero gradient, the ghost cell copies the interior
        P[0][j]   = P[1][j];
        P[N+1][j] = P[N][j];
    }

    // 2. HORIZONTAL walls (bottom and top)
    for (int i = 1; i <= N; ++i) {
        // Velocity V (normal to the horizontal wall): the node lies EXACTLY on the wall
        v[i][1]   = 0.0; // South wall (floor)
        v[i][M+1] = 0.0; // North wall (lid)

        // Velocity U (parallel to the horizontal wall):
        // South wall (fixed floor): u = 0
        u[i][0]   = -u[i][1];

        // North wall (moving lid): u = U_ref
        // (u_ghost + u_interior) / 2 = U_ref  =>  u_ghost = 2*U_ref - u_interior
        u[i][M+1] = 2.0 * U_ref - u[i][M];

        // Pressure (Neumann dp/dn = 0)
        P[i][0]   = P[i][1];
        P[i][M+1] = P[i][M];
    }
}

// =========================================================
// FUNCTION: PREDICTOR VELOCITY (steps 1 and 2)
// =========================================================
void compute_predictor_velocity(int N, int M, double dx, double dt, double rho, double mu,
                                const vector<vector<double>>& u,
                                const vector<vector<double>>& v,
                                vector<vector<double>>& u_p,
                                vector<vector<double>>& v_p) {

    // Assumes a uniform square mesh: dx = dy
    double dVol = dx * dx; // 2D cell volume (area = dx*dy = dx^2)

    // ---------------------------------------------------------
    // 1. Staggered-X mesh (horizontal axis, velocity U)
    // ---------------------------------------------------------
    // Start at i=2 because u[1][j] is the left wall (already 0 from the BCs)
    for (int i = 2; i <= N; ++i) {
        for (int j = 1; j <= M; ++j) {

            // A. Interpolation to the faces (CDS scheme - arithmetic mean)
            double u_e = 0.5 * (u[i+1][j] + u[i][j]);   // East face
            double u_w = 0.5 * (u[i-1][j] + u[i][j]);   // West face
            double u_n = 0.5 * (u[i][j+1] + u[i][j]);   // North face
            double u_s = 0.5 * (u[i][j-1] + u[i][j]);   // South face

            // THE TRICKY PART: vertical velocity v on the North/South faces of U.
            // There is no V node at the centre of U's North face, so the 2
            // adjacent V nodes are interpolated.
            double v_n = 0.5 * (v[i-1][j+1] + v[i][j+1]);
            double v_s = 0.5 * (v[i-1][j]   + v[i][j]);

            // B. Convective term: -[m_e*u_e - m_w*u_w + m_n*u_n - m_s*u_s]
            // Note: the mass flux m_e is (rho * u_e * face_area). face_area = dx.
            double conv_x = -rho * dx * (u_e*u_e - u_w*u_w + v_n*u_n - v_s*u_s);

            // C. Diffusive term: mu * [ (u_E - u_P)/dx * A_e ... ]
            // Since A_e = dx and the distance is dx, this reduces to the classic
            // finite-difference cross stencil:
            double diff_x = mu * (u[i+1][j] + u[i-1][j] + u[i][j+1] + u[i][j-1] - 4.0 * u[i][j]);

            // D. Total operator R(u)
            double R_u = conv_x + diff_x;

            // E. Explicit time integration (Euler)
            u_p[i][j] = u[i][j] + (dt / (rho * dVol)) * R_u;
        }
    }

    // ---------------------------------------------------------
    // 2. Staggered-Y mesh (vertical axis, velocity V)
    // ---------------------------------------------------------
    // Start at j=2 because v[i][1] is the bottom wall (floor)
    for (int i = 1; i <= N; ++i) {
        for (int j = 2; j <= M; ++j) {

            // A. Interpolation to the faces (CDS)
            double v_e = 0.5 * (v[i+1][j] + v[i][j]);
            double v_w = 0.5 * (v[i-1][j] + v[i][j]);
            double v_n = 0.5 * (v[i][j+1] + v[i][j]);
            double v_s = 0.5 * (v[i][j-1] + v[i][j]);

            // THE TRICKY PART, REVERSED: velocity U on the East/West faces of V
            double u_e = 0.5 * (u[i+1][j] + u[i+1][j-1]);
            double u_w = 0.5 * (u[i][j]   + u[i][j-1]);

            // B. Convective term
            double conv_y = -rho * dx * (u_e*v_e - u_w*v_w + v_n*v_n - v_s*v_s);

            // C. Diffusive term
            double diff_y = mu * (v[i+1][j] + v[i-1][j] + v[i][j+1] + v[i][j-1] - 4.0 * v[i][j]);

            // D. Total operator R(v)
            double R_v = conv_y + diff_y;

            // E. Time integration
            v_p[i][j] = v[i][j] + (dt / (rho * dVol)) * R_v;
        }
    }
}

// =========================================================
// FUNCTION: SOLVE PRESSURE (step 3 - Poisson equation with Gauss-Seidel)
// =========================================================
void solve_pressure(int N, int M, double dx, double dt, double rho, double Tol,
                    const vector<vector<double>>& u_p,
                    const vector<vector<double>>& v_p,
                    vector<vector<double>>& P) {

    // Matrix holding the source term b_P (the divergence)
    auto b_P = create_matrix(N, M);

    // ---------------------------------------------------------
    // 1. Mesh coefficients (uniform, dx=dy)
    // ---------------------------------------------------------
    // a_E = East_area / distance_E = dx / dx = 1.0; (same for all of them)
    double a_E = 1.0;
    double a_W = 1.0;
    double a_N = 1.0;
    double a_S = 1.0;

    // ---------------------------------------------------------
    // 2. Source term b_P (divergence of u_p, v_p)
    // ---------------------------------------------------------
    double factor = rho / dt; // Note: the derivation divides by dt with a -1 sign;
                              // both are grouped here to match the matrix equation.

    for(int i = 1; i <= N; ++i) {
        for(int j = 1; j <= M; ++j) {
            // On a staggered mesh the face velocities are read directly from the matrices
            double u_east  = u_p[i+1][j];
            double u_west  = u_p[i][j];
            double v_north = v_p[i][j+1];
            double v_south = v_p[i][j];

            // Net mass flux "created" by the cell (hence the negative sign in the theory)
            // b_P = -(1/dt) * [ (rho*u_e)*Area ... ] -> with dVol simplified:
            b_P[i][j] = -factor * dx * (u_east - u_west + v_north - v_south);
        }
    }

    // ---------------------------------------------------------
    // 3. ITERATIVE GAUSS-SEIDEL SOLVER
    // ---------------------------------------------------------
    double error = 1.0;
    int iter = 0;
    int max_iter = 10000;

    while(error > Tol && iter < max_iter) {
        error = 0.0;

        for(int i = 1; i <= N; ++i) {
            for(int j = 1; j <= M; ++j) {

                // --- Initialise the coefficients for this node ---
                double aE_node = a_E; double aW_node = a_W;
                double aN_node = a_N; double aS_node = a_S;

                // --- Apply the boundary conditions (dp/dn = 0)
                // At a wall, the neighbour coefficient towards the wall is set to zero.
                if (i == 1) aW_node = 0.0; // Left wall
                if (i == N) aE_node = 0.0; // Right wall
                if (j == 1) aS_node = 0.0; // Bottom wall
                if (j == M) aN_node = 0.0; // Top wall

                // The central coefficient is always the sum of the ACTIVE neighbours
                double aP_node = aE_node + aW_node + aN_node + aS_node;

                // --- Gauss-Seidel update ---
                // a_P * P_P = a_E * P_E + a_W * P_W + a_N * P_N + a_S * P_S + b_P
                double P_new = (aE_node * P[i+1][j] +
                                aW_node * P[i-1][j] +
                                aN_node * P[i][j+1] +
                                aS_node * P[i][j-1] +
                                b_P[i][j]) / aP_node;

                // Track the maximum error over the whole mesh to decide when to stop
                double diff = abs(P_new - P[i][j]);
                if (diff > error) {
                    error = diff;
                }

                // Update the pressure IMMEDIATELY (this is what makes it Gauss-Seidel)
                P[i][j] = P_new;
            }
        }

        // --- Pressure reference condition ---
        // Since every wall has dp/dn=0, there are infinitely many solutions.
        // A central node is fixed to 0.0 to "pin" the pressure.
        double P_ref = P[N/2][M/2];
        for(int i = 1; i <= N; ++i) {
            for(int j = 1; j <= M; ++j) {
                P[i][j] = P[i][j] - P_ref; // Shift the whole field so the centre is 0
            }
        }

        // --- Update the ghost cells to satisfy dp/dn=0 outside the iterative sweep ---
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

    cout << "Gauss-Seidel converged in " << iter << " iterations. Error: " << error << endl;
}

// =========================================================
// FUNCTION: CORRECT VELOCITIES (step 4)
// =========================================================
void correct_velocities(int N, int M, double dx, double dt, double rho,
                        const vector<vector<double>>& u_p,
                        const vector<vector<double>>& v_p,
                        const vector<vector<double>>& P,
                        vector<vector<double>>& u,
                        vector<vector<double>>& v) {

    // Correct velocity U (staggered-X mesh)
    for (int i = 2; i <= N; ++i) {
        for (int j = 1; j <= M; ++j) {
            u[i][j] = u_p[i][j] - (dt / (rho * dx)) * (P[i][j] - P[i-1][j]);
        }
    }

    // Correct velocity V (staggered-Y mesh)
    for (int i = 1; i <= N; ++i) {
        for (int j = 2; j <= M; ++j) {
            v[i][j] = v_p[i][j] - (dt / (rho * dx)) * (P[i][j] - P[i][j-1]);
        }
    }
}

// =========================================================
// FUNCTION: EXPORT RESULTS TO CSV
// =========================================================
void export_results(int N, int M, double dx, double dy, double Re,
                    const vector<vector<double>>& u,
                    const vector<vector<double>>& v) {

    // Dynamic file name, e.g. "results_LDC_Re100.csv"
    string filename = "results_LDC_Re" + to_string((int)Re) + ".csv";
    ofstream file(filename);

    for (int i = 1; i <= N; ++i) {
        for (int j = 1; j <= M; ++j) {
            double x_pos = (i - 0.5) * dx;
            double y_pos = (j - 0.5) * dy;
            double u_centre = 0.5 * (u[i][j] + u[i+1][j]);
            double v_centre = 0.5 * (v[i][j] + v[i][j+1]);
            double magnitude = sqrt(u_centre*u_centre + v_centre*v_centre);

            file << x_pos << "," << y_pos << ","
                 << u_centre << "," << v_centre << "," << magnitude << "\n";
        }
    }
    file.close();
    cout << "File '" << filename << "' written!" << endl;
}

int main() {
    int N = 60; int M = 60; double L = 1.0; double dx = L / N;
    double U_ref = 1.0; double rho = 1.0;

    // --- LIST OF REYNOLDS NUMBERS TO SIMULATE ---
    vector<double> re_list = {100.0, 400.0, 1000.0};

    for (double Re : re_list) {
        cout << "\n=============================================" << endl;
        cout << "STARTING SIMULATION FOR Re = " << Re << endl;
        cout << "=============================================" << endl;

        double mu = (rho * U_ref * L) / Re;

        // RESET EVERYTHING TO ZERO FOR EACH REYNOLDS NUMBER
        auto u = create_matrix(N, M); auto v = create_matrix(N, M); auto P = create_matrix(N, M);
        auto u_p = create_matrix(N, M); auto v_p = create_matrix(N, M);
        auto u_old = u; auto v_old = v;

        apply_boundary_conditions(N, M, U_ref, u, v, P);

        double t = 0.0; double t_final = 50.0; double dt = 0.001;
        double tol_pressure = 1e-5; double tol_steady = 1e-6;
        bool steady_state = false; int global_iteration = 0;

        while (t < t_final && !steady_state) {
            global_iteration++;

            compute_predictor_velocity(N, M, dx, dt, rho, mu, u, v, u_p, v_p);
            solve_pressure(N, M, dx, dt, rho, tol_pressure, u_p, v_p, P);
            correct_velocities(N, M, dx, dt, rho, u_p, v_p, P, u, v);
            apply_boundary_conditions(N, M, U_ref, u, v, P);

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

            if (max_diff < tol_steady) steady_state = true;

            double dt_c = 0.35 * dx / (max_vel + 1e-9);
            double dt_d = 0.20 * (dx * dx) / (mu / rho);
            dt = min(dt_c, dt_d);
            t += dt;
        }

        cout << ">>> STEADY STATE (Re=" << Re << ") REACHED AT t = " << t << "s" << endl;

        // Export the file with Re in its name
        export_results(N, M, dx, dx, Re, u, v);
    }

    return 0;
}
