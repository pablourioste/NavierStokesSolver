#include "mesh/Mesh_config.h"
#include "PressurePoissonSystem.hpp"
#include "solver/sparse/solver_class_sparse.hpp"
#include "MomentumSolver.hpp"
#include "EnergySolver.hpp"
#include "Buoyancy.hpp"
#include "postprocessing/ObstacleForces.hpp"
#include "mainUtils.hpp"
#include <iostream>
#include <fstream>
#include <string>
#include <algorithm>
#include <cmath>
#include <vector>
#include <sstream>
#include <filesystem>
#include <iomanip>

using namespace std;

// Results are classified as output/<scenario>/:
//   mesh/                     mesh_P/U/V.vtk (+ obstacle_mask.vtk)
//   <Re_N | Ra_N>/fields/     VTK fields for ParaView / the plotters
//   <Re_N | Ra_N>/data/       CSV profiles, forces and benchmark metrics
//   <Re_N | Ra_N>/figures/    written by plotters/*.py, never by the solver
const string OUTPUT_ROOT = "output/";

// Which physical problem main() solves this run — selected via
// input/scenario.txt. Additive: LID_DRIVEN_CAVITY and
// DIFFERENTIALLY_HEATED_CAVITY preserve every code path exactly as before
// this enum existed (see energy_on derivation below); SQUARE_CYLINDER is the
// new third branch.
enum class Scenario { LID_DRIVEN_CAVITY, DIFFERENTIALLY_HEATED_CAVITY, SQUARE_CYLINDER };

// Reads input/scenario.txt (a single "scenario = <name>" line) BEFORE any
// other config, since every other input file now lives under
// input/<scenario>/. Falls back to "lid_driven_cavity" if the file is
// missing, preserving the historical default behavior.
static std::string read_scenario_name(const std::string& filename) {
    std::ifstream file(filename);
    std::string line;
    while (std::getline(file, line)) {
        if (line.empty() || line[0] == '#') continue;
        std::stringstream ss(line);
        std::string key, eq, val;
        if (ss >> key >> eq >> val && key == "scenario") return val;
    }
    return "lid_driven_cavity";
}

// Reads the Square Cylinder's obstacle.txt ("key = value" lines) into the
// 4 bounds mark_solid_block() needs. Returns false (bounds left at 0) if the
// file is missing — only ever called for the SQUARE_CYLINDER scenario.
static bool read_obstacle_bounds(const std::string& filename,
                                  double& x0, double& x1, double& y0, double& y1) {
    std::ifstream file(filename);
    if (!file.is_open()) return false;
    std::string line;
    while (std::getline(file, line)) {
        if (line.empty() || line[0] == '#') continue;
        std::stringstream ss(line);
        std::string key, eq; double val;
        if (!(ss >> key >> eq >> val)) continue;
        if (key == "obstacle_x0") x0 = val;
        else if (key == "obstacle_x1") x1 = val;
        else if (key == "obstacle_y0") y0 = val;
        else if (key == "obstacle_y1") y1 = val;
    }
    return true;
}

static std::vector<int> read_reynolds_cases(const std::string& filename) {
    std::vector<int> cases;
    std::ifstream file(filename);
    if (file.is_open()) {
        std::string line;
        while (std::getline(file, line)) {
            if (line.empty() || line[0] == '#') continue;
            for (char& c : line) {
                if (c == ',') c = ' ';
            }
            std::stringstream ss(line);
            int re_val = 0;
            while (ss >> re_val) {
                if (re_val > 0) cases.push_back(re_val);
            }
        }
    }

    if (cases.empty()) {
        for (int re = 100; re <= 1000; re += 100) {
            cases.push_back(re);
        }
        std::cout << "No se encontraron Reynolds en '" << filename
                  << "'. Usando bateria por defecto: 100..1000 (paso 100)\n";
    }
    return cases;
}

// Reads a list of Rayleigh numbers (doubles, e.g. 1e5) for the Differentially
// Heated Cavity sweep. Falls back to the de Vahl Davis benchmark set.
static std::vector<double> read_rayleigh_cases(const std::string& filename) {
    std::vector<double> cases;
    std::ifstream file(filename);
    if (file.is_open()) {
        std::string line;
        while (std::getline(file, line)) {
            if (line.empty() || line[0] == '#') continue;
            for (char& c : line) if (c == ',') c = ' ';
            std::stringstream ss(line);
            double ra = 0.0;
            while (ss >> ra) {
                if (ra > 0.0) cases.push_back(ra);
            }
        }
    }
    if (cases.empty()) {
        cases = {1e3, 1e4, 1e5, 1e6};
        std::cout << "No se encontraron Rayleigh en '" << filename
                  << "'. Usando bateria por defecto: 1e3,1e4,1e5,1e6\n";
    }
    return cases;
}

int main() {

    // =========================================================================
    // 0. SCENARIO SELECTION — read first, drives which input/<scenario>/
    //    subfolder everything else below is read from.
    // =========================================================================
    const std::string scenario_name = read_scenario_name("input/scenario.txt");
    const std::string INPUT_DIR     = "input/" + scenario_name + "/";
    const std::string OUTPUT_DIR = OUTPUT_ROOT + scenario_name + "/";
    const std::string MESH_DIR   = OUTPUT_DIR + "mesh/";
    std::filesystem::create_directories(MESH_DIR);
    std::cout << "Active scenario: " << scenario_name
              << "  (reading from " << INPUT_DIR << ")\n";

    Scenario scenario = Scenario::LID_DRIVEN_CAVITY;
    if (scenario_name == "differentially_heated_cavity") scenario = Scenario::DIFFERENTIALLY_HEATED_CAVITY;
    // square_cylinder and its benchmark variants (e.g. square_cylinder_breuer)
    else if (scenario_name.rfind("square_cylinder", 0) == 0) scenario = Scenario::SQUARE_CYLINDER;
    // else stays LID_DRIVEN_CAVITY — matches the historical default behavior.

    // =========================================================================
    // 1. MESH GENERATION
    // =========================================================================
    MeshConfig myConfig;
    myConfig.get_data(INPUT_DIR       + "data.txt");
    myConfig.get_stretch(INPUT_DIR    + "stretch.txt");
    myConfig.get_sim_params(INPUT_DIR + "sim_params.txt");

    GridData pMesh = myConfig.generate_mesh();
    GridData uMesh = myConfig.generate_u_mesh(pMesh);
    GridData vMesh = myConfig.generate_v_mesh(pMesh);

    myConfig.define_boundaries(pMesh, INPUT_DIR + "boundaries.txt");
    myConfig.build_connectivity(pMesh);
    myConfig.build_connectivity(uMesh);
    myConfig.build_connectivity(vMesh);

    // Square Cylinder only: mark the immersed solid obstacle AFTER
    // build_connectivity (needs bc_*/n_* already derived) and BEFORE the
    // pressure matrix is assembled below (it must see is_solid already set).
    // Never called for the other two scenarios — pMesh.cells[*].is_solid
    // stays false everywhere for them, by construction.
    if (scenario == Scenario::SQUARE_CYLINDER) {
        double ox0 = 0.0, ox1 = 0.0, oy0 = 0.0, oy1 = 0.0;
        if (!read_obstacle_bounds(INPUT_DIR + "obstacle.txt", ox0, ox1, oy0, oy1)) {
            std::cerr << "Error: could not read " << INPUT_DIR << "obstacle.txt\n";
            return 1;
        }
        myConfig.mark_solid_block(pMesh, ox0, ox1, oy0, oy1);
    }

    myConfig.export_mesh_to_VTK(pMesh, MESH_DIR + "mesh_P.vtk");
    myConfig.export_mesh_to_VTK(uMesh, MESH_DIR + "mesh_U.vtk");
    myConfig.export_mesh_to_VTK(vMesh, MESH_DIR + "mesh_V.vtk");


    // =========================================================================
    // 2. CALCULATIONS  (TODO)
    // =========================================================================
    // Pre-processing steps before the solver, e.g.:
    //   - Compute initial conditions
    //   - Set up field arrays (u, v, p)
    //   - Compute diffusion/convection coefficients

    PressurePoissonSystem poisson(pMesh, uMesh, vMesh);
    auto triplets = poisson.assembleMatrix();
    // Closed-cavity scenarios (LDC, DHC) have no pressure reference (all
    // boundaries are homogeneous Neumann) and need pinReferencePressure to
    // regularize the singular system. The Square Cylinder has a real OUTLET
    // with a Dirichlet pressure row (see LaplaceOperator), which already
    // regularizes the system on its own — pinning it too would be redundant
    // and would incorrectly zero a node that isn't the designated reference.
    if (scenario != Scenario::SQUARE_CYLINDER) {
        // pinReferencePressure needs a dummy b here; b[0]=0 is re-applied each step.
        Vector dummy_b(pMesh.num_active_nodes, 0.0);
        poisson.pinReferencePressure(triplets, dummy_b);
    }

    CSRMatrix A = CSRMatrix::fromTriplets(pMesh.num_active_nodes, triplets);
    PCG PCG_solver;
    CG CG_solver;
    GS GS_solver;
    SOR SOR_solver;

    // Simulation parameters (base from input/sim_params.txt):
    double U_lid     = myConfig.U_lid;
    double dt        = myConfig.dt;
    double rho       = myConfig.rho;
    int    max_steps = myConfig.max_steps;
    double ss_tol    = myConfig.ss_tol;
    std::filesystem::create_directories(OUTPUT_DIR);

    // --- Problem selection: lid-driven cavity vs differentially heated cavity ---
    // energy_on keeps its exact historical meaning (drives EnergySolver +
    // Buoyancy below) — now cross-derived from `scenario` instead of
    // solve_energy directly, but still requires solve_energy!=0 so a
    // misconfigured sim_params.txt can't silently enable the energy equation.
    const bool   energy_on = (scenario == Scenario::DIFFERENTIALLY_HEATED_CAVITY) &&
                             (myConfig.solve_energy != 0);
    const double Pr        = myConfig.Pr;
    const double U_inf     = myConfig.U_inf;

    // Unified case list: Reynolds numbers (lid / cylinder) or Rayleigh numbers (heated).
    std::vector<double> case_values;
    if (scenario == Scenario::DIFFERENTIALLY_HEATED_CAVITY) {
        case_values = read_rayleigh_cases(INPUT_DIR + "rayleigh_cases.txt");
        myConfig.define_temperature_boundaries(pMesh, INPUT_DIR + "boundaries_T.txt");
        std::cout << "\n*** Energy equation ON — Differentially Heated Cavity "
                  << "(Ra sweep, Pr = " << Pr << ") ***\n";
    } else {
        // Lid-Driven Cavity and Square Cylinder both sweep Reynolds numbers
        // from the same per-scenario reynolds_cases.txt file.
        for (int re : read_reynolds_cases(INPUT_DIR + "reynolds_cases.txt"))
            case_values.push_back(static_cast<double>(re));
    }


    // =========================================================================
    // 4. POST-PROCESSING — export solution fields to VTK
    // =========================================================================

     // Writes a scalar field (one value per cell) on a rectilinear grid to VTK.
    auto writeFieldVTK = [](const GridData& grid,
                            const Vector& field,
                            const std::string& fieldName,
                            const std::string& filename)
    {
        const int N = grid.N_cells_x;
        const int M = grid.M_cells_y;

        std::ofstream file(filename);
        if (!file.is_open()) {
            std::cerr << "Error: could not open " << filename << "\n";
            return;
        }

        file << "# vtk DataFile Version 2.0\n";
        file << fieldName << "\n";
        file << "ASCII\n";
        file << "DATASET RECTILINEAR_GRID\n";
        file << "DIMENSIONS " << (N + 1) << " " << (M + 1) << " 1\n";

        file << "X_COORDINATES " << (N + 1) << " double\n";
        for (int i = 0; i <= N; ++i) file << grid.xvc[i] << " ";
        file << "\n";

        file << "Y_COORDINATES " << (M + 1) << " double\n";
        for (int j = 0; j <= M; ++j) file << grid.yvc[j] << " ";
        file << "\n";

        file << "Z_COORDINATES 1 double\n0.0\n";

        file << "\nCELL_DATA " << (N * M) << "\n";
        file << "SCALARS " << fieldName << " double 1\n";
        file << "LOOKUP_TABLE default\n";
        // VTK cell order: j-major (bottom row first), i minor
        for (int j = 1; j <= M; ++j)
            for (int i = 1; i <= N; ++i)
                file << field[(j - 1) * N + (i - 1)] << "\n";

        file.close();
        std::cout << "Field '" << fieldName << "' exported to: " << filename << "\n";
    };

    if (scenario == Scenario::SQUARE_CYLINDER) {
        // Debug visual check (validation checkpoint 1): export the solid
        // mask as its own cell-data scalar, independent of any flow case, so
        // the obstacle footprint can be confirmed visually in ParaView.
        Vector mask(pMesh.num_active_nodes, 0.0);
        for (int k = 0; k < pMesh.num_active_nodes; ++k)
            mask[k] = pMesh.cells[k].is_solid ? 1.0 : 0.0;
        writeFieldVTK(pMesh, mask, "obstacle_mask", MESH_DIR + "obstacle_mask.vtk");
    }

    for (double caseVal : case_values) {
        // Lid mode: caseVal = Re, nu = U*L/Re. Heated mode: caseVal = Ra, nu = Pr.
        // Cylinder mode: caseVal = Re, nu = U_inf*obstacle_L/Re.
        const int    Re = energy_on ? 0 : static_cast<int>(std::lround(caseVal));
        const double Ra = energy_on ? caseVal : 0.0;
        double nu;
        if (energy_on) {
            nu = Pr;
        } else if (scenario == Scenario::SQUARE_CYLINDER) {
            const double obstacle_L = pMesh.obstacle_xmax - pMesh.obstacle_xmin;
            nu = U_inf * obstacle_L / caseVal;
        } else {
            nu = U_lid * myConfig.L_domain / caseVal;
        }

        cout << "\n============================================================\n";
        if (energy_on)
            cout << "Running case Ra = " << Ra << "  Pr = " << Pr
                 << "  nu(=Pr) = " << nu << "  dt = " << dt
                 << "  max_steps = " << max_steps << "\n";
        else
            cout << "Running case Re = " << Re
                 << "  U_lid = " << U_lid
                 << "  nu = " << nu
                 << "  dt = " << dt
                 << "  max_steps = " << max_steps << "\n";

        // Initial fields per case
        Vector u_star(uMesh.num_active_nodes, 0.0);
        Vector v_star(vMesh.num_active_nodes, 0.0);
        Vector p_sol(pMesh.num_active_nodes, 0.0);
        Vector u_old(uMesh.num_active_nodes, 0.0);
        Vector b_rhs(pMesh.num_active_nodes, 0.0);

        MomentumSolver momentum(pMesh, uMesh, vMesh, nu);
        momentum.applyBoundaryConditions(u_star, v_star);

        // Energy field (heated cavity only). theta initialised to 0.5 (mid).
        EnergySolver energy(pMesh, uMesh, vMesh, /*alpha=*/1.0,
                            ConvectionScheme::UPWIND);
        Vector T(pMesh.num_active_nodes, 0.5);

        // Square Cylinder only: output directory created up front (needed to
        // open forces.csv before the time loop starts) and drag/lift
        // recorded incrementally every step.
        const bool is_cylinder = (scenario == Scenario::SQUARE_CYLINDER);
        std::ostringstream cyl_dir_stream;
        cyl_dir_stream << OUTPUT_DIR << "Re_" << Re;
        const std::string cyl_data_dir = cyl_dir_stream.str() + "/data/";
        ObstacleForces obstacleForces(pMesh, uMesh, vMesh, nu, rho);
        std::ofstream forcesFile;
        std::vector<double> force_time, force_Cl;
        if (is_cylinder) {
            std::filesystem::create_directories(cyl_data_dir);
            forcesFile.open(cyl_data_dir + "forces.csv");
            forcesFile << "step,time,Cd,Cl\n";
        }
        const double dynPress = 0.5 * rho * U_inf * U_inf *
                                 (pMesh.obstacle_xmax - pMesh.obstacle_xmin) * myConfig.W_depth;

        bool diverged = false;
        for (int step = 0; step < max_steps; ++step) {
            u_old = u_star;

            momentum.advance(u_star, v_star, dt);
            // Re-impose the velocity BCs on u* before the pressure solve: the
            // predictor also updates the nodes ON an immersed obstacle's
            // faces, and a non-zero u* there would enter the divergence of
            // the adjacent fluid cells as a flux the PPE cannot remove (the
            // interface is a homogeneous-Neumann pressure boundary), leaving
            // a mass-conservation error once the BCs are reset after the
            // projection. Also refreshes the zero-gradient outlet with u*.
            // A no-op for walls with fixed values (cavity scenarios).
            momentum.applyBoundaryConditions(u_star, v_star);
            if (energy_on)
                addBuoyancy(T, v_star, dt, Ra, Pr, pMesh, vMesh);  // inject into v*
            poisson.setVelocityFields(u_star, v_star);

            b_rhs = poisson.computeRHS(dt, rho);
            if (!is_cylinder) b_rhs[0] = 0.0;  // only needed when pinReferencePressure was applied

            // Square Cylinder meshes are highly stretched near the obstacle
            // (high cell aspect ratio), which raises the Laplacian's condition
            // number well beyond what plain CG handles well in the iteration
            // cap below — use Jacobi-preconditioned CG (PCG) for this
            // scenario only. LDC/DHC keep plain CG, exactly as before.
            p_sol = is_cylinder ? PCG_solver.solve(A, b_rhs, p_sol)
                                : CG_solver.solve(A, b_rhs, p_sol);
            poisson.projectVelocity(p_sol, dt, rho, u_star, v_star);
            momentum.applyBoundaryConditions(u_star, v_star);

            if (energy_on)
                energy.advance(T, u_star, v_star, dt);  // transport theta

            if (is_cylinder && dynPress > 0.0) {
                auto [Fx, Fy] = obstacleForces.computeForces(p_sol, u_star, v_star);
                const double Cd = Fx / dynPress;
                const double Cl = Fy / dynPress;
                forcesFile << step << "," << (step * dt) << "," << Cd << "," << Cl << "\n";
                force_time.push_back(step * dt);
                force_Cl.push_back(Cl);
            }

            if (step % 100 == 0) {
                double du_max = 0.0;
                bool finite = true;
                for (int i = 0; i < (int)u_star.size(); ++i) {
                    const double d = std::abs(u_star[i] - u_old[i]);
                    if (!std::isfinite(d)) finite = false;   // std::max would drop a NaN
                    du_max = std::max(du_max, d);
                }
                double ss_err = du_max / dt;
                const std::string tag = energy_on ? ("Ra=" + std::to_string((long long)Ra))
                                                  : ("Re=" + std::to_string(Re));
                cout << tag << "  Step " << step
                     << "  ||du||_inf/dt = " << ss_err << "\n";
                if (!finite) {
                    std::cerr << "ERROR: " << tag << " diverged (non-finite velocity) at step "
                              << step << " — reduce dt. No results written for this case.\n";
                    diverged = true;
                    break;
                }
                if (step > 0 && ss_err < ss_tol) {
                    cout << "Steady state reached at step " << step
                         << " for " << tag << "\n";
                    break;
                }
            }
        }

        if (is_cylinder) forcesFile.close();
        if (diverged) continue;   // never write NaN fields as if they were results

        std::ostringstream re_dir;
        if (energy_on) re_dir << OUTPUT_DIR << "Ra_" << (long long)Ra;
        else           re_dir << OUTPUT_DIR << "Re_" << Re;
        const std::string fields_dir = re_dir.str() + "/fields/";
        const std::string case_dir   = re_dir.str() + "/data/";   // CSV outputs
        std::filesystem::create_directories(fields_dir);
        std::filesystem::create_directories(case_dir);

        writeFieldVTK(pMesh, p_sol,   "pressure",   fields_dir + "pressure.vtk");
        writeFieldVTK(uMesh, u_star,  "u_velocity", fields_dir + "u_velocity.vtk");
        writeFieldVTK(vMesh, v_star,  "v_velocity", fields_dir + "v_velocity.vtk");
        if (energy_on)
            writeFieldVTK(pMesh, T,   "temperature", fields_dir + "temperature.vtk");

        if (is_cylinder) {
            const double St = estimate_shedding_frequency(force_time, force_Cl) *
                               (pMesh.obstacle_xmax - pMesh.obstacle_xmin) / U_inf;
            std::ofstream sf(case_dir + "strouhal.csv");
            sf << "Re,St\n" << Re << "," << St << "\n";
            std::cout << "\n--- Square Cylinder  Re=" << Re << " ---\n"
                       << "  Strouhal (peak-detection estimate, UNVALIDATED against "
                       << "literature this round) = " << St << "\n";
        }

    // =========================================================================
    // 5. VALIDATION — Ghia et al. (1982) benchmark, Re = 100
    //    u(y) along vertical centreline  (x = L/2)
    //    v(x) along horizontal centreline (y = H/2)
    //
    //    Ghia reference data normalised: y/H, u/U_lid  |  x/L, v/U_lid
    //    Make sure Re = U_lid*L/nu = 100 (set nu = L/100 in data.txt)
    // =========================================================================
        if (!energy_on && !is_cylinder) {
        const Vector& u_final = u_star;
        const Vector& v_final = v_star;

        const int Nu = uMesh.N_cells_x;   // columns in U mesh (= N+1 for staggered)
        const int Mu = uMesh.M_cells_y;
        const int Nv = vMesh.N_cells_x;
        const int Mv = vMesh.M_cells_y;

        // Column index closest to x = L/2 in the U mesh
        const double x_mid = myConfig.L_domain * 0.5;
        int ic_u = 0;
        double best = 1e30;
        for (int k = 0; k < Nu; ++k) {
            double dx = std::abs(uMesh.cells[k].x - x_mid);
            if (dx < best) { best = dx; ic_u = k % Nu; }
        }

        // Row index closest to y = H/2 in the V mesh
        const double y_mid = myConfig.H_domain * 0.5;
        int jc_v = 0;
        best = 1e30;
        for (int l = 0; l < Mv; ++l) {
            double dy = std::abs(vMesh.cells[l * Nv].y - y_mid);
            if (dy < best) { best = dy; jc_v = l; }
        }

        // --- u(y) along x = L/2 ---
        {
            std::ofstream f(case_dir + "centerline_u.csv");
            f << "y,y_norm,u,u_norm\n";
            for (int j = 0; j < Mu; ++j) {
                int idx = j * Nu + ic_u;
                double y     = uMesh.cells[idx].y;
                double u_val = u_final[idx];
                f << y << ","
                  << y / myConfig.H_domain << ","
                  << u_val << ","
                  << u_val             // already normalised if U_lid=1
                  << "\n";
            }
            std::cout << "Validation: u(y) centreline -> " << case_dir << "centerline_u.csv\n";
        }

        // --- v(x) along y = H/2 ---
        {
            std::ofstream f(case_dir + "centerline_v.csv");
            f << "x,x_norm,v,v_norm\n";
            for (int i = 0; i < Nv; ++i) {
                int idx = jc_v * Nv + i;
                double x     = vMesh.cells[idx].x;
                double v_val = v_final[idx];
                f << x << ","
                  << x / myConfig.L_domain << ","
                  << v_val << ","
                  << v_val
                  << "\n";
            }
            std::cout << "Validation: v(x) centreline -> " << case_dir << "centerline_v.csv\n";
        }

        // --- Ghia et al. (1982) reference — selected Re column ---
        // Full tables: input/ghia_reference.txt
        // Re columns available: 100, 400, 1000, 3200, 5000, 7500, 10000
        const int ghia_re_list[]   = {100, 400, 1000, 3200, 5000, 7500, 10000};
        const int ghia_ncols       = 7;
        int ghia_col = 0;
        {   // find nearest Re in list
            double best_err = 1e18;
            for (int c = 0; c < ghia_ncols; ++c) {
                double err = std::abs(static_cast<double>(Re) - ghia_re_list[c]);
                if (err < best_err) { best_err = err; ghia_col = c; }
            }
        }
        std::cout << "\n--- Ghia et al. Re=" << ghia_re_list[ghia_col]
                  << " reference (u along y/H) ---\n";
        std::cout << "y/H      u/U_lid\n";
        // TABLE I: [row][0]=y/H, [row][1..7]=u for each Re
        const double ghia_u_table[][8] = {
            {1.0000,  1.00000,  1.00000,  1.00000,  1.00000,  1.00000,  1.00000,  1.00000},
            {0.9766,  0.84123,  0.75837,  0.65928,  0.53236,  0.48223,  0.47244,  0.47221},
            {0.9688,  0.78871,  0.68439,  0.57492,  0.48296,  0.46120,  0.47048,  0.47783},
            {0.9609,  0.73722,  0.61756,  0.51117,  0.46547,  0.45992,  0.47323,  0.48070},
            {0.9531,  0.68717,  0.55892,  0.46604,  0.45863,  0.46036,  0.47048,  0.47804},
            {0.8516,  0.23151,  0.29093,  0.33304,  0.33556,  0.33754,  0.34228,  0.34635},
            {0.7344,  0.00332,  0.16256,  0.18719,  0.20087,  0.20580,  0.20591,  0.20673},
            {0.6172, -0.13641,  0.02135,  0.05702,  0.08344,  0.08183,  0.08342,  0.08344},
            {0.5000, -0.20581, -0.11477, -0.06080, -0.03111, -0.03039, -0.03432, -0.03581},
            {0.4531, -0.21090, -0.17119, -0.10648, -0.07540, -0.07404, -0.07503, -0.07496},
            {0.2813, -0.15662, -0.32726, -0.27805, -0.23186, -0.22855, -0.23176, -0.23186},
            {0.1719, -0.10150, -0.24299, -0.38289, -0.32709, -0.33050, -0.32393, -0.32709},
            {0.1016, -0.06434, -0.14612, -0.29730, -0.38000, -0.37633, -0.38324, -0.38000},
            {0.0703, -0.04775, -0.10338, -0.22220, -0.41657, -0.41810, -0.43643, -0.41657},
            {0.0625, -0.04192, -0.09266, -0.20196, -0.42768, -0.43590, -0.43048, -0.42768},
            {0.0547, -0.03717, -0.08186, -0.18109, -0.43025, -0.43025, -0.41487, -0.43025},
            {0.0000,  0.00000,  0.00000,  0.00000,  0.00000,  0.00000,  0.00000,  0.00000},
        };
        for (auto& row : ghia_u_table)
            std::cout << row[0] << "   " << row[1 + ghia_col] << "\n";

        std::cout << "\n--- Ghia et al. Re=" << ghia_re_list[ghia_col]
                  << " reference (v along x/L) ---\n";
        std::cout << "x/L      v/U_lid\n";
        // TABLE II: [row][0]=x/L, [row][1..7]=v for each Re
        const double ghia_v_table[][8] = {
            {1.0000,  0.00000,  0.00000,  0.00000,  0.00000,  0.00000,  0.00000,  0.00000},
            {0.9688, -0.05906, -0.12146, -0.21388, -0.39017, -0.49774, -0.53858, -0.54302},
            {0.9609, -0.07391, -0.15663, -0.27669, -0.47425, -0.55069, -0.55216, -0.52987},
            {0.9531, -0.08864, -0.19254, -0.33714, -0.52357, -0.55408, -0.52347, -0.49099},
            {0.9453, -0.10313, -0.22847, -0.39188, -0.54053, -0.52876, -0.48590, -0.45863},
            {0.9063, -0.16914, -0.23827, -0.51500, -0.44307, -0.41442, -0.41050, -0.41496},
            {0.8594, -0.22445, -0.44993, -0.42665, -0.37401, -0.36214, -0.36213, -0.36737},
            {0.8047, -0.24533, -0.38598, -0.31966, -0.31184, -0.30018, -0.30448, -0.30719},
            {0.5000,  0.05454,  0.05186,  0.02526,  0.00999,  0.00945,  0.00824,  0.00831},
            {0.2344,  0.17527,  0.30174,  0.32235,  0.28188,  0.27280,  0.27348,  0.27224},
            {0.2266,  0.17507,  0.30203,  0.33075,  0.29030,  0.28066,  0.28117,  0.28003},
            {0.1563,  0.16077,  0.28124,  0.37095,  0.37119,  0.35368,  0.35060,  0.35070},
            {0.0938,  0.12317,  0.22965,  0.32627,  0.42768,  0.42951,  0.41824,  0.41487},
            {0.0781,  0.10890,  0.20920,  0.30353,  0.41906,  0.43648,  0.43564,  0.43124},
            {0.0703,  0.10091,  0.19713,  0.29012,  0.40917,  0.43329,  0.44030,  0.43733},
            {0.0625,  0.09233,  0.18360,  0.27485,  0.39560,  0.42447,  0.43979,  0.43983},
            {0.0000,  0.00000,  0.00000,  0.00000,  0.00000,  0.00000,  0.00000,  0.00000},
        };
        for (auto& row : ghia_v_table)
            std::cout << row[0] << "   " << row[1 + ghia_col] << "\n";
        }  // end lid-cavity (Ghia) validation
        else if (is_cylinder) {
            // Square Cylinder: forces.csv and strouhal.csv were already
            // written above (per-step recording + end-of-case estimate) —
            // nothing further needed here. See docs/scenarios/square_cylinder.md
            // for the current validation status.
        }
        else {
        // =====================================================================
        // VALIDATION — de Vahl Davis (1983) Differentially Heated Cavity
        //   u_max on the vertical mid-plane (x = L/2)
        //   v_max on the horizontal mid-plane (y = H/2)
        //   average Nusselt number on the hot (west) wall
        // =====================================================================
        const int Nu_ = uMesh.N_cells_x, Mu_ = uMesh.M_cells_y;
        const int Nv_ = vMesh.N_cells_x, Mv_ = vMesh.M_cells_y;

        // u_max along the vertical centreline x = L/2
        const double x_mid = myConfig.L_domain * 0.5;
        int ic_u = 0; double best = 1e30;
        for (int k = 0; k < Nu_; ++k) {
            double d = std::abs(uMesh.cells[k].x - x_mid);
            if (d < best) { best = d; ic_u = k; }
        }
        double u_max = 0.0;
        for (int j = 0; j < Mu_; ++j)
            u_max = std::max(u_max, std::abs(u_star[j * Nu_ + ic_u]));

        // v_max along the horizontal centreline y = H/2
        const double y_mid = myConfig.H_domain * 0.5;
        int jc_v = 0; best = 1e30;
        for (int l = 0; l < Mv_; ++l) {
            double d = std::abs(vMesh.cells[l * Nv_].y - y_mid);
            if (d < best) { best = d; jc_v = l; }
        }
        double v_max = 0.0;
        for (int i = 0; i < Nv_; ++i)
            v_max = std::max(v_max, std::abs(v_star[jc_v * Nv_ + i]));

        // Average Nusselt on the hot west wall:  Nu = (1/H) * sum_j q_j * dy_j
        // local flux q_j = (theta_wall - theta_P) / (x_P - x_wall)
        const int N = pMesh.N_cells_x;
        double Nu_avg = 0.0;
        for (int j = 1; j <= pMesh.M_cells_y; ++j) {
            const int    k     = (j - 1) * N + 0;            // P cell (i=1, j)
            const double Twall = pMesh.bound_west_T.empty() ? 1.0
                                                            : pMesh.bound_west_T[j-1].value;
            const double dist  = pMesh.cells[k].x - pMesh.xvc[0];
            const double dy    = pMesh.yvc[j] - pMesh.yvc[j-1];
            Nu_avg += (Twall - T[k]) / dist * dy;
        }
        Nu_avg /= myConfig.H_domain;

        std::cout << "\n--- de Vahl Davis benchmark  Ra=" << Ra << " ---\n";
        std::cout << "  u_max (x=L/2) = " << u_max << "\n";
        std::cout << "  v_max (y=H/2) = " << v_max << "\n";
        std::cout << "  Nu_avg (hot wall) = " << Nu_avg << "\n";

        std::ofstream f(case_dir + "dhc_metrics.csv");
        f << "Ra,Pr,u_max,v_max,Nu_avg\n"
          << Ra << "," << Pr << "," << u_max << "," << v_max << "," << Nu_avg << "\n";
        std::cout << "Validation: metrics -> " << case_dir << "dhc_metrics.csv\n";
        }  // end heated-cavity validation
    } // end case loop

    //run_mesh_inspector(myConfig, pMesh, uMesh, vMesh);

    return 0;
}
