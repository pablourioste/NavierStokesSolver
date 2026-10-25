#include <iostream>
#include <fstream>
#include <string>
#include <sstream>   // stringstream, used to parse "key = value" lines
#include <algorithm> // For 'remove' and 'find_if'
#include <cctype>    // For 'isspace'
#include <vector>
#include <cmath>     // For math functions such as sinh, abs
#include <iomanip>   // For formatted console output

#include "Mesh_config.h"
using namespace std;


string bc_type_to_string(BCType t) {
    switch(t) {
        case WALL_ADIABATIC:   return "WALL (Adiabatic)";
        case WALL_FIXED_VALUE: return "WALL (Fixed Val)";
        case INLET:            return "INLET";
        case OUTLET:           return "OUTLET";
        case SYMMETRY:         return "SYMMETRY";
        default:               return "UNDEFINED";
    }
}

MeshConfig::MeshConfig() 
    : 
    // 1. Initialise dimensions
    L_domain(3.0),
    H_domain(3.0),
    W_depth(0.1),
    Re(100.0),
    U_lid(1.0),
    dt(1e-3),
    rho(1.0),
    max_steps(1000000),
    ss_tol(1e-6),
    solve_energy(0),
    Pr(0.71),

    // 2. Initialise positions (std::pair)
    // NW = North-West: min x, max y | NE = North-East: max x, max y
    // SW = South-West: min x, min y | SE = South-East: max x, min y
    node_NW(1.5, 2.0),
    node_SE(2.0, 1.5),
    node_NE(2.0, 2.0),
    node_SW(1.5, 1.5),

    // 3. Initialise stretching
    beta_x_west(2.0),
    beta_x_east(-2.0),
    beta_y_north(-2.0),
    beta_y_south(2.0),
    beta_x_center(0.0),
    beta_y_center(0.0),

    // 4. Initialise number of cells
    N_cells_east(50),
    M_cells_east(50),
    N_cells_west(50),
    M_cells_west(50),
    N_cells_center(50),
    M_cells_center(50)
    {
    // The body only reports creation: all values are set in the initialiser list
    cout << "MeshConfig object created with default parameters." << endl;
}
// ==========================================
// END OF CONSTRUCTOR
// ==========================================
/**
 * @brief Generates 1D face coordinates with stretching at an arbitrary 'focus_point'.
 * @param focus_point Physical coordinate (in [0, length]) to center the refinement.
 */




void MeshConfig::generate_stretched_coords(vector<double>& vc, int num_cells, double length, 
                               double beta,double offset) {

    if (abs(beta) < 1e-6) { // Uniform mesh case (no changes)
        cout<<"Generating uniform mesh in 1D."<<endl;
        cout<<"Offset applied: "<<offset<<endl;
        double delta = length / static_cast<double>(num_cells);
        for (int i = 0; i <= num_cells; ++i) {
            vc[i] = delta * static_cast<double>(i);
            vc[i] += offset; // Apply offset to position the segment correctly
        }
        return;
    } 
    
    double physical_focus;

    if (beta < 0.0) {
        
        // Beta value is negative so we have to refine at the beggining (x'=0).
        physical_focus=0.0;
        cout<<"Generating stretched mesh in 1D with refinement at the beginning (x=0)."<<endl;
    } else {

        // Beta value is positive so we have to refine at the end (x'=L).
        physical_focus=length;
        cout<<"Generating stretched mesh in 1D with refinement at the end (x=L)."<<endl;
    }

        cout<<"Offset applied: "<<offset<<endl;
        // 1. Convert the physical focus [0, L] to a logical focus [-1, 1]
        // Formula: xi = 2*(x/L) - 1
        double xi_focus = 2.0 * (physical_focus / length) - 1.0;

        // 2. Calculate the min/max values of the *shifted* sinh function
        // to re-normalize it.
        double val_min = sinh(beta * (-1.0 - xi_focus));
        double val_max = sinh(beta * (1.0 - xi_focus));
        double val_range = val_max - val_min;

        // Protection against division by zero if beta or the range is tiny
        if (std::abs(val_range) < 1e-10) {
             double delta = length / static_cast<double>(num_cells);
             for (int i = 0; i <= num_cells; ++i) {
                 vc[i] = delta * static_cast<double>(i);
             }
             return;
        }

        // 3. Generate the coordinates
        for (int i = 0; i <= num_cells; ++i) {
            // Uniform logical coordinate xi (no changes)
            double xi = -1.0 + 2.0 * static_cast<double>(i) / static_cast<double>(num_cells);
            
            // Current 'sinh' value, shifted by the focus
            double val_current = sinh(beta * (xi - xi_focus));
            
            // Re-normalize the value to the [-1, 1] range
            // g_xi = ( (val_current - val_min) / val_range ) * 2.0 - 1.0;
            double g_xi = ((val_current - val_min) / val_range) * 2.0 - 1.0;
            
            // Map the normalized g_xi [-1, 1] to the physical [0, L]
            // Formula: x = (L/2) * (1 + g_xi)
            vc[i] = (length / 2.0) * (1.0 + g_xi);
            vc[i] += offset; // Apply offset to position the segment correctly
        }
    
    
}


/**
 * @brief Geometry function, with refinement at (x_focus, y_focus).
 */
GridData MeshConfig::generate_mesh() const {
   
    // --- 1. Create the GridData object to return ---
    GridData data;

    // Total number of cells in each direction
    int N_cells = N_cells_west + N_cells_center + N_cells_east;
    int M_cells = M_cells_west + M_cells_center + M_cells_east;

    data.N_cells_x = N_cells;
    data.M_cells_y = M_cells;
    // Data resizing:
    data.xvc.resize(N_cells + 2); 
    data.yvc.resize(M_cells + 2);
    data.x.resize(N_cells + 2); // Including ghost cells
    data.y.resize(M_cells +2); // Including ghost cells
    data.Se.resize(N_cells + 2, vector<double>(M_cells + 2, 0.0));
    data.Sw.resize(N_cells + 2, vector<double>(M_cells + 2, 0.0));
    data.Sn.resize(N_cells + 2, vector<double>(M_cells + 2, 0.0));
    data.Ss.resize(N_cells + 2, vector<double>(M_cells + 2, 0.0));
    data.V.resize(N_cells + 2, vector<double>(M_cells + 2, 0.0));
    

    // Definition of auxialiary variables of vc[i] and vc[j]
    vector<double> xvc_west(N_cells_west + 2);
    vector<double> xvc_center(N_cells_center + 2);
    vector<double> xvc_east(N_cells_east + 2);
    vector<double> yvc_south(M_cells_west + 2);
    vector<double> yvc_center(M_cells_center + 2);
    vector<double> yvc_north(M_cells_east + 2);

    // Computation of Lengths and Heights of each region
    double L_west = (node_NW.first - 0.0);
    double L_center = (node_NE.first - node_NW.first);
    double L_east = (L_domain - node_NE.first);
    double H_south = (node_SW.second - 0.0);
    double H_center = (node_NW.second - node_SW.second);
    double H_north = (H_domain - node_NW.second);
    
    // --- 2. Generate 1D Coordinates ---
  
    // X direction
    cout<<"Generating X coordinates..."<<endl;
    generate_stretched_coords(xvc_west, N_cells_west, L_west, beta_x_west,0.0);
    generate_stretched_coords(xvc_center, N_cells_center, L_center, beta_x_center,node_NW.first);
    generate_stretched_coords(xvc_east, N_cells_east, L_east, beta_x_east,node_NE.first);
    
    // Y direction
    cout<<"Generating Y coordinates..."<<endl;
    generate_stretched_coords(yvc_south, M_cells_west, H_south, beta_y_south,0.0);
    generate_stretched_coords(yvc_center, M_cells_center, H_center, beta_y_center,node_SW.second);
    generate_stretched_coords(yvc_north, M_cells_east, H_north, beta_y_north,node_NW.second);
    
    
    // --- Merging the segments into the full grid coordinates ---
    // X direction
    for (int i = 1; i <= N_cells_west; ++i) {
        data.xvc[i] = xvc_west[i];
        //cout<<data.xvc[i]<<endl;
    }

    for (int i = 1; i <= N_cells_center; ++i) {
        data.xvc[i + N_cells_west] = xvc_center[i];
       // cout<<data.xvc[i + N_cells_west]<<endl;
    }
    for (int i = 1; i <= N_cells_east; ++i) {
        data.xvc[i + N_cells_west + N_cells_center] = xvc_east[i];
        //cout<<data.xvc[i + N_cells_west + N_cells_center]<<endl;
    }
    data.xvc[N_cells + 1] = L_domain;

    cout<<"yvc merging..."<<endl;
    // Y direction
    for (int j = 1; j <= M_cells_west; ++j) {
        data.yvc[j] = yvc_south[j];
       // cout<<data.yvc[j]<<endl;
    }
    for (int j = 1; j <= M_cells_center; ++j) {
        data.yvc[j + M_cells_west] = yvc_center[j];
       // cout<<data.yvc[j + M_cells_west]<<endl;
    }
    for (int j = 1; j <= M_cells_east; ++j) {
        data.yvc[j + M_cells_west + M_cells_center] = yvc_north[j];
       // cout<<data.yvc[j + M_cells_west + M_cells_center]<<endl;
    }
    data.yvc[M_cells + 1] = H_domain;

    // --- Compute cell-center coordinates (including ghost cells) ---
    // X direction: internal centers + mirror ghost cells outside domain
    for (int i = 1; i <= N_cells; ++i)
        data.x[i] = (data.xvc[i] + data.xvc[i-1]) / 2.0;
    data.x[0]         = 2.0 * data.xvc[0]       - data.x[1];
    data.x[N_cells+1] = 2.0 * data.xvc[N_cells] - data.x[N_cells];

    // Y direction
    for (int j = 1; j <= M_cells; ++j)
        data.y[j] = (data.yvc[j] + data.yvc[j-1]) / 2.0;
    data.y[0]         = 2.0 * data.yvc[0]       - data.y[1];
    data.y[M_cells+1] = 2.0 * data.yvc[M_cells] - data.y[M_cells];

    // --- Compute face areas (needed by the solver) ---
    for (int j = 1; j <= M_cells; ++j) {
        for (int i = 1; i <= N_cells; ++i) {
            double dy = data.yvc[j] - data.yvc[j-1];
            double dx = data.xvc[i] - data.xvc[i-1];
            data.Se[i][j] = data.Sw[i][j] = dy * W_depth; // East/West: normal in X
            data.Sn[i][j] = data.Ss[i][j] = dx * W_depth; // North/South: normal in Y
            data.V[i][j]  = dx * dy * W_depth;             // Cell volume
        }
    }

    // --- 1. GENERATE BOUNDARIES (SURFACES) ---
    cout << "Generating boundaries..." << endl;

    // West boundary (face i=0)
    for (int j = 1; j <= M_cells; ++j) {
        BoundaryFace face;
        face.x = data.xvc[0]; 
        face.y = (data.yvc[j] + data.yvc[j-1]) / 2.0; 
        face.Area = (data.yvc[j] - data.yvc[j-1]) * W_depth;
        // Connects to the first internal column (the ID is handled later)
        data.bound_west.push_back(face);
    }

    // East boundary (face i=N_cells)
    for (int j = 1; j <= M_cells; ++j) {
        BoundaryFace face;
        face.x = data.xvc[N_cells]; 
        face.y = (data.yvc[j] + data.yvc[j-1]) / 2.0; 
        face.Area = (data.yvc[j] - data.yvc[j-1]) * W_depth;
        data.bound_east.push_back(face);
    }

    // South boundary (face j=0)
    for (int i = 1; i <= N_cells; ++i) {
        BoundaryFace face;
        face.x = (data.xvc[i] + data.xvc[i-1]) / 2.0;
        face.y = data.yvc[0];
        face.Area = (data.xvc[i] - data.xvc[i-1]) * W_depth;
        data.bound_south.push_back(face);
    }

    // North boundary (face j=M_cells)
    for (int i = 1; i <= N_cells; ++i) {
        BoundaryFace face;
        face.x = (data.xvc[i] + data.xvc[i-1]) / 2.0;
        face.y = data.yvc[M_cells];
        face.Area = (data.xvc[i] - data.xvc[i-1]) * W_depth;
        data.bound_north.push_back(face);
    }

    // --- 2. GENERATE INTERNAL ELEMENTS (VOLUMES) ---
    cout << "Generating internal volume elements..." << endl;
    
    int current_id = 0;
    
    // Note: iterate ONLY over the internal fluid cells
    for (int j = 1; j <= M_cells; ++j) {
        for (int i = 1; i <= N_cells; ++i) {
            InternalCell cell;
            
            cell.global_id = current_id++; // 0, 1, 2...
            
            // Centroids
            cell.x = (data.xvc[i] + data.xvc[i-1]) / 2.0;
            cell.y = (data.yvc[j] + data.yvc[j-1]) / 2.0;
            
            // Geometric properties
            double dx = data.xvc[i] - data.xvc[i-1];
            double dy = data.yvc[j] - data.yvc[j-1];
            cell.Vol = dx * dy * W_depth;
            
            // Add to the list
            data.cells.push_back(cell);
        }
    }
    
    data.num_active_nodes = current_id;
    return data;
}

GridData MeshConfig::generate_u_mesh(const GridData& pMesh) const {
    GridData uMesh;

    // 1. Basic configuration
    // The U mesh has the same Y as P, but a different X
    uMesh.M_cells_y = pMesh.M_cells_y;
    uMesh.yvc = pMesh.yvc; // Copy the exact Y faces of P

    // The U mesh has N + 1 velocity points (boundaries included),
    // so it has N + 1 control cells in X
    uMesh.N_cells_x = pMesh.N_cells_x + 1;

    // 2. Generate the new X coordinates (xvc for U)
    // The U cell limits are the CENTRES of the P cells
    uMesh.xvc.push_back(pMesh.xvc[0]); // Left limit (wall)

    for (int i = 1; i <= pMesh.N_cells_x; ++i) {
        // Average of the P faces gives the P centre
        double p_center = (pMesh.xvc[i] + pMesh.xvc[i-1]) / 2.0;
        uMesh.xvc.push_back(p_center);
    }
    
    uMesh.xvc.push_back(pMesh.xvc[pMesh.N_cells_x]); // Right limit (wall)

    // 3. Generate the elements (cells) for U
    // Note: this is required for export_mesh_to_VTK to work
    int current_id = 0;
    
    // Iterate over the new dimensions
    // j goes up to M (same as P), i goes up to N+1 (there is one extra U)
    for (int j = 1; j <= uMesh.M_cells_y; ++j) {
        for (int i = 1; i <= uMesh.N_cells_x; ++i) {
            InternalCell cell;
            cell.global_id = current_id++;

            // Limits of the current U cell
            double x_left = uMesh.xvc[i-1];
            double x_right = uMesh.xvc[i];
            double y_bottom = uMesh.yvc[j-1];
            double y_top = uMesh.yvc[j];

            // Centroid of the U cell (where the velocity u_i,j lives)
            cell.x = (x_left + x_right) / 2.0;
            cell.y = (y_bottom + y_top) / 2.0;

            // Volume
            double dx = x_right - x_left;
            double dy = y_top - y_bottom;
            cell.Vol = dx * dy * W_depth;

            uMesh.cells.push_back(cell);
        }
    }
    
    uMesh.num_active_nodes = current_id;
    return uMesh;
}

GridData MeshConfig::generate_v_mesh(const GridData& pMesh) const {
    GridData vMesh;

    // 1. Basic configuration
    // The V mesh has the same X as P, but a different Y
    vMesh.N_cells_x = pMesh.N_cells_x;
    vMesh.xvc = pMesh.xvc; // Copy the exact X faces of P

    // The V mesh has M + 1 velocity points
    vMesh.M_cells_y = pMesh.M_cells_y + 1;

    // 2. Generate the new Y coordinates (yvc for V)
    // The V cell limits are the CENTRES of the P cells
    vMesh.yvc.push_back(pMesh.yvc[0]); // South wall

    for (int j = 1; j <= pMesh.M_cells_y; ++j) {
        double p_center = (pMesh.yvc[j] + pMesh.yvc[j-1]) / 2.0;
        vMesh.yvc.push_back(p_center);
    }

    vMesh.yvc.push_back(pMesh.yvc[pMesh.M_cells_y]); // North wall

    // 3. Generate the elements (cells) for V
    int current_id = 0;

    for (int j = 1; j <= vMesh.M_cells_y; ++j) {
        for (int i = 1; i <= vMesh.N_cells_x; ++i) {
            InternalCell cell;
            cell.global_id = current_id++;

            double x_left = vMesh.xvc[i-1];
            double x_right = vMesh.xvc[i];
            double y_bottom = vMesh.yvc[j-1];
            double y_top = vMesh.yvc[j];

            cell.x = (x_left + x_right) / 2.0;
            cell.y = (y_bottom + y_top) / 2.0;

            double dx = x_right - x_left;
            double dy = y_top - y_bottom;
            cell.Vol = dx * dy * W_depth;

            vMesh.cells.push_back(cell);
        }
    }

    vMesh.num_active_nodes = current_id;
    return vMesh;
}

// Face directions, for readability
enum Direction { EAST = 0, WEST = 1, NORTH = 2, SOUTH = 3 };

int MeshConfig::get_global_id(int i, int j, int N_total) {
    return (j - 1) * N_total + (i - 1);
}

void MeshConfig::build_connectivity(GridData& data) const {

    int N = data.N_cells_x;
    int M = data.M_cells_y;
    int total = N * M;

    // Initialise the connectivity map and the inverse index maps
    data.id_to_i.resize(total);
    data.id_to_j.resize(total);
    data.connectivity.resize(total, std::vector<int>(4, -1));

    // Helper lambda giving the 1D index from 2D (i,j)
    // Lexicographic order (loop over j, then i), as in generate_mesh
    auto get_idx = [&](int i, int j) {
        return (j - 1) * N + (i - 1);
    };

    for (int j = 1; j <= M; ++j) {
        for (int i = 1; i <= N; ++i) {

            int current_idx = get_idx(i, j);
            InternalCell& cell = data.cells[current_idx];

            // Inverse maps
            data.id_to_i[current_idx] = i;
            data.id_to_j[current_idx] = j;

            // --- EAST NEIGHBOUR ---
            if (i < N) {
                cell.n_east = get_idx(i + 1, j);
                cell.bc_east = false;
            } else {
                cell.n_east = -1; // No internal neighbour
                cell.bc_east = true; // Touches the East boundary
            }
            data.connectivity[current_idx][EAST] = cell.n_east;

            // --- WEST NEIGHBOUR ---
            if (i > 1) {
                cell.n_west = get_idx(i - 1, j);
                cell.bc_west = false;
            } else {
                cell.n_west = -1;
                cell.bc_west = true; // Touches the West boundary
            }
            data.connectivity[current_idx][WEST] = cell.n_west;

            // --- NORTH NEIGHBOUR ---
            if (j < M) {
                cell.n_north = get_idx(i, j + 1);
                cell.bc_north = false;
            } else {
                cell.n_north = -1;
                cell.bc_north = true; // Touches the North boundary
            }
            data.connectivity[current_idx][NORTH] = cell.n_north;

            // --- SOUTH NEIGHBOUR ---
            if (j > 1) {
                cell.n_south = get_idx(i, j - 1);
                cell.bc_south = false;
            } else {
                cell.n_south = -1;
                cell.bc_south = true; // Touches the South boundary
            }
            data.connectivity[current_idx][SOUTH] = cell.n_south;
        }
    }

    cout << "Connectivity built. Active internal nodes: " << data.num_active_nodes << endl;
}

// ---------------------------------------------------------------------------
// mark_solid_block — immersed solid obstacle (Square Cylinder scenario)
//
// Marks every P-cell whose centroid lies inside [x0,x1]x[y0,y1] as is_solid,
// and re-derives bc_east/west/north/south (+ n_east/west/north/south = -1)
// at every fluid/solid interface, EXACTLY the way build_connectivity already
// flags the outer domain edge — so LaplaceOperator's existing
// "if (!cell.bc_east) ..." branches automatically treat the obstacle as a
// zero-flux wall with no further code changes there.
//
// Must run AFTER build_connectivity(pMesh) and BEFORE the pressure matrix is
// assembled. Never called for the Lid-Driven Cavity or DHC scenarios.
// ---------------------------------------------------------------------------
void MeshConfig::mark_solid_block(GridData& pMesh, double x0, double x1,
                                   double y0, double y1) const {
    const int N = pMesh.N_cells_x;
    const int M = pMesh.M_cells_y;

    // Snap requested bounds to the nearest existing face coordinate, so the
    // obstacle boundary coincides exactly with P-cell faces (required for an
    // exact no-slip ghost-image method at the surface).
    auto snap = [](const std::vector<double>& vc, double target) {
        double best = vc[0];
        double best_d = std::abs(vc[0] - target);
        for (double v : vc) {
            double d = std::abs(v - target);
            if (d < best_d) { best_d = d; best = v; }
        }
        if (best_d > 1e-6) {
            cout << "WARNING mark_solid_block: requested coordinate " << target
                 << " snapped to " << best << " (distance " << best_d << ")\n";
        }
        return best;
    };

    pMesh.obstacle_xmin = snap(pMesh.xvc, x0);
    pMesh.obstacle_xmax = snap(pMesh.xvc, x1);
    pMesh.obstacle_ymin = snap(pMesh.yvc, y0);
    pMesh.obstacle_ymax = snap(pMesh.yvc, y1);

    auto get_idx = [&](int i, int j) { return (j - 1) * N + (i - 1); };

    // 1. Flag solid cells by centroid position.
    int n_solid = 0;
    for (int j = 1; j <= M; ++j) {
        for (int i = 1; i <= N; ++i) {
            InternalCell& cell = pMesh.cells[get_idx(i, j)];
            if (cell.x > pMesh.obstacle_xmin && cell.x < pMesh.obstacle_xmax &&
                cell.y > pMesh.obstacle_ymin && cell.y < pMesh.obstacle_ymax) {
                cell.is_solid = true;
                ++n_solid;
            }
        }
    }

    // 2. Re-derive bc_*/n_* at every fluid/solid interface (both sides).
    for (int j = 1; j <= M; ++j) {
        for (int i = 1; i <= N; ++i) {
            InternalCell& cell = pMesh.cells[get_idx(i, j)];

            if (i < N && pMesh.cells[get_idx(i + 1, j)].is_solid != cell.is_solid) {
                cell.bc_east = true; cell.n_east = -1;
            }
            if (i > 1 && pMesh.cells[get_idx(i - 1, j)].is_solid != cell.is_solid) {
                cell.bc_west = true; cell.n_west = -1;
            }
            if (j < M && pMesh.cells[get_idx(i, j + 1)].is_solid != cell.is_solid) {
                cell.bc_north = true; cell.n_north = -1;
            }
            if (j > 1 && pMesh.cells[get_idx(i, j - 1)].is_solid != cell.is_solid) {
                cell.bc_south = true; cell.n_south = -1;
            }
        }
    }

    cout << "mark_solid_block: obstacle snapped to x=[" << pMesh.obstacle_xmin
         << "," << pMesh.obstacle_xmax << "] y=[" << pMesh.obstacle_ymin << ","
         << pMesh.obstacle_ymax << "], " << n_solid << " solid P-cells marked.\n";
}

void MeshConfig::export_mesh_to_gnuplot(const GridData& data, const std::string& filename) const {
    // Use the grid's own dimensions (critical for staggered U/V meshes)
    int N_cells = data.N_cells_x;
    int M_cells = data.M_cells_y;

    int N = N_cells;
    int M = M_cells;
    
    ofstream outfile(filename);
    if (!outfile.is_open()) {
        cerr << "Error: Could not open file " << filename << std::endl;
        return;
    }

    // --- BLOCK 0: Physical mesh (wireframe) ---
    // Draws the face lines, showing the physical domain.
    outfile << "# Block 0: Mesh Lines (xvc, yvc faces)\n";
    for (int i = 0; i <= N; ++i) { 
        outfile << data.xvc[i] << " " << data.yvc[0] << "\n";
        outfile << data.xvc[i] << " " << data.yvc[M] << "\n";
        outfile << "\n";
    }
    for (int j = 0; j <= M; ++j) { 
        outfile << data.xvc[0] << " " << data.yvc[j] << "\n";
        outfile << data.xvc[N] << " " << data.yvc[j] << "\n";
        outfile << "\n";
    }

    // --- BLOCK 1: All nodes (centres + boundaries) ---
    outfile << "\n\n# Block 1: All Nodes (Internal + Ghost)\n";

    // Loop over i = 0..N+1 (and j = 0..M+1) so the ghost nodes are included

    for (int i = 0; i <= N + 1; ++i) {
        for (int j = 0; j <= M + 1; ++j) {
            outfile << data.x[i] << " " << data.y[j] << "\n";
        }
    }

    outfile.close();
    cout << "Mesh exported successfully to: " << filename << std::endl;
}
 
/**
 * @brief Exports the mesh (grid) to a .vtk file for ParaView.
 * This is the VTK equivalent of 'export_mesh_to_gnuplot'.
 */
void MeshConfig::export_mesh_to_VTK(const GridData& grid, const std::string& filename) const 
{
    // Use the grid's own dimensions (critical for staggered U/V meshes)
    int N_cells = grid.N_cells_x;
    int M_cells = grid.M_cells_y;

    ofstream file(filename);
    if (!file.is_open()) {
        cerr << "Error: Could not open file " << filename << endl;
        return;
    }

    // --- 1. VTK Header ---
    file << "# vtk DataFile Version 2.0\n";
    file << "MeshConfig Generated Grid\n";
    file << "ASCII\n";
    file << "DATASET RECTILINEAR_GRID\n";

    // --- 2. Grid Dimensions (Number of NODES/POINTS) ---
    // A grid with N_cells has (N_cells + 1) nodes (or faces).
    // DIMENSIONS <num_x_points> <num_y_points> <num_z_points>
    file << "DIMENSIONS " << (N_cells + 1) << " " << (M_cells + 1) << " 1\n";

    // --- 3. X-Coordinates (the N_cells+1 face positions) ---
    file << "X_COORDINATES " << (N_cells + 1) << " double\n";
    for (int i = 0; i <= N_cells; ++i) {
        file << grid.xvc[i] << " ";
    }
    file << "\n";

    // --- 4. Y-Coordinates (the M_cells+1 face positions) ---
    file << "Y_COORDINATES " << (M_cells + 1) << " double\n";
    for (int j = 0; j <= M_cells; ++j) {
        file << grid.yvc[j] << " ";
    }
    file << "\n";

    // --- 5. Z-Coordinates (1 point for 2D) ---
    file << "Z_COORDINATES 1 double\n";
    file << "0.0\n";

    // --- 6. Add Cell Data (so ParaView can color the cells) ---
    // This allows you to visualize the cells, not just the wireframe.
    file << "\nCELL_DATA " << (N_cells * M_cells) << "\n";
    file << "SCALARS cell_ID double 1\n";
    file << "LOOKUP_TABLE default\n";
    // Write a unique ID for each cell
    for (int j = 1; j <= M_cells; ++j) {
        for (int i = 1; i <= N_cells; ++i) {
            // Use the internal 1-based indices
            file << (i + (j - 1) * N_cells) << "\n";
        }
    }

    file.close();
    cout << "Mesh exported successfully to: " << filename << endl;
}

/**
 * @brief Reads a configuration file line by line and updates the parameters.
 */
void MeshConfig::get_data(const std::string& filename) {
    ifstream file(filename);
    if (!file.is_open()) {
        cerr << "Error: could not open data file: " << filename << endl;
        return;
    }

    string line, key, value, eq;
    while (getline(file, line)) {
        // Skip comments (#) and empty lines
        if (line.empty() || line[0] == '#') {
            continue;
        }

        stringstream ss(line);
        // Expects the format "key = value"
        if (ss >> key >> eq >> value && eq == "=") {
            // Strip spaces from the key
            key.erase(remove(key.begin(), key.end(), ' '), key.end());
            update_parameter(key, value);
        }
    }
    file.close();
    cout << "Data file '" << filename << "' loaded." << endl;
}

/**
 * @brief Reads the stretching configuration file line by line and updates the parameters.
 * (Same as get_data; both call the same update_parameter helper)
 */
void MeshConfig::get_stretch(const std::string& filename) {
    ifstream file(filename);
    if (!file.is_open()) {
        cerr << "Error: could not open stretching file: " << filename << endl;
        return;
    }

    string line, key, value, eq;
    while (getline(file, line)) {
        if (line.empty() || line[0] == '#') {
            continue;
        }

        stringstream ss(line);
        if (ss >> key >> eq >> value && eq == "=") {
            key.erase(remove(key.begin(), key.end(), ' '), key.end());
            update_parameter(key, value);
        }
    }
    file.close();
    cout << "Stretching file '" << filename << "' loaded." << endl;
}

/**
 * @brief Reads the physical/numerical simulation parameters (Re, U_lid, dt, rho, ...).
 */
void MeshConfig::get_sim_params(const std::string& filename) {
    ifstream file(filename);
    if (!file.is_open()) {
        cerr << "Error: could not open parameter file: " << filename << endl;
        return;
    }

    string line, key, value, eq;
    while (getline(file, line)) {
        if (line.empty() || line[0] == '#') continue;
        stringstream ss(line);
        if (ss >> key >> eq >> value && eq == "=") {
            key.erase(remove(key.begin(), key.end(), ' '), key.end());
            update_parameter(key, value);
        }
    }
    file.close();
    cout << "Simulation parameters '" << filename << "' loaded." << endl;
}

/**
 * @brief Helper that assigns a value to a class member.
 * std::stod = string to double
 * std::stoi = string to integer
 */
void MeshConfig::update_parameter(const std::string& key, const std::string& value) {
    
    // --- "Data" parameters (read by get_data / get_sim_params) ---
    if (key == "L_domain") {
        L_domain = stod(value);
    } else if (key == "H_domain") {
        H_domain = stod(value);
    } else if (key == "W_depth") {
        W_depth = stod(value);
    } else if (key == "Re") {
        Re = stod(value);
    } else if (key == "U_lid") {
        U_lid = stod(value);
    } else if (key == "U_inf") {
        U_inf = stod(value);
    } else if (key == "dt") {
        dt = stod(value);
    } else if (key == "rho") {
        rho = stod(value);
    } else if (key == "max_steps") {
        max_steps = stoi(value);
    } else if (key == "ss_tol") {
        ss_tol = stod(value);
    } else if (key == "solve_energy") {
        solve_energy = stoi(value);
    } else if (key == "Pr") {
        Pr = stod(value);
    }
    // std::pair members are read as separate _x and _y components
    else if (key == "node_NW_x") {
        node_NW.first = stod(value);
    } else if (key == "node_NW_y") {
        node_NW.second = stod(value);
    } else if (key == "node_SE_x") {
        node_SE.first = stod(value);
    } else if (key == "node_SE_y") {
        node_SE.second = stod(value);
    } else if (key == "node_NE_x") {
        node_NE.first = stod(value);
    } else if (key == "node_NE_y") {
        node_NE.second = stod(value);
    } else if (key == "node_SW_x") {
        node_SW.first = stod(value);
    } else if (key == "node_SW_y") {
        node_SW.second = stod(value);
    }

    // --- "Stretching/meshing" parameters (read by get_stretch) ---
    else if (key == "beta_x_west") {
        beta_x_west = stod(value);
    } else if (key == "beta_x_east") {
        beta_x_east = stod(value);
    } else if (key == "beta_y_north") {
        beta_y_north = stod(value);
    } else if (key == "beta_y_south") {
        beta_y_south = stod(value);
    } else if (key == "beta_x_center") {
        beta_x_center = stod(value);
    } else if (key == "beta_y_center") {
        beta_y_center = stod(value);
    }
    
    // --- Cell counts ---
    else if (key == "N_cells_east") {
        N_cells_east = stoi(value);
    } else if (key == "M_cells_east") {
        M_cells_east = stoi(value);
    } else if (key == "N_cells_west") {
        N_cells_west = stoi(value);
    } else if (key == "M_cells_west") {
        M_cells_west = stoi(value);
    } else if (key == "N_cells_center") {
        N_cells_center = stoi(value);
    } else if (key == "M_cells_center") {
        M_cells_center = stoi(value);
    }
    
    // --- Unknown key ---
    else {
        cerr << "Warning: unknown key in configuration file: " << key << endl;
    }
}

// Helper printing a separator line
void print_line() {
    cout << "-------------------------------------------------------------" << endl;
}


// --- MESH INSPECTOR ---
void MeshConfig::inspect_mesh(const GridData& mesh) const {
    
    // 1. PRINT SUMMARY
    cout << "\n=== MESH SUMMARY ===" << endl;
    cout << "Logical dimensions (Nx x My): " << mesh.N_cells_x << " x " << mesh.M_cells_y << endl;
    cout << "Total internal nodes (unknowns): " << mesh.num_active_nodes << endl;
    cout << "Boundary faces:" << endl;
    cout << "  - West:  " << mesh.bound_west.size() << " faces" << endl;
    cout << "  - East:  " << mesh.bound_east.size() << " faces" << endl;
    cout << "  - South: " << mesh.bound_south.size() << " faces" << endl;
    cout << "  - North: " << mesh.bound_north.size() << " faces" << endl;
    print_line();

    // 2. INTERACTIVE LOOP
    char option;
    do {
        cout << "\nSelect what to inspect:" << endl;
        cout << " [I]nternal node (by ID)" << endl;
        cout << " [B]oundaries (list the data of one wall)" << endl;
        cout << " [Q]uit" << endl;
        cout << "Option: ";
        cin >> option;
        option = toupper(option);

        if (option == 'I') {
            int node_id;
            cout << "\nEnter internal node ID (0 to " << mesh.num_active_nodes - 1 << "): ";
            cin >> node_id;

            if (node_id >= 0 && node_id < mesh.num_active_nodes) {
                const InternalCell& cell = mesh.cells[node_id];

                print_line();
                cout << " INTERNAL NODE DETAILS, ID: " << cell.global_id << endl;
                cout << " Centroid (x, y): (" << cell.x << ", " << cell.y << ")" << endl;
                cout << " Volume: " << cell.Vol << " m3" << endl;
                print_line();
                cout << " CONNECTIVITY (NEIGHBOURS):" << endl;

                // Local printing lambda
                auto check_neighbor = [](string dir, int n_id, bool is_bc) {
                    cout << "   " << left << setw(10) << dir << ": ";
                    if (is_bc) {
                        cout << "[BOUNDARY CONDITION / WALL]" << endl;
                    } else if (n_id != -1) {
                        cout << "Internal node ID " << n_id << endl;
                    } else {
                        cout << "Error (not connected)" << endl;
                    }
                };

                check_neighbor("EAST",  cell.n_east,  cell.bc_east);
                check_neighbor("WEST",  cell.n_west,  cell.bc_west);
                check_neighbor("NORTH", cell.n_north, cell.bc_north);
                check_neighbor("SOUTH", cell.n_south, cell.bc_south);
                print_line();

            } else {
                cout << "Error: ID out of range." << endl;
            }

        } else if (option == 'B') {
            char wall;
            cout << "Which wall? (W: West, E: East, S: South, N: North): ";
            cin >> wall;
            wall = toupper(wall);

            const vector<BoundaryFace>* selected_bound = nullptr;
            string wall_name = "";

            if (wall == 'W') { selected_bound = &mesh.bound_west; wall_name = "WEST"; }
            else if (wall == 'E') { selected_bound = &mesh.bound_east; wall_name = "EAST"; }
            else if (wall == 'S') { selected_bound = &mesh.bound_south; wall_name = "SOUTH"; }
            else if (wall == 'N') { selected_bound = &mesh.bound_north; wall_name = "NORTH"; }

            if (selected_bound) {
                cout << "\n--- FACES ON THE " << wall_name << " WALL ---" << endl;
                cout << left << setw(8)  << "Idx"
                             << setw(12) << "x"
                             << setw(12) << "y"
                             << setw(12) << "Area"
                             << setw(20) << "BC type"
                             << setw(10) << "Value" << endl;
                print_line();
                int idx = 0;
                for (const auto& face : *selected_bound) {
                    cout << left << setw(8)  << idx
                                 << setw(12) << face.x
                                 << setw(12) << face.y
                                 << setw(12) << face.Area
                                 << setw(20) << bc_type_to_string(face.type)
                                 << setw(10) << face.value << endl;
                    idx++;
                }
                print_line();
            } else {
                cout << "Invalid wall." << endl;
            }
        }

    } while (option != 'Q');
}


// Parses a BC type string from the input file into the enum
BCType string_to_bc_enum(const string& s) {
    if (s == "INLET") return INLET;
    if (s == "INLET_PARABOLIC") return INLET;   // per-face profile, see define_boundaries
    if (s == "OUTLET") return OUTLET;
    if (s == "WALL") return WALL_ADIABATIC; // Default
    if (s == "WALL_FIXED_VALUE") return WALL_FIXED_VALUE;
    if (s == "WALL_ADIABATIC") return WALL_ADIABATIC;
    if (s == "WALL_ISOTHERMAL") return WALL_ISOTHERMAL;
    if (s == "SYMMETRY") return SYMMETRY;
    return UNDEFINED;
}

void MeshConfig::define_boundaries(GridData& mesh, const std::string& filename) {
    ifstream file(filename);
    if (!file.is_open()) {
        cerr << "Error: could not open " << filename << endl;
        return;
    }

    cout << "--- Loading boundary conditions from " << filename << " ---" << endl;

    string line;
    while (getline(file, line)) {
        // Skip comments (#) and empty lines
        if (line.empty() || line[0] == '#') continue;

        stringstream ss(line);
        string boco_name, face_dir, type_str;
        int idx_start, idx_end;

        // Format: BOCO_NAME FACE START END TYPE [VALUE]
        // Example: INLET WEST 0 50 INLET 1.0
        ss >> boco_name >> face_dir >> idx_start >> idx_end >> type_str;
        double val = 0.0;
        ss >> val; // If the line has no value, val stays 0.0

        // 1. Select the right vector for the direction
        vector<BoundaryFace>* target_vec = nullptr;
        int max_limit = 0;

        if (face_dir == "WEST") { 
            target_vec = &mesh.bound_west; 
            max_limit = mesh.M_cells_y; // West/East depend on Y (j)
        } 
        else if (face_dir == "EAST") { 
            target_vec = &mesh.bound_east; 
            max_limit = mesh.M_cells_y; 
        }
        else if (face_dir == "SOUTH") { 
            target_vec = &mesh.bound_south; 
            max_limit = mesh.N_cells_x; // South/North depend on X (i)
        }
        else if (face_dir == "NORTH") { 
            target_vec = &mesh.bound_north; 
            max_limit = mesh.N_cells_x; 
        }

        // 2. Apply the conditions to the selected range
        if (target_vec) {
            // Range guard
            if (idx_end > max_limit) idx_end = max_limit; 
            if (idx_start < 0) idx_start = 0;

            BCType type_enum = string_to_bc_enum(type_str);

            // INLET_PARABOLIC (WEST/EAST only): fully developed channel profile
            // u(y) = 4 U_max (y - y0)(y1 - y) / H^2, with U_max = val on the axis.
            const bool parabolic = (type_str == "INLET_PARABOLIC");
            if (parabolic && face_dir != "WEST" && face_dir != "EAST") {
                cerr << "Error: INLET_PARABOLIC only supports WEST/EAST" << endl;
            }
            const double y0 = mesh.yvc.front();
            const double y1 = mesh.yvc.back();
            const double H  = y1 - y0;

            for (int k = idx_start; k < idx_end; ++k) {
                // The boundary vectors have an exact size, so check k
                if (k < (int)target_vec->size()) {
                    BoundaryFace& f = (*target_vec)[k];
                    f.type  = type_enum;
                    f.value = parabolic ? 4.0 * val * (f.y - y0) * (y1 - f.y) / (H * H)
                                        : val;
                    f.boco  = boco_name;
                }
            }
            cout << " -> Boco [" << boco_name << "] " << type_str << " on " << face_dir
                 << " [" << idx_start << "-" << idx_end << "]" << endl;
        } else {
            cerr << "Error: unknown direction " << face_dir << endl;
        }
    }
    file.close();
}

// ---------------------------------------------------------------------------
// define_temperature_boundaries
//
// Populates mesh.bound_*_T with temperature BCs, parallel to the velocity BC
// lists. Geometry (positions/areas) is copied from the existing bound_* lists;
// only the BC type/value are overwritten from the file. Defaults to adiabatic
// for any face not mentioned. Same file format as boundaries.txt.
// ---------------------------------------------------------------------------
void MeshConfig::define_temperature_boundaries(GridData& mesh, const std::string& filename) {
    // 1. Initialise T-BC lists as adiabatic copies of the velocity-BC geometry.
    auto init_adiabatic = [](const std::vector<BoundaryFace>& src,
                             std::vector<BoundaryFace>& dst) {
        dst = src;                       // copy face positions / areas / sizes
        for (BoundaryFace& f : dst) {
            f.type  = WALL_ADIABATIC;    // default: zero heat flux
            f.value = 0.0;
            f.boco  = "adiabatic";
        }
    };
    init_adiabatic(mesh.bound_west,  mesh.bound_west_T);
    init_adiabatic(mesh.bound_east,  mesh.bound_east_T);
    init_adiabatic(mesh.bound_south, mesh.bound_south_T);
    init_adiabatic(mesh.bound_north, mesh.bound_north_T);

    ifstream file(filename);
    if (!file.is_open()) {
        cerr << "Error: could not open " << filename
             << " (temperature conditions). Using adiabatic walls by default." << endl;
        return;
    }

    cout << "--- Loading temperature boundary conditions from " << filename << " ---" << endl;

    string line;
    while (getline(file, line)) {
        if (line.empty() || line[0] == '#') continue;

        stringstream ss(line);
        string boco_name, face_dir, type_str;
        int idx_start, idx_end;
        ss >> boco_name >> face_dir >> idx_start >> idx_end >> type_str;
        double val = 0.0;
        ss >> val;

        vector<BoundaryFace>* target_vec = nullptr;
        int max_limit = 0;
        if (face_dir == "WEST")       { target_vec = &mesh.bound_west_T;  max_limit = mesh.M_cells_y; }
        else if (face_dir == "EAST")  { target_vec = &mesh.bound_east_T;  max_limit = mesh.M_cells_y; }
        else if (face_dir == "SOUTH") { target_vec = &mesh.bound_south_T; max_limit = mesh.N_cells_x; }
        else if (face_dir == "NORTH") { target_vec = &mesh.bound_north_T; max_limit = mesh.N_cells_x; }

        if (target_vec) {
            if (idx_end > max_limit) idx_end = max_limit;
            if (idx_start < 0) idx_start = 0;

            BCType type_enum = string_to_bc_enum(type_str);
            for (int k = idx_start; k < idx_end; ++k) {
                if (k < (int)target_vec->size()) {
                    (*target_vec)[k].type  = type_enum;
                    (*target_vec)[k].value = val;
                    (*target_vec)[k].boco  = boco_name;
                }
            }
            cout << " -> T-Boco [" << boco_name << "] " << type_str << " on " << face_dir
                 << " [" << idx_start << "-" << idx_end << "] val=" << val << endl;
        } else {
            cerr << "Error: unknown direction " << face_dir << endl;
        }
    }
    file.close();
}