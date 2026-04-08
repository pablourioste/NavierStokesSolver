#include <iostream>
#include <fstream>
#include <string>
#include <sstream> // (La necesitamos para 'stod' y 'stoi' en update_parameter, pero no para leer)
#include <algorithm> // Para 'remove' y 'find_if'
#include <cctype>    // Para 'isspace'
#include <vector>
#include <cmath>    // Para funciones matemáticas como sinh, abs
#include <iomanip> // Para formateo bonito

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
    // 1. Inicializamos Dimensiones
    L_domain(3.0),
    H_domain(3.0),
    W_depth(0.1),
    Re(100.0),
    U_lid(1.0),
    dt(1e-3),
    rho(1.0),
    max_steps(50000),
    ss_tol(1e-6),

    // 2. Inicializamos Posiciones (std::pair)
    // NW = North-West: min x, max y | NE = North-East: max x, max y
    // SW = South-West: min x, min y | SE = South-East: max x, min y
    node_NW(1.5, 2.0),
    node_SE(2.0, 1.5),
    node_NE(2.0, 2.0),
    node_SW(1.5, 1.5),

    // 3. Inicializamos Stretching
    beta_x_west(2.0),
    beta_x_east(-2.0),
    beta_y_north(-2.0),
    beta_y_south(2.0),
    beta_x_center(0.0),
    beta_y_center(0.0),

    // 4. Inicializamos Número de Celdas
    N_cells_east(50),
    M_cells_east(50),
    N_cells_west(50),
    M_cells_west(50),
    N_cells_center(50),
    M_cells_center(50)
    {
    // El cuerpo queda vacío porque solo estamos inicializando valores
    cout << "MeshConfig object created with default parameters." << endl;
}
// ==========================================
// FIN DEL CONSTRUCTOR
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

    // --- 1. GENERAR FRONTERAS (SUPERFICIES) ---
    cout << "Generando Fronteras..." << endl;

    // Frontera Oeste (cara i=0)
    for (int j = 1; j <= M_cells; ++j) {
        BoundaryFace face;
        face.x = data.xvc[0]; 
        face.y = (data.yvc[j] + data.yvc[j-1]) / 2.0; 
        face.Area = (data.yvc[j] - data.yvc[j-1]) * W_depth;
        // Conecta con la primera columna interna (el ID se gestiona luego)
        data.bound_west.push_back(face);
    }

    // Frontera Este (cara i=N_cells)
    for (int j = 1; j <= M_cells; ++j) {
        BoundaryFace face;
        face.x = data.xvc[N_cells]; 
        face.y = (data.yvc[j] + data.yvc[j-1]) / 2.0; 
        face.Area = (data.yvc[j] - data.yvc[j-1]) * W_depth;
        data.bound_east.push_back(face);
    }

    // Frontera Sur (cara j=0)
    for (int i = 1; i <= N_cells; ++i) {
        BoundaryFace face;
        face.x = (data.xvc[i] + data.xvc[i-1]) / 2.0;
        face.y = data.yvc[0];
        face.Area = (data.xvc[i] - data.xvc[i-1]) * W_depth;
        data.bound_south.push_back(face);
    }

    // Frontera Norte (cara j=M_cells)
    for (int i = 1; i <= N_cells; ++i) {
        BoundaryFace face;
        face.x = (data.xvc[i] + data.xvc[i-1]) / 2.0;
        face.y = data.yvc[M_cells];
        face.Area = (data.xvc[i] - data.xvc[i-1]) * W_depth;
        data.bound_north.push_back(face);
    }

    // --- 2. GENERAR ELEMENTOS INTERNOS (VOLÚMENES) ---
    cout << "Generando Elementos de Volumen Internos..." << endl;
    
    int current_id = 0;
    
    // Nota: Iteramos SOLO sobre las celdas de fluido internas
    for (int j = 1; j <= M_cells; ++j) {
        for (int i = 1; i <= N_cells; ++i) {
            InternalCell cell;
            
            cell.global_id = current_id++; // 0, 1, 2...
            
            // Centroides
            cell.x = (data.xvc[i] + data.xvc[i-1]) / 2.0;
            cell.y = (data.yvc[j] + data.yvc[j-1]) / 2.0;
            
            // Propiedades geométricas
            double dx = data.xvc[i] - data.xvc[i-1];
            double dy = data.yvc[j] - data.yvc[j-1];
            cell.Vol = dx * dy * W_depth;
            
            // Agregamos a la lista
            data.cells.push_back(cell);
        }
    }
    
    data.num_active_nodes = current_id;
    return data;
}

GridData MeshConfig::generate_u_mesh(const GridData& pMesh) const {
    GridData uMesh;

    // 1. Configuración Básica
    // La malla U tiene la misma Y que P, pero diferente X
    uMesh.M_cells_y = pMesh.M_cells_y;
    uMesh.yvc = pMesh.yvc; // Copiamos las caras Y exactas de P

    // La malla U tiene N + 1 puntos de velocidad (incluyendo fronteras)
    // Por tanto, tendrá N + 1 celdas de control en X
    uMesh.N_cells_x = pMesh.N_cells_x + 1;

    // 2. Generar nuevas coordenadas X (xvc para U)
    // Los límites de las celdas U son los CENTROS de las celdas P
    uMesh.xvc.push_back(pMesh.xvc[0]); // Límite izquierdo (pared)

    for (int i = 1; i <= pMesh.N_cells_x; ++i) {
        // Promedio entre caras de P para obtener el centro de P
        double p_center = (pMesh.xvc[i] + pMesh.xvc[i-1]) / 2.0;
        uMesh.xvc.push_back(p_center);
    }
    
    uMesh.xvc.push_back(pMesh.xvc[pMesh.N_cells_x]); // Límite derecho (pared)

    // 3. Generar Elementos (Celdas) para U
    // Nota: Esto es crucial para que export_to_VTK funcione
    int current_id = 0;
    
    // Iteramos sobre las nuevas dimensiones
    // j va hasta M (igual que P), i va hasta N+1 (porque hay una U más)
    for (int j = 1; j <= uMesh.M_cells_y; ++j) {
        for (int i = 1; i <= uMesh.N_cells_x; ++i) {
            InternalCell cell;
            cell.global_id = current_id++;

            // Límites de la celda U actual
            double x_left = uMesh.xvc[i-1];
            double x_right = uMesh.xvc[i];
            double y_bottom = uMesh.yvc[j-1];
            double y_top = uMesh.yvc[j];

            // Centroide de la celda U (Donde vive la velocidad u_i,j)
            cell.x = (x_left + x_right) / 2.0;
            cell.y = (y_bottom + y_top) / 2.0;

            // Volumen
            double dx = x_right - x_left;
            double dy = y_top - y_bottom;
            cell.Vol = dx * dy * W_depth; // Asumiendo que W_depth es accesible

            uMesh.cells.push_back(cell);
        }
    }
    
    uMesh.num_active_nodes = current_id;
    return uMesh;
}

GridData MeshConfig::generate_v_mesh(const GridData& pMesh) const {
    GridData vMesh;

    // 1. Configuración Básica
    // La malla V tiene la misma X que P, pero diferente Y
    vMesh.N_cells_x = pMesh.N_cells_x;
    vMesh.xvc = pMesh.xvc; // Copiamos las caras X exactas de P

    // La malla V tiene M + 1 puntos de velocidad
    vMesh.M_cells_y = pMesh.M_cells_y + 1;

    // 2. Generar nuevas coordenadas Y (yvc para V)
    // Los límites de las celdas V son los CENTROS de las celdas P
    vMesh.yvc.push_back(pMesh.yvc[0]); // Pared Sur

    for (int j = 1; j <= pMesh.M_cells_y; ++j) {
        double p_center = (pMesh.yvc[j] + pMesh.yvc[j-1]) / 2.0;
        vMesh.yvc.push_back(p_center);
    }

    vMesh.yvc.push_back(pMesh.yvc[pMesh.M_cells_y]); // Pared Norte

    // 3. Generar Elementos (Celdas) para V
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

// Enumeración para facilitar lectura
enum Direction { EAST = 0, WEST = 1, NORTH = 2, SOUTH = 3 };

int MeshConfig::get_global_id(int i, int j, int N_total) {
    return (j - 1) * N_total + (i - 1);
}

void MeshConfig::build_connectivity(GridData& data) const {

    int N = data.N_cells_x;
    int M = data.M_cells_y;
    int total = N * M;

    // Inicializar los mapas de conectividad e índices inversos
    data.id_to_i.resize(total);
    data.id_to_j.resize(total);
    data.connectivity.resize(total, std::vector<int>(4, -1));

    // Lambda auxiliar para obtener índice 1D desde 2D (i,j)
    // Asumiendo orden lexicográfico (bucle j, luego i) como en generate_mesh
    auto get_idx = [&](int i, int j) {
        return (j - 1) * N + (i - 1);
    };

    for (int j = 1; j <= M; ++j) {
        for (int i = 1; i <= N; ++i) {

            int current_idx = get_idx(i, j);
            InternalCell& cell = data.cells[current_idx];

            // Mapas inversos
            data.id_to_i[current_idx] = i;
            data.id_to_j[current_idx] = j;

            // --- VECINO ESTE ---
            if (i < N) {
                cell.n_east = get_idx(i + 1, j);
                cell.bc_east = false;
            } else {
                cell.n_east = -1; // No hay vecino interno
                cell.bc_east = true; // Toca la Frontera Este
            }
            data.connectivity[current_idx][EAST] = cell.n_east;

            // --- VECINO OESTE ---
            if (i > 1) {
                cell.n_west = get_idx(i - 1, j);
                cell.bc_west = false;
            } else {
                cell.n_west = -1;
                cell.bc_west = true; // Toca la Frontera Oeste
            }
            data.connectivity[current_idx][WEST] = cell.n_west;

            // --- VECINO NORTE ---
            if (j < M) {
                cell.n_north = get_idx(i, j + 1);
                cell.bc_north = false;
            } else {
                cell.n_north = -1;
                cell.bc_north = true; // Toca la Frontera Norte
            }
            data.connectivity[current_idx][NORTH] = cell.n_north;

            // --- VECINO SUR ---
            if (j > 1) {
                cell.n_south = get_idx(i, j - 1);
                cell.bc_south = false;
            } else {
                cell.n_south = -1;
                cell.bc_south = true; // Toca la Frontera Sur
            }
            data.connectivity[current_idx][SOUTH] = cell.n_south;
        }
    }

    cout << "Conectividad construida. Nodos Internos Activos: " << data.num_active_nodes << endl;
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

    // --- BLOCK 0: Malla física (Wireframe) ---
    // Esto dibuja las líneas de las caras. Lo dejamos igual para ver el dominio físico.
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

    // --- BLOCK 1: Todos los Nodos (Centros + Fronteras) ---
    outfile << "\n\n# Block 1: All Nodes (Internal + Ghost)\n";
    
    // CORRECCIÓN AQUÍ:
    // Antes: for (int i = 1; i <= N; ...
    // Ahora: for (int i = 0; i <= N + 1; ... (Incluye 0 y N+1)
    
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
 * @brief Lee un archivo de configuración línea por línea y actualiza los parámetros.
 */
void MeshConfig::get_data(const std::string& filename) {
    ifstream file(filename);
    if (!file.is_open()) {
        cerr << "Error: No se pudo abrir el archivo de datos: " << filename << endl;
        return;
    }

    string line, key, value, eq;
    while (getline(file, line)) {
        // Ignorar comentarios (#) y líneas vacías
        if (line.empty() || line[0] == '#') {
            continue;
        }

        stringstream ss(line);
        // Espera el formato "clave = valor"
        if (ss >> key >> eq >> value && eq == "=") {
            // Limpiar espacios de la clave
            key.erase(remove(key.begin(), key.end(), ' '), key.end());
            update_parameter(key, value);
        }
    }
    file.close();
    cout << "Archivo de datos '" << filename << "' cargado." << endl;
}

/**
 * @brief Lee un archivo de configuración línea por línea y actualiza los parámetros.
 * (Esta función es idéntica a get_data, pero llama al mismo helper)
 */
void MeshConfig::get_stretch(const std::string& filename) {
    ifstream file(filename);
    if (!file.is_open()) {
        cerr << "Error: No se pudo abrir el archivo de estiramiento: " << filename << endl;
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
    cout << "Archivo de estiramiento '" << filename << "' cargado." << endl;
}

/**
 * @brief Lee parámetros físicos/numéricos de la simulación (Re, U_lid, dt, rho, ...).
 */
void MeshConfig::get_sim_params(const std::string& filename) {
    ifstream file(filename);
    if (!file.is_open()) {
        cerr << "Error: No se pudo abrir el archivo de parámetros: " << filename << endl;
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
    cout << "Parámetros de simulación '" << filename << "' cargados." << endl;
}

/**
 * @brief Función 'helper' que asigna un valor a un miembro de la clase.
 * std::stod = string to double
 * std::stoi = string to integer
 */
void MeshConfig::update_parameter(const std::string& key, const std::string& value) {
    
    // --- Parámetros de "Datos" (leídos por get_data) ---
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
    } else if (key == "dt") {
        dt = stod(value);
    } else if (key == "rho") {
        rho = stod(value);
    } else if (key == "max_steps") {
        max_steps = stoi(value);
    } else if (key == "ss_tol") {
        ss_tol = stod(value);
    }
    // Para std::pair, es más fácil leerlos como componentes _x y _y
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

    // --- Parámetros de "Estiramiento/Mallado" (leídos por get_stretch) ---
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
    
    // --- Cuentas de Celdas ---
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
    
    // --- Manejar clave desconocida ---
    else {
        cerr << "Advertencia: Clave desconocida en el archivo de configuración: " << key << endl;
    }
}

// Función auxiliar para imprimir líneas separadoras
void print_line() {
    cout << "-------------------------------------------------------------" << endl;
}


// --- IMPLEMENTACIÓN DEL INSPECTOR ---
void MeshConfig::inspect_mesh(const GridData& mesh) const {
    
    // 1. IMPRIMIR RESUMEN
    cout << "\n=== RESUMEN DE LA MALLA ===" << endl;
    cout << "Dimensiones logicas (Nx x My): " << mesh.N_cells_x << " x " << mesh.M_cells_y << endl;
    cout << "Total Nodos Internos (Incognitas): " << mesh.num_active_nodes << endl;
    cout << "Caras de Frontera:" << endl;
    cout << "  - Oeste (West):  " << mesh.bound_west.size() << " caras" << endl;
    cout << "  - Este (East):   " << mesh.bound_east.size() << " caras" << endl;
    cout << "  - Sur (South):   " << mesh.bound_south.size() << " caras" << endl;
    cout << "  - Norte (North): " << mesh.bound_north.size() << " caras" << endl;
    print_line();

    // 2. BUCLE INTERACTIVO
    char opcion;
    do {
        cout << "\nSelecciona que deseas inspeccionar:" << endl;
        cout << " [I]nspeccionar un Nodo Interno (Por ID)" << endl;
        cout << " [F]ronteras (Listar datos de una pared)" << endl;
        cout << " [S]alir" << endl;
        cout << "Opcion: ";
        cin >> opcion;
        opcion = toupper(opcion);

        if (opcion == 'I') {
            int node_id;
            cout << "\nIntroduce ID del nodo interno (0 a " << mesh.num_active_nodes - 1 << "): ";
            cin >> node_id;

            if (node_id >= 0 && node_id < mesh.num_active_nodes) {
                const InternalCell& cell = mesh.cells[node_id];

                print_line();
                cout << " DETALLES DEL NODO INTERNO ID: " << cell.global_id << endl;
                cout << " Centroide (x, y): (" << cell.x << ", " << cell.y << ")" << endl;
                cout << " Volumen: " << cell.Vol << " m3" << endl;
                print_line();
                cout << " CONECTIVIDAD (VECINOS):" << endl;

                // Lambda local para imprimir
                auto check_neighbor = [](string dir, int n_id, bool is_bc) {
                    cout << "   " << left << setw(10) << dir << ": ";
                    if (is_bc) {
                        cout << "[CONDICION DE FRONTERA / PARED]" << endl;
                    } else if (n_id != -1) {
                        cout << "Nodo Interno ID " << n_id << endl;
                    } else {
                        cout << "Error (No conectado)" << endl;
                    }
                };

                check_neighbor("ESTE",  cell.n_east,  cell.bc_east);
                check_neighbor("OESTE", cell.n_west,  cell.bc_west);
                check_neighbor("NORTE", cell.n_north, cell.bc_north);
                check_neighbor("SUR",   cell.n_south, cell.bc_south);
                print_line();

            } else {
                cout << "Error: ID fuera de rango." << endl;
            }

        } else if (opcion == 'F') {
            char pared;
            cout << "Cual pared? (W: Oeste, E: Este, S: Sur, N: Norte): ";
            cin >> pared;
            pared = toupper(pared);

            const vector<BoundaryFace>* selected_bound = nullptr;
            string nombre_pared = "";

            if (pared == 'W') { selected_bound = &mesh.bound_west; nombre_pared = "OESTE"; }
            else if (pared == 'E') { selected_bound = &mesh.bound_east; nombre_pared = "ESTE"; }
            else if (pared == 'S') { selected_bound = &mesh.bound_south; nombre_pared = "SUR"; }
            else if (pared == 'N') { selected_bound = &mesh.bound_north; nombre_pared = "NORTE"; }

            if (selected_bound) {
                cout << "\n--- LISTADO DE CARAS EN PARED " << nombre_pared << " ---" << endl;
                cout << left << setw(8)  << "Idx"
                             << setw(12) << "x"
                             << setw(12) << "y"
                             << setw(12) << "Area"
                             << setw(20) << "Tipo BC"
                             << setw(10) << "Valor" << endl;
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
                cout << "Pared no valida." << endl;
            }
        }

    } while (opcion != 'S');
}


// Función para parsear el string del archivo a Enum
BCType string_to_bc_enum(const string& s) {
    if (s == "INLET") return INLET;
    if (s == "OUTLET") return OUTLET;
    if (s == "WALL") return WALL_ADIABATIC; // Por defecto
    if (s == "WALL_FIXED_VALUE") return WALL_FIXED_VALUE;
    if (s == "SYMMETRY") return SYMMETRY;
    return UNDEFINED;
}

void MeshConfig::define_boundaries(GridData& mesh, const std::string& filename) {
    ifstream file(filename);
    if (!file.is_open()) {
        cerr << "Error: No se pudo abrir " << filename << endl;
        return;
    }

    cout << "--- Cargando Condiciones de Frontera desde " << filename << " ---" << endl;

    string line;
    while (getline(file, line)) {
        // Ignorar comentarios (#) o líneas vacías
        if (line.empty() || line[0] == '#') continue;

        stringstream ss(line);
        string boco_name, face_dir, type_str;
        int idx_start, idx_end;

        // Leer formato: BOCO_NAME CARA INICIO FIN TIPO [VALOR]
        // Ejemplo: INLET WEST 0 50 INLET 1.0
        ss >> boco_name >> face_dir >> idx_start >> idx_end >> type_str;
        double val = 0.0;
        ss >> val; // Si no hay valor en la línea, val queda en 0.0

        // 1. Seleccionar el vector correcto según la dirección
        vector<BoundaryFace>* target_vec = nullptr;
        int max_limit = 0;

        if (face_dir == "WEST") { 
            target_vec = &mesh.bound_west; 
            max_limit = mesh.M_cells_y; // West/East dependen de Y (j)
        } 
        else if (face_dir == "EAST") { 
            target_vec = &mesh.bound_east; 
            max_limit = mesh.M_cells_y; 
        }
        else if (face_dir == "SOUTH") { 
            target_vec = &mesh.bound_south; 
            max_limit = mesh.N_cells_x; // South/North dependen de X (i)
        }
        else if (face_dir == "NORTH") { 
            target_vec = &mesh.bound_north; 
            max_limit = mesh.N_cells_x; 
        }

        // 2. Aplicar las condiciones al rango seleccionado
        if (target_vec) {
            // Protección de rango
            if (idx_end > max_limit) idx_end = max_limit; 
            if (idx_start < 0) idx_start = 0;

            BCType type_enum = string_to_bc_enum(type_str);

            for (int k = idx_start; k < idx_end; ++k) {
                // Como los vectores de frontera tienen tamaño exacto, verificamos k
                if (k < (int)target_vec->size()) {
                    (*target_vec)[k].type  = type_enum;
                    (*target_vec)[k].value = val;
                    (*target_vec)[k].boco  = boco_name;
                }
            }
            cout << " -> Boco [" << boco_name << "] " << type_str << " en " << face_dir
                 << " [" << idx_start << "-" << idx_end << "]" << endl;
        } else {
            cerr << "Error: Direccion desconocida " << face_dir << endl;
        }
    }
    file.close();
}