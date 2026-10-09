#ifndef MESH_CONFIG_H
#define MESH_CONFIG_H

#include <vector>
#include <string>
#include <utility> // Required for std::pair

struct InternalCell {
    int global_id; // El ID en el sistema lineal (0 a N_active-1)
    double x, y;   // Centroide
    double Vol;    // Volumen
    // Vecinos (índices de otras InternalCells)
    int n_east, n_west, n_north, n_south; 
    
    // Banderas para saber si un vecino es frontera
    bool bc_east, bc_west, bc_north, bc_south; 


};
// --- A. NUEVO ENUM PARA TIPOS DE FRONTERA ---
enum BCType {
    WALL_ADIABATIC,
    WALL_FIXED_VALUE,
    INLET,
    OUTLET,
    SYMMETRY,
    WALL_ISOTHERMAL,   // Dirichlet temperature wall (T = value); energy equation only
    UNDEFINED
};

struct BoundaryFace {
    double x, y;         // Posición del centro de la cara
    double Area;         // Área superficial
    int neighbor_id;     // ID de la celda interna que toca esta frontera
    BCType type = UNDEFINED;
    double value = 0.0;  // Valor asociado a la condición (ej. velocidad de entrada)
    std::string boco;    // Nombre del grupo de condición de contorno (boco)
};


struct GridData {
    std::vector<double> xvc, x, yvc, y;
    std::vector<std::vector<double>> V, Se, Sw, Sn, Ss;

    // --- NUEVO: SISTEMA DE CONECTIVIDAD ---
    
    // Matriz de conectividad: [ID_Celda][Direccion]
    // Direcciones: 0=East, 1=West, 2=North, 3=South
    // Si el valor es -1, significa que es una frontera (boundary)
    
    // Matriz de conectividad: [ID_NODO][CARA] -> Devuelve el ID del vecino
    std::vector<std::vector<int>> connectivity;

    // Mapa inverso: Dado un ID, recuperamos sus i, j (útil para debug o post-proceso)
    std::vector<int> id_to_i;
    std::vector<int> id_to_j;
    
    // Total de nodos activos
    int num_active_nodes;

    // 1. Vector de Punteros: Cada nodo apunta a una condición. 
    //    Si es nullptr, es un nodo interno (FLUIDO/DOMINIO).
    // 2. ELEMENTOS INTERNOS (Los que tienen Volumen)
    // Podemos almacenarlos en un vector 1D para el solver lineal
    std::vector<InternalCell> cells; 

    // 3. SUPERFICIES DE FRONTERA (Los que son solo áreas)
    std::vector<BoundaryFace> bound_west;
    std::vector<BoundaryFace> bound_east;
    std::vector<BoundaryFace> bound_south;
    std::vector<BoundaryFace> bound_north;

    // --- Temperature boundary conditions (energy equation only) ---
    // Parallel to the velocity bound_* lists above, populated only when the
    // energy equation is solved (MeshConfig::define_temperature_boundaries).
    // Left empty for the lid-driven cavity, so velocity BCs are untouched.
    std::vector<BoundaryFace> bound_west_T;
    std::vector<BoundaryFace> bound_east_T;
    std::vector<BoundaryFace> bound_south_T;
    std::vector<BoundaryFace> bound_north_T;
    
    // Mapa auxiliar por si necesitas acceso rápido (opcional)
    int N_cells_x, M_cells_y;

    // 2. Vector de IDs (Enteros): Redundante pero necesario para VTK y Paraview
    //    ya que Paraview no lee punteros de C++.
    std::vector<int> node_material_id;
};



class MeshConfig {
public: 
    // --- CONSTRUCTOR ---
    // Aquí declaramos que existe un constructor que se encargará de poner los valores
    MeshConfig(); 

    // --- Grid Dimensions (Parameters) ---
    double L_domain;
    double H_domain;
    double W_depth;

    // --- Simulation Parameters (read from sim_params.txt) ---
    double Re;        // Reynolds number  (Re = U_lid * L / nu)
    double U_lid;     // Lid velocity [m/s]
    double dt;        // Time step [s]
    double rho;       // Fluid density [kg/m³]
    int    max_steps; // Maximum number of time steps
    double ss_tol;    // Steady-state tolerance ||du||_inf / dt

    // --- Energy equation / Boussinesq (Differentially Heated Cavity) ---
    int    solve_energy; // 0 = momentum+pressure only (lid cavity); 1 = solve energy
    double Pr;           // Prandtl number (momentum diffusion = Pr in Ra-Pr scaling)

    // --- Object Positioning on the grid ---
    std::pair<double, double> node_NW;
    std::pair<double, double> node_SE;
    std::pair<double, double> node_NE;
    std::pair<double, double> node_SW;

    // --- Stretching Parameters ---
    double beta_x_west;  
    double beta_x_east;
    double beta_y_north;
    double beta_y_south;
    double beta_x_center;
    double beta_y_center;

    // --- Number of Cells ---
    int N_cells_east;
    int M_cells_east;
    int N_cells_west;
    int M_cells_west;
    int N_cells_center;
    int M_cells_center;
    
    // --- Member Functions (Actions) ---
    void get_data(const std::string& filename);
    void get_stretch(const std::string& filename);
    void get_sim_params(const std::string& filename);

    void define_boundaries(GridData& mesh, const std::string& filename);
    // Populates the parallel bound_*_T temperature BC lists from a file with the
    // same format as boundaries.txt (types WALL_ISOTHERMAL / WALL_ADIABATIC).
    // Geometry is copied from the existing bound_* lists; must be called after
    // generate_mesh(). Only needed when solve_energy == 1.
    void define_temperature_boundaries(GridData& mesh, const std::string& filename);
    GridData generate_mesh() const;
    GridData generate_u_mesh(const GridData& mesh) const;
    GridData generate_v_mesh(const GridData& mesh) const;

    void export_mesh_to_gnuplot(const GridData& data, const std::string& filename) const;
    void export_mesh_to_VTK(const GridData& grid, const std::string& filename) const;
    /**
     * @brief Genera la matriz de conectividad y los mapas de ID global.
     * Debe llamarse DESPUÉS de generate_mesh().
     */
    void build_connectivity(GridData& data) const;

    //void apply_boundary(GridData& data, BoundaryCondition* bc, 
                       // std::function<bool(int i, int j, double x, double y)> selector) const;
    
    // Función auxiliar estática para convertir (i,j) a ID Global
    // Útil si necesitas consultar un ID desde fuera
    static int get_global_id(int i, int j, int N_cells_east_total);

    void inspect_mesh(const GridData& mesh) const;

private:
    static void generate_stretched_coords(std::vector<double>& vc, int num_cells, double length, 
                                          double beta, double offset);
    void update_parameter(const std::string& key, const std::string& value);
};

#endif // MESH_CONFIG_H