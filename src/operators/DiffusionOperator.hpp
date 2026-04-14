#pragma once
#include "mesh/Mesh_config.h"
#include "solver/sparse/solver_class_sparse.hpp"  // Vector

// =============================================================================
// DiffusionOperator
//
// Discretiza el término de difusión  ν∇²u  y  ν∇²v  de las ecuaciones de
// momentum sobre las mallas escalonadas U y V, usando diferencias finitas de
// 2.º orden en malla no uniforme.
//
// Para un nodo U interior (k=2..N, j=1..M):
//
//   d²u/dx² = 2/(dx_e + dx_w) * [(u_E - u_P)/dx_e − (u_P - u_W)/dx_w]
//   d²u/dy² = 2/(dy_n + dy_s) * [(u_N - u_P)/dy_n − (u_P - u_S)/dy_s]
//   resultado = ν * (d²u/dx² + d²u/dy²)
//
// En las paredes (j=1 o j=M) se aplica el método imagen (celda fantasma):
//   u_ghost = 2·u_wall − u_P    →  (u_ghost + u_P)/2 = u_wall  ✓
//   dist    = 2·(y_P − y_wall)   →  distancia correcta a la cara de frontera
//
// Los valores de celda fantasma deben proporcionarse al constructor.
// Para paredes no-slip (velocidad = 0): ghost = −u_P → se pasa u_wall = 0.
// Para la tapa (NORTH, u_wall = U_lid): ghost = 2·U_lid − u_P.
//
// NOTA: este operador solo contiene el término DIFUSIVO. La convección
// permanece en MomentumSolver.
// =============================================================================
class DiffusionOperator {
public:
    // nu              : viscosidad cinemática
    // u_ghost_north/south: valor prescrito en la pared N/S para u
    // v_ghost_west/east  : valor prescrito en la pared O/E para v
    DiffusionOperator(const GridData& pMesh,
                      const GridData& uMesh,
                      const GridData& vMesh,
                      double nu,
                      double u_ghost_north,
                      double u_ghost_south,
                      double v_ghost_west,
                      double v_ghost_east);

    // ν∇²u en todos los nodos U. Solo se calculan los nodos interiores (k=2..N).
    // Los nodos de frontera (k=1, k=N+1) se ponen a 0.0 en el vector resultado.
    void applyU(const Vector& u, Vector& result) const;

    // ν∇²v en todos los nodos V. Solo se calculan los nodos interiores (l=2..M).
    // Los nodos de frontera (l=1, l=M+1) se ponen a 0.0 en el vector resultado.
    void applyV(const Vector& v, Vector& result) const;

    // Núcleos de cálculo por nodo — públicos para que MomentumSolver pueda
    // delegar la parte difusiva de uRHS/vRHS sin duplicar la fórmula.
    double diffuseU(const Vector& u, int k, int j) const;
    double diffuseV(const Vector& v, int i, int l) const;

private:
    const GridData& pMesh_;
    const GridData& uMesh_;
    const GridData& vMesh_;
    double nu_;
    double u_ghost_north_, u_ghost_south_;
    double v_ghost_west_,  v_ghost_east_;

    inline int uIdx(int k, int j) const {
        return (j - 1) * uMesh_.N_cells_x + (k - 1);
    }
    inline int vIdx(int i, int l) const {
        return (l - 1) * vMesh_.N_cells_x + (i - 1);
    }

};
