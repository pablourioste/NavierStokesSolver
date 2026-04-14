#pragma once
#include "mesh/Mesh_config.h"
#include "solver/sparse/solver_class_sparse.hpp"  // Vector

// =============================================================================
// DivergenceOperator
//
// Discretiza la divergencia  ∇·u  integrando el flujo neto sobre cada celda P
// mediante el método de volúmenes finitos (FVM).
//
// Para la celda P con índices lógicos (i,j):
//
//   flux = u_E * Se  −  u_W * Sw  +  v_N * Sn  −  v_S * Ss
//
//   Convención de normal exterior:
//     Cara Este  (n̂=+x): +u_E * Se
//     Cara Oeste (n̂=−x): −u_W * Sw
//     Cara Norte (n̂=+y): +v_N * Sn
//     Cara Sur   (n̂=−y): −v_S * Ss
//
//   Mapeado de índices escalonados para P(i,j):
//     u_E = u[uIdx(i+1, j)],  u_W = u[uIdx(i,   j)]
//     v_N = v[vIdx(i,   j+1)], v_S = v[vIdx(i,   j)]
//
// El resultado se divide por el volumen de celda V[i][j], dando
// la divergencia por unidad de volumen (dimensiones de [1/tiempo] si
// u tiene dimensiones de velocidad). Esto permite comparar directamente
// con ∇·u analítico en las pruebas MMS.
//
// IMPORTANTE: sin factores ρ/dt — operador puramente geométrico.
// =============================================================================
class DivergenceOperator {
public:
    DivergenceOperator(const GridData& pMesh,
                       const GridData& uMesh,
                       const GridData& vMesh);

    // Calcula ∇·u discreto por unidad de volumen en cada celda P.
    // Entrada: u (tamaño uMesh.num_active_nodes), v (tamaño vMesh.num_active_nodes)
    // Salida:  vector de tamaño pMesh.num_active_nodes
    Vector apply(const Vector& u, const Vector& v) const;

private:
    const GridData& pMesh_;
    const GridData& uMesh_;
    const GridData& vMesh_;

    inline int pIdx(int i, int j) const {
        return (j - 1) * pMesh_.N_cells_x + (i - 1);
    }
    inline int uIdx(int k, int j) const {
        return (j - 1) * uMesh_.N_cells_x + (k - 1);
    }
    inline int vIdx(int i, int l) const {
        return (l - 1) * vMesh_.N_cells_x + (i - 1);
    }
};
