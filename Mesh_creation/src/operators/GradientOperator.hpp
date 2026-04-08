#pragma once
#include "mesh/Mesh_config.h"
#include "solver/sparse/solver_class_sparse.hpp"  // Vector

// =============================================================================
// GradientOperator
//
// Discretiza el gradiente de presión  ∇p  de la malla P a las mallas U/V,
// usando diferencias finitas centradas de 2.º orden sobre la malla escalonada.
//
// Para un nodo U interior (k=2..N, j=1..M):
//   (∂p/∂x)|_{k,j} = (p[pIdx(k,j)] − p[pIdx(k-1,j)]) / (x_P[k] − x_P[k-1])
//   donde x_P[k] es el centro de la celda P con índice lógico k.
//
// Para un nodo V interior (i=1..N, l=2..M):
//   (∂p/∂y)|_{i,l} = (p[pIdx(i,l)] − p[pIdx(i,l-1)]) / (y_P[l] − y_P[l-1])
//
// IMPORTANTE: este operador es puramente geométrico, sin factores dt/ρ.
// La corrección de velocidad en PressurePoissonSystem::projectVelocity()
// aplica adicionalmente el factor −(dt/ρ).
//
// Los nodos de frontera (k=1, k=N+1, l=1, l=M+1) se ponen a 0.0 en la salida.
// =============================================================================
class GradientOperator {
public:
    GradientOperator(const GridData& pMesh,
                     const GridData& uMesh,
                     const GridData& vMesh);

    // Calcula el gradiente discreto de p.
    // Entrada:   p      — campo de presión en nodos P, tamaño pMesh.num_active_nodes
    // Salidas:   grad_u — ∂p/∂x en nodos U interiores, tamaño uMesh.num_active_nodes
    //            grad_v — ∂p/∂y en nodos V interiores, tamaño vMesh.num_active_nodes
    void apply(const Vector& p, Vector& grad_u, Vector& grad_v) const;

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
