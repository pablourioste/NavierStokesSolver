#pragma once
#include "mesh/Mesh_config.h"
#include "PressurePoissonSystem.hpp"           // Triplet
#include "solver/sparse/solver_class_sparse.hpp"  // CSRMatrix, Vector

// =============================================================================
// LaplaceOperator
//
// Discretiza el operador de Laplace  ∇²p  sobre la malla de presión P mediante
// el método de volúmenes finitos (FVM).
//
// Para cada celda P con índices lógicos (i,j):
//   a_E = Se[i][j] / (x[i+1] - x[i])      ← contribución cara Este
//   a_W = Sw[i][j] / (x[i]   - x[i-1])    ← contribución cara Oeste
//   a_N = Sn[i][j] / (y[j+1] - y[j])      ← contribución cara Norte
//   a_S = Ss[i][j] / (y[j]   - y[j-1])    ← contribución cara Sur
//   a_P = +(a_E + a_W + a_N + a_S)         ← diagonal (positivo)
//
// Condición de contorno Neumann homogénea (∂p/∂n = 0) por defecto: las caras
// de frontera se ignoran (coeficiente = 0).
//
// La matriz resultante es SPD → apta para CG sin precondicionador.
// Si todas las caras son Neumann el sistema es singular: usar pinReference().
// =============================================================================
class LaplaceOperator {
public:
    // pMesh debe sobrevivir a este objeto.
    explicit LaplaceOperator(const GridData& pMesh);

    // Construye la matriz Laplaciana esparcida A en formato de triplets.
    // Equivale a PressurePoissonSystem::assembleMatrix().
    std::vector<Triplet> assembleToTriplets() const;

    // Construye la matriz en formato CSR directamente.
    CSRMatrix assemble() const;

    // Producto matriz-vector  y = A·p  (sin resolver el sistema).
    // Útil para pruebas MMS: comparar (A·φ)/V con ∇²φ analítico.
    Vector apply(const Vector& p) const;

    // Regularización para sistemas todo-Neumann: fija p[0] = 0.
    // Igual que PressurePoissonSystem::pinReferencePressure().
    void pinReference(std::vector<Triplet>& triplets, Vector& b) const;

private:
    const GridData& pMesh_;

    inline int pIdx(int i, int j) const {
        return (j - 1) * pMesh_.N_cells_x + (i - 1);
    }
};
