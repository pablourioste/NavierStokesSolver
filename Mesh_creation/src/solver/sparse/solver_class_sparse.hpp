#pragma once
#include <vector>
#include "../PressurePoissonSystem.hpp"
// Definimos los alias para los tipos de datos
using Vector = std::vector<double>;

struct CSRMatrix {
    int n;
    std::vector<double> values;
    std::vector<int> col_indices;
    std::vector<int> row_ptr;
    std::vector<double> diagonal;

    // This MUST be inside the struct to be called as CSRMatrix::fromDense
    static CSRMatrix fromDense(const std::vector<std::vector<double>>& dense) {
        CSRMatrix csr;
        csr.n = dense.size();
        csr.row_ptr.push_back(0);
        csr.diagonal.resize(csr.n);

        for (int i = 0; i < csr.n; ++i) {
            for (int j = 0; j < csr.n; ++j) {
                if (dense[i][j] != 0.0) {
                    csr.values.push_back(dense[i][j]);
                    csr.col_indices.push_back(j);
                }
                if (i == j) {
                    csr.diagonal[i] = dense[i][j];
                }
            }
            csr.row_ptr.push_back(static_cast<int>(csr.values.size()));
        }
        return csr;
    }

    // Añadir dentro de struct CSRMatrix
static CSRMatrix fromTriplets(int n, const std::vector<Triplet>& triplets) {
    CSRMatrix csr;
    csr.n = n;
    csr.row_ptr.assign(n + 1, 0);
    csr.diagonal.resize(n, 0.0);

    // 1. Contar elementos por fila
    for (const auto& t : triplets) {
        if (t.value != 0.0) csr.row_ptr[t.row + 1]++;
    }

    // 2. Acumular punteros de fila
    for (int i = 0; i < n; ++i) csr.row_ptr[i + 1] += csr.row_ptr[i];

    // 3. Llenar valores y columnas
    csr.values.resize(csr.row_ptr[n]);
    csr.col_indices.resize(csr.row_ptr[n]);
    std::vector<int> current_pos = csr.row_ptr;

    for (const auto& t : triplets) {
        if (t.value != 0.0) {
            int pos = current_pos[t.row]++;
            csr.values[pos] = t.value;
            csr.col_indices[pos] = t.col;
            if (t.row == t.col) csr.diagonal[t.row] = t.value;
        }
    }
    return csr;
}
};


// ---------------------------------------------------------------------------
// Free helper functions (implemented in solver_sparse.cpp)
// ---------------------------------------------------------------------------
Vector multiply(const CSRMatrix& A, const Vector& v);   // y = A·x
Vector subtract(const Vector& a, const Vector& b);       // a - b
Vector adittion(const Vector& a, const Vector& b);       // a + b
double dot(const Vector& u, const Vector& v);            // u·v

// --- CLASE BASE ---
class LinearSolver {
public:
    virtual ~LinearSolver() = default;
    
    // Declaración de la función virtual pura
    virtual Vector solve(const CSRMatrix& A, const Vector& b) const = 0;
};

// --- CLASES HEREDADAS ---

class GS : public LinearSolver {
public:
    Vector solve(const CSRMatrix& A, const Vector& b) const override;
    // Warm-start: uses x0 as initial guess (same as reference Lid_driven_structured.cpp)
    Vector solve(const CSRMatrix& A, const Vector& b, const Vector& x0) const;
};

class Jacobi : public LinearSolver {
public:
    Vector solve(const CSRMatrix& A, const Vector& b) const override;
    Vector solve(const CSRMatrix& A, const Vector& b, const Vector& x0) const;
};

class SOR : public LinearSolver {
public:
    explicit SOR(double omega = 1.5) : omega_(omega) {}
    Vector solve(const CSRMatrix& A, const Vector& b) const override;
    Vector solve(const CSRMatrix& A, const Vector& b, const Vector& x0) const;
private:
    double omega_;
};

class CG : public LinearSolver {
public:
    Vector solve(const CSRMatrix& A, const Vector& b) const override;
    // Warm-start overload: starts from x0 instead of zero
    Vector solve(const CSRMatrix& A, const Vector& b, const Vector& x0) const;
};

class PCG : public LinearSolver {
public:
    Vector solve(const CSRMatrix& A, const Vector& b) const override;
    // Warm-start + Jacobi preconditioning: M = diag(A), faster for non-uniform grids
    Vector solve(const CSRMatrix& A, const Vector& b, const Vector& x0) const;
};