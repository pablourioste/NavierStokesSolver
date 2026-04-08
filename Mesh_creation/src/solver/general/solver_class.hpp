#pragma once
#include <vector>

// Definimos los alias para los tipos de datos
using Matrix = std::vector<std::vector<double>>;
using Vector = std::vector<double>;


// --- CLASE BASE ---
class LinearSolver {
public:
    virtual ~LinearSolver() = default;
    
    // Declaración de la función virtual pura
    virtual Vector solve(const Matrix& A, const Vector& b) const = 0;
};

// --- CLASES HEREDADAS ---

class GS : public LinearSolver {
public:
    // Solo declaramos que sobrescribiremos esta función, la implementaremos en el .cpp
    Vector solve(const Matrix& A, const Vector& b) const override;
};

class Jacobi : public LinearSolver {
public:
    Vector solve(const Matrix& A, const Vector& b) const override;
};

class SOR : public LinearSolver {
public:
    Vector solve(const Matrix& A, const Vector& b) const override;
private:
    double getOptimalOmega(const Matrix& A) const;
};

class CG : public LinearSolver {
public:
    Vector solve(const Matrix& A, const Vector& b) const override;
};

class PCG : public LinearSolver {
public:
    Vector solve(const Matrix& A, const Vector& b) const override;
};