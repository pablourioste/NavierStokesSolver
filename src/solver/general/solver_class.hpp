#pragma once
#include <vector>

// Type aliases for the data types
using Matrix = std::vector<std::vector<double>>;
using Vector = std::vector<double>;


// --- BASE CLASS ---
class LinearSolver {
public:
    virtual ~LinearSolver() = default;
    
    // Pure virtual solve function
    virtual Vector solve(const Matrix& A, const Vector& b) const = 0;
};

// --- DERIVED CLASSES ---

class GS : public LinearSolver {
public:
    // Overridden here, implemented in the .cpp
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