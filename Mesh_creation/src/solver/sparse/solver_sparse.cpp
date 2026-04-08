#include "solver_class_sparse.hpp"
#include <iostream>
#include <cmath>
#include <stdexcept>
#include <vector>
#include <algorithm>
#include <numeric>
using namespace std;

// ---------------------------------------------------------------------------
// Jacobi — warm-start core (all variants delegate here)
// ---------------------------------------------------------------------------
static Vector jacobi_core(const CSRMatrix& A, const Vector& b, const Vector& x0) {
    int n = b.size();
    Vector x = x0;
    Vector x_new(n);

    const double b_norm = sqrt(dot(b, b));
    const double tol_abs = (b_norm > 0.0) ? 1e-6 * b_norm : 1e-10;
    const int max_iter = std::max(500, (int)(10.0 * sqrt((double)n)));

    if (b_norm == 0.0) return x;

    for (int iter = 1; iter <= max_iter; ++iter) {
        for (int i = 0; i < n; ++i) {
            double sum = 0.0;
            for (int k = A.row_ptr[i]; k < A.row_ptr[i+1]; ++k) {
                int j = A.col_indices[k];
                if (j != i) sum += A.values[k] * x[j];
            }
            x_new[i] = (b[i] - sum) / A.diagonal[i];
        }
        // Convergence check on ||r||
        double r_norm = 0.0;
        for (int i = 0; i < n; ++i) {
            double Ax_i = 0.0;
            for (int k = A.row_ptr[i]; k < A.row_ptr[i+1]; ++k)
                Ax_i += A.values[k] * x_new[A.col_indices[k]];
            double ri = b[i] - Ax_i;
            r_norm += ri * ri;
        }
        x = x_new;
        if (sqrt(r_norm) < tol_abs) return x;
    }
    return x;
}

Vector Jacobi::solve(const CSRMatrix& A, const Vector& b) const {
    return jacobi_core(A, b, Vector(b.size(), 0.0));
}
Vector Jacobi::solve(const CSRMatrix& A, const Vector& b, const Vector& x0) const {
    return jacobi_core(A, b, x0);
}


// ---------------------------------------------------------------------------
// Gauss-Seidel — warm-start core (same structure as Lid_driven_structured.cpp)
//
// The reference code (Lid_driven_structured.cpp) uses the previous time step's
// pressure as initial guess (P is modified in-place, acting as warm start).
// Our cold-start version (x0=0) needs O(N²) iterations for Poisson; with warm
// start from the previous step, only O(1)–O(N) sweeps suffice in steady regime.
// ---------------------------------------------------------------------------
static Vector gs_core(const CSRMatrix& A, const Vector& b, const Vector& x0) {
    int n = b.size();
    Vector x = x0;

    const double b_norm = sqrt(dot(b, b));
    const double tol_abs = (b_norm > 0.0) ? 1e-6 * b_norm : 1e-10;
    // Without warm start GS needs O(N²) sweeps; cap high enough to converge
    const int max_iter = std::max(500, (int)(10.0 * sqrt((double)n)));

    if (b_norm == 0.0) return x;

    for (int iter = 1; iter <= max_iter; ++iter) {
        // One GS sweep — uses updated values immediately (in-place, like reference)
        for (int i = 0; i < n; ++i) {
            double sum = 0.0;
            for (int k = A.row_ptr[i]; k < A.row_ptr[i+1]; ++k) {
                int j = A.col_indices[k];
                if (j != i) sum += A.values[k] * x[j];
            }
            x[i] = (b[i] - sum) / A.diagonal[i];
        }

        // Convergence: ||b - Ax|| < tol_abs
        double r_norm = 0.0;
        for (int i = 0; i < n; ++i) {
            double Ax_i = 0.0;
            for (int k = A.row_ptr[i]; k < A.row_ptr[i+1]; ++k)
                Ax_i += A.values[k] * x[A.col_indices[k]];
            double ri = b[i] - Ax_i;
            r_norm += ri * ri;
        }
        if (sqrt(r_norm) < tol_abs) return x;
    }
    cout << "GS WARNING: did not converge in " << max_iter << " iters.\n";
    return x;
}

Vector GS::solve(const CSRMatrix& A, const Vector& b) const {
    return gs_core(A, b, Vector(b.size(), 0.0));
}
Vector GS::solve(const CSRMatrix& A, const Vector& b, const Vector& x0) const {
    return gs_core(A, b, x0);
}

// ---------------------------------------------------------------------------
// SOR — warm-start core
// ω ∈ (1,2): over-relaxation. For Poisson on N×N grid the optimal value is
//   ω_opt = 2 / (1 + sin(π/(N+1)))  ≈ 2 - 2π/N for large N.
// The constructor default ω=1.5 is a safe, generally good choice.
// Without warm start SOR is still O(N²) iters; with warm start ~O(N).
// ---------------------------------------------------------------------------
static Vector sor_core(const CSRMatrix& A, const Vector& b,
                       const Vector& x0, double omega) {
    int n = b.size();
    Vector x = x0;

    const double b_norm = sqrt(dot(b, b));
    const double tol_abs = (b_norm > 0.0) ? 1e-6 * b_norm : 1e-10;
    const int max_iter = std::max(500, (int)(10.0 * sqrt((double)n)));

    if (b_norm == 0.0) return x;

    for (int iter = 1; iter <= max_iter; ++iter) {
        for (int i = 0; i < n; ++i) {
            double L = 0.0, U = 0.0;
            for (int k = A.row_ptr[i]; k < A.row_ptr[i+1]; ++k) {
                int j = A.col_indices[k];
                if (j < i) L += A.values[k] * x[j];       // already updated
                else if (j > i) U += A.values[k] * x[j];  // not yet updated
            }
            double gs_val = (b[i] - L - U) / A.diagonal[i];
            x[i] = (1.0 - omega) * x[i] + omega * gs_val;
        }

        // Convergence: ||b - Ax||
        double r_norm = 0.0;
        for (int i = 0; i < n; ++i) {
            double Ax_i = 0.0;
            for (int k = A.row_ptr[i]; k < A.row_ptr[i+1]; ++k)
                Ax_i += A.values[k] * x[A.col_indices[k]];
            double ri = b[i] - Ax_i;
            r_norm += ri * ri;
        }
        if (sqrt(r_norm) < tol_abs) return x;
    }
    cout << "SOR WARNING: did not converge in " << max_iter << " iters.\n";
    return x;
}

Vector SOR::solve(const CSRMatrix& A, const Vector& b) const {
    return sor_core(A, b, Vector(b.size(), 0.0), omega_);
}
Vector SOR::solve(const CSRMatrix& A, const Vector& b, const Vector& x0) const {
    return sor_core(A, b, x0, omega_);
}

Vector multiply(const CSRMatrix& A, const Vector& v) {
    Vector res(A.n, 0.0);
    for (int i = 0; i < A.n; i++) {
        // row_ptr[i] to row_ptr[i+1] gives the range of non-zeros in row i
        for (int k = A.row_ptr[i]; k < A.row_ptr[i+1]; k++) {
            int col = A.col_indices[k];
            res[i] += A.values[k] * v[col];
        }
    }
    return res;
}


Vector subtract(const Vector& a, const Vector& b) {
    Vector res(a.size());
    for (size_t i = 0; i < a.size(); i++) {
        res[i] = a[i] - b[i];
    }
    return res;
}

Vector adittion(const Vector& a, const Vector& b) {
    Vector res(a.size());
    for (size_t i = 0; i < a.size(); i++) {
        res[i] = a[i] + b[i];
    }
    return res;
}

// Helper: Dot Product (u . v)
double dot(const Vector& u, const Vector& v) {
    return std::inner_product(u.begin(), u.end(), v.begin(), 0.0);
}

// Warm-start overload: delegates to main solve with x0 as initial guess
Vector CG::solve(const CSRMatrix& A, const Vector& b, const Vector& x0) const {
    int n = b.size();
    Vector x = x0;
    Vector r(n), p(n), Ap(n);   // pre-allocated — no heap alloc inside the loop

    // r = b - A*x0
    for (int i = 0; i < n; ++i) {
        double Ax_i = 0.0;
        for (int k = A.row_ptr[i]; k < A.row_ptr[i+1]; ++k)
            Ax_i += A.values[k] * x[A.col_indices[k]];
        r[i] = b[i] - Ax_i;
        p[i] = r[i];
    }

    double r_dot = dot(r, r);
    const double b_norm = sqrt(dot(b, b));
    const double tol_abs = (b_norm > 0.0) ? 1e-6 * b_norm : 1e-10;
    // CG converges in O(sqrt(kappa)) ~ O(N) iters for Poisson; cap at 5*sqrt(n)
    const int max_iter = std::max(200, (int)(5.0 * sqrt((double)n)));

    if (sqrt(r_dot) < tol_abs) return x;

    for (int iter = 1; iter <= max_iter; ++iter) {
        // Ap = A * p  (in-place, no allocation)
        for (int i = 0; i < n; ++i) {
            double s = 0.0;
            for (int k = A.row_ptr[i]; k < A.row_ptr[i+1]; ++k)
                s += A.values[k] * p[A.col_indices[k]];
            Ap[i] = s;
        }
        double pAp = dot(p, Ap);
        if (abs(pAp) < 1e-300) break;
        double alpha = r_dot / pAp;
        for (int i = 0; i < n; ++i) { x[i] += alpha * p[i]; r[i] -= alpha * Ap[i]; }
        double r_dot_new = dot(r, r);
        if (sqrt(r_dot_new) < tol_abs) return x;
        double beta = r_dot_new / r_dot;
        for (int i = 0; i < n; ++i) p[i] = r[i] + beta * p[i];
        r_dot = r_dot_new;
    }
    cout << "CG WARNING: did not converge in " << max_iter
         << " iters. ||r||/||b|| = " << sqrt(r_dot)/b_norm << "\n";
    return x;
}

Vector CG::solve(const CSRMatrix& A, const Vector& b) const {
    int n = b.size();
    Vector x(n, 0.0);
    Vector r = subtract(b, multiply(A, x));  // r = b - A*x (x=0, so r = b)
    Vector p = r;                             // search direction
    double r_dot = dot(r, r);                 // r^T r
    const double b_norm = sqrt(dot(b, b));
    const double tol_abs = (b_norm > 0.0) ? 1e-6 * b_norm : 1e-10;
    int max_iter = 10 * n;                    // O(N) is typical for Poisson

    if (sqrt(r_dot) < tol_abs) {
        cout << "CG: initial residual already below tolerance." << endl;
        return x;
    }

    for (int iter = 1; iter <= max_iter; ++iter) {
        Vector Ap = multiply(A, p);
        double pAp = dot(p, Ap);

        if (abs(pAp) < 1e-300) {
            cout << "CG: breakdown (p^T A p ~ 0) at iter " << iter << endl;
            break;
        }

        double alpha = r_dot / pAp;

        for (int i = 0; i < n; i++) {
            x[i] += alpha * p[i];
            r[i] -= alpha * Ap[i];
        }

        double r_dot_new = dot(r, r);

        if (sqrt(r_dot_new) < tol_abs) {
            return x;
        }

        double beta = r_dot_new / r_dot;
        for (int i = 0; i < n; i++)
            p[i] = r[i] + beta * p[i];

        r_dot = r_dot_new;
    }

    cout << "CG WARNING: did not converge in " << max_iter
         << " iters. ||r||/||b|| = " << sqrt(r_dot)/b_norm << "\n";
    return x;
}

// Warm-start PCG with Jacobi (diagonal) preconditioning.
// Standard PCG iteration: M = diag(A), z = M^{-1} r.
// For non-uniform meshes the diagonal varies → Jacobi lowers condition number.
// All working vectors pre-allocated (no heap inside the loop).
Vector PCG::solve(const CSRMatrix& A, const Vector& b, const Vector& x0) const {
    int n = A.n;
    Vector x = x0;
    Vector r(n), z(n), p(n), Ap(n);

    // r = b - A*x0;  z = M^{-1} r;  p = z
    for (int i = 0; i < n; ++i) {
        double Ax_i = 0.0;
        for (int k = A.row_ptr[i]; k < A.row_ptr[i+1]; ++k)
            Ax_i += A.values[k] * x[A.col_indices[k]];
        r[i] = b[i] - Ax_i;
        z[i] = r[i] / A.diagonal[i];
        p[i] = z[i];
    }

    double rz = dot(r, z);
    const double b_norm = sqrt(dot(b, b));
    const double tol_abs = (b_norm > 0.0) ? 1e-6 * b_norm : 1e-10;
    const int max_iter = std::max(200, (int)(5.0 * sqrt((double)n)));

    if (sqrt(dot(r, r)) < tol_abs) return x;

    for (int iter = 1; iter <= max_iter; ++iter) {
        // Ap = A * p
        for (int i = 0; i < n; ++i) {
            double s = 0.0;
            for (int k = A.row_ptr[i]; k < A.row_ptr[i+1]; ++k)
                s += A.values[k] * p[A.col_indices[k]];
            Ap[i] = s;
        }
        double pAp = dot(p, Ap);
        if (abs(pAp) < 1e-300) break;
        double alpha = rz / pAp;
        for (int i = 0; i < n; ++i) { x[i] += alpha * p[i]; r[i] -= alpha * Ap[i]; }
        if (sqrt(dot(r, r)) < tol_abs) return x;
        for (int i = 0; i < n; ++i) z[i] = r[i] / A.diagonal[i];
        double rz_new = dot(r, z);
        double beta = rz_new / rz;
        for (int i = 0; i < n; ++i) p[i] = z[i] + beta * p[i];
        rz = rz_new;
    }
    cout << "PCG WARNING: did not converge in " << max_iter
         << " iters. ||r||/||b|| = " << sqrt(dot(r,r))/b_norm << "\n";
    return x;
}

Vector PCG::solve(const CSRMatrix& A, const Vector& b) const {
    int n = A.n;
    Vector x(n, 0.0);
    Vector D_inv_sqrt(n);
    for (int i = 0; i < n; i++) {
        D_inv_sqrt[i] = 1.0 / sqrt(A.diagonal[i]); //
    }

    // Computation of preconditioned matrix

    Vector r = subtract(b,multiply(A,x));
    Vector w = r;
    Vector v = w;
    for (int i = 0; i < n; i++) {
    w[i] = D_inv_sqrt[i] * r[i];
    v[i] = D_inv_sqrt[i] * w[i];
    }
    Vector v_new =w;
    double alpha=dot(w,w);
    double beta = 0.0;
    int max_iter = 100;
    double tol = 1e-6;


    for (int iter = 1; iter <= max_iter; ++iter) {
        if(sqrt(dot(v,v))<tol){cout<<"Procedure was successful"<<endl;return x;};
        
        Vector u = multiply(A, v);
        double t = alpha / dot(v, u);
        

        for (int i = 0; i < n; i++) {
            x[i] = x[i] + t * v[i];
            r[i] = r[i] - t * u[i];
            w[i] = D_inv_sqrt[i]*r[i];
            
        }
        beta=dot(w,w);

        if (abs(beta) < tol) {
            cout << "CG converged in " << iter << " iterations." << endl;
            return x;
        }

        double s = beta / alpha;
        for (int i = 0; i < n; i++) {
            v_new[i] = D_inv_sqrt[i]*w[i] + s * v[i];
        }
        alpha=beta;
        v = v_new;
    }
    
    cout << "CG did not converge within max iterations." << endl;
    return x;
}
