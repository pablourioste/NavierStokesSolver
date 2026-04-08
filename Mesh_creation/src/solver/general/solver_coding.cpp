#include "solver_class.hpp"
#include <iostream>
#include <cmath>
#include <stdexcept>
#include <vector>
#include <algorithm>
#include <numeric>
using namespace std;

// Jacobi Iteration Solver
Vector Jacobi::solve(const Matrix& A, const Vector& b) const {
    int n = b.size();
    Vector x(n, 0.0);
    Vector x_prev(n,0.0);
    int max_iter =100;
    double tol=1e-6;

   for (int iter = 1; iter <= max_iter; ++iter) {
        x_prev=x;
        for (int i = 0; i < n; i++) {
            double term_A=0.0;

            if (abs(A[i][i]) < 1e-15) {
                cerr << "Error: Zero diagonal at row " << i << endl;
                return x; 
            }

            for (int j=0; j<n; j++){
                if (j!=i){ term_A+= A[i][j]*x_prev[j];}
            }

            x[i]=1/(A[i][i])*(-term_A+b[i]);
        }
            double max_error =0.0;
        for (int i=0; i<n; i++){
            max_error=max(max_error,abs(x[i]-x_prev[i]));
        }
        if (max_error < tol) {
        cout << "We break the loop because convergence has been achieved" << endl;
        cout << "Number of iteration needed: "<<iter<<endl;
        break; 
        }

        if (iter>max_iter){cout<<"The procedure was not succesful"<<endl;}
        
    }  

  
    return x;
}


// Gauss Seidel (GS) iterative solver.
Vector GS::solve(const Matrix& A, const Vector& b) const {
    int n = b.size();
    Vector x(n, 0.0);
    Vector x_prev(n,0.0);
    int max_iter =100;
    double tol=1e-6;

   for (int iter = 1; iter <= max_iter; ++iter) {
        x_prev=x;
        for (int i = 0; i < n; i++) {
            double term_A=0.0;
            double term_B=0.0;

            if (abs(A[i][i]) < 1e-15) {
                cerr << "Error: Zero diagonal at row " << i << endl;
                return x; 
            }

            for (int j=0; j<i; j++){// Unknowns already updated
            term_A+= A[i][j]*x[j];
            }

            for (int j=i+1; j<n;j++){ // Unknown still not updated from previous iteration
            term_B+=A[i][j]*x_prev[j];
            }
           
            x[i]=1/(A[i][i])*(-term_A-term_B+b[i]);
        }
    

        double max_error =0.0;
        for (int i=0; i<n; i++){
            max_error=max(max_error,abs(x[i]-x_prev[i]));
        }
        if (max_error < tol) {
        cout << "We break the loop because convergence has been achieved" << endl;
        cout << "Number of iteration needed: "<<iter<<endl;
        break; 
        }

        if (iter>max_iter){cout<<"The procedure was not succesful"<<endl;}
        
    }  

  
    return x;
}

Vector SOR::solve(const Matrix& A, const Vector& b) const {
    int n = b.size();
    Vector x(n, 0.0);
    Vector x_prev(n,0.0);
    int max_iter =100;
    double tol=1e-6;

    //double omega=getOptimalOmega(A);
    double omega=1.25;
    cout<<"Optimal relaxation factor: "<<omega<<endl;
    for (int iter = 1; iter <= max_iter; ++iter) {
        x_prev=x;
        for (int i = 0; i < n; i++) {
            double term_A=0.0;
            double term_B=0.0;

            if (abs(A[i][i]) < 1e-15) {
                cerr << "Error: Zero diagonal at row " << i << endl;
                return x; 
            }

            for (int j=0; j<i; j++){// Unknowns already updated
            term_A+= A[i][j]*x[j];
            }

            for (int j=i+1; j<n;j++){ // Unknown still not updated from previous iteration
            term_B+=A[i][j]*x_prev[j];
            }
           
            x[i]=(1-omega)*x_prev[i]+omega/(A[i][i])*(-term_A-term_B+b[i]);
        }
    

        double max_error =0.0;
        for (int i=0; i<n; i++){
            max_error=max(max_error,abs(x[i]-x_prev[i]));
        }
        if (max_error < tol) {
        cout << "We break the loop because convergence has been achieved" << endl;
        cout << "Number of iteration needed: "<<iter<<endl;
        break; 
        }

        if (iter>max_iter){cout<<"The procedure was not succesful"<<endl;}
        
    }  
    return x;
}



double SOR::getOptimalOmega(const Matrix& A) const{
    int n = A.size();
    std::vector<double> v(n, 1.0); // Start with a guess vector
    std::vector<double> next_v(n, 0.0);
    double rho_prev = 0.0;
    double rho = 0.0;
    
    // Power Method to find the spectral radius of the Jacobi matrix
    for (int iter = 0; iter < 100; ++iter) {
        // Compute next_v = T_j * v
        for (int i = 0; i < n; ++i) {
            double sum = 0.0;
            for (int j = 0; j < n; ++j) {
                if (i != j) {
                    sum += (-A[i][j] / A[i][i]) * v[j]; // Definition of T_j
                }
            }
            next_v[i] = sum;
        }

        // Calculate the norm for the eigenvalue estimate
        double norm = 0.0;
        for (double val : next_v) norm += val * val;
        norm = sqrt(norm);

        // Normalize and estimate rho
        rho_prev = rho;
        rho = norm; // Simplified estimate for the largest eigenvalue
        
        for (int i = 0; i < n; ++i) v[i] = next_v[i] / norm;

        if (std::abs(rho - rho_prev) < 1e-6) break;
    }

    // Step 2 from your image: Apply the optimal omega formula
    double rho_sq = pow(rho, 2);
    if (rho_sq >= 1.0) return 1.0; // SOR only converges for rho < 1
    
    return 2.0 / (1.0 + sqrt(1.0 - rho_sq));
}


Vector multiply(const Matrix& A, const Vector& v) {
    int n = A.size();
    Vector res(n, 0.0);
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            res[i] += A[i][j] * v[j];
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

Vector CG::solve(const Matrix& A, const Vector& b) const {
    int n = b.size();
    Vector x(n, 0.0);
    Vector r = subtract(b,multiply(A,x));
    Vector w = r;
    Vector v = w;
    Vector v_new =w;
    double alpha=dot(w,w);
    double beta = 0.0;
    int max_iter = 10000;
    double tol = 1e-6;
    
  
    

    for (int iter = 1; iter <= max_iter; ++iter) {
        if(sqrt(dot(v,v))<tol){cout<<"Procedure was successful"<<endl;return x;};
        
        Vector u = multiply(A, v);
        double t = alpha / dot(v, u);
        

        for (int i = 0; i < n; i++) {
            x[i] = x[i] + t * v[i];
            r[i] = r[i] - t * u[i];
            w[i] = r[i];
            
        }
        beta=dot(w,w);

        if (abs(beta) < tol) {
            cout << "CG converged in " << iter << " iterations." << endl;
            return x;
        }

        double s = beta / alpha;
        for (int i = 0; i < n; i++) {
            v_new[i] = w[i] + s * v[i];
        }
        alpha=beta;
        v = v_new;
    }
    
    cout << "CG did not converge within max iterations." << endl;
    return x;
}

Vector PCG::solve(const Matrix& A, const Vector& b) const {
    int n = b.size();
    Vector x(n, 0.0);
    Vector D_inv_sqrt(n);
    for (int i = 0; i < n; i++) {
        D_inv_sqrt[i] = 1.0 / sqrt(A[i][i]); //
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
    int max_iter = 10000;
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
