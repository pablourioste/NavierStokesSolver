#include <iostream>
#include <vector>
#include <memory>
#include <chrono>
#include <iomanip>
#include <string>
#include "solver_class_sparse.hpp" 

using namespace std;
using namespace std::chrono;

void printVector(const Vector& v) {
    cout << "[ ";
    for (size_t i = 0; i < v.size(); ++i) {
        cout << fixed << setprecision(7) << v[i] << (i == v.size() - 1 ? "" : ", ");
    }
    cout << " ]";
}

int main() {
    int N = 300;
        vector<vector<double>> denseA(N, vector<double>(N, 0.0));
        Vector b(N);

    for (int i = 0; i < N; i++) {
        // Diagonal principal potente para asegurar convergencia (SPD)
        denseA[i][i] = 4.0; 
        
        // Sub-diagonal y Super-diagonal (Simetría)
        if (i > 0) denseA[i][i-1] = -1.0;
        if (i < N-1) denseA[i][i+1] = -1.0;
        
        // Acoplamientos lejanos (para que no sea solo una línea)
        if (i > 4) denseA[i][i-5] = -0.5;
        if (i < N-5) denseA[i][i+5] = -0.5;

        // Vector b (por ejemplo, b[i] = i)
        b[i] = (double)i;
    }

    

    // 2. Rearrange into Sparse (CSR) Configuration
    CSRMatrix A = CSRMatrix::fromDense(denseA);

  

    // 3. Setup Solvers
    struct SolverEntry {
        string name;
        unique_ptr<LinearSolver> ptr;
    };

    vector<SolverEntry> solvers;
    solvers.push_back({"Jacobi", make_unique<Jacobi>()});
    solvers.push_back({"Gauss-Seidel", make_unique<GS>()});
    solvers.push_back({"SOR", make_unique<SOR>()});
    solvers.push_back({"CG", make_unique<CG>()}); 
    solvers.push_back({"PCG", make_unique<PCG>()});

    cout << "===============================================================" << endl;
    cout << "           SPARSE (CSR) LINEAR SYSTEM SOLVER COMPARISON        " << endl;
    cout << "===============================================================" << endl;
    
    cout << left << setw(15) << "Solver" 
         << setw(18) << "Time (ms)" 
         << "Solution Vector" << endl;
    cout << "---------------------------------------------------------------" << endl;

    for (auto& entry : solvers) {
        auto start = high_resolution_clock::now();
        
        // Use the new CSR solve method
        Vector x = entry.ptr->solve(A, b);
        
        auto stop = high_resolution_clock::now();
        auto duration = duration_cast<microseconds>(stop - start);

        cout << left << setw(15) << entry.name 
             << setw(18) << (duration.count() / 1000.0);
        //printVector(x);
        cout << endl;
    }

    cout << "===============================================================" << endl;

    return 0;
}