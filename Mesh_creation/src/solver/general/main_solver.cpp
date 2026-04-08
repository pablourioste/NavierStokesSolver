#include <iostream>
#include <vector>
#include <memory>
#include <chrono>
#include <iomanip>
#include <string>
#include "solver_class.hpp" 

using namespace std;
using namespace std::chrono;

// Helper to print the vector in a compact format
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
     Matrix A=denseA;

    //Matrix A = {
    //{ 60.0,  15.5,  20.1,  12.3,   5.1},
    //{ 12.5,  80.0, -10.2,  18.4, -11.0},
    //{  9.1, -12.3,  90.0,  22.5, -15.2},
    //{ 14.2,  11.1,  19.8, 120.0,  40.4},
    //{  5.5, -10.1, -22.2,  34.1, 700.0}
    //};

    //Vector b = {1.0, 2.0, 3.0, 4.0, 5.0};

    // List of solvers to test
    struct SolverEntry {
        string name;
        unique_ptr<LinearSolver> ptr;
    };

    vector<SolverEntry> solvers;
    solvers.push_back({"Jacobi", make_unique<Jacobi>()});
    solvers.push_back({"Gauss-Seidel", make_unique<GS>()});
    solvers.push_back({"SOR", make_unique<SOR>()});
    solvers.push_back({"CG", make_unique<CG>()});   // Uncomment when implemented
    solvers.push_back({"PCG", make_unique<PCG>()}); // Uncomment when implemented

    cout << "===============================================================" << endl;
    cout << "                 LINEAR SYSTEM SOLVER COMPARISON               " << endl;
    cout << "===============================================================" << endl;
    
    // Table Header
    cout << left << setw(15) << "Solver" 
         << setw(18) << "Time (ms)" 
         << "Solution Vector" << endl;
    cout << "---------------------------------------------------------------" << endl;

    for (auto& entry : solvers) {
        // Start timing
        auto start = high_resolution_clock::now();
        
        // Solve (The iteration count will be printed by your internal class methods)
        Vector x = entry.ptr->solve(A, b);
        
        // End timing
        auto stop = high_resolution_clock::now();
        auto duration = duration_cast<microseconds>(stop - start);

        // Print table row
        cout << left << setw(15) << entry.name 
             << setw(18) << (duration.count() / 1000.0);
        //printVector(x);
        cout << endl;
    }

    cout << "===============================================================" << endl;

    return 0;
}