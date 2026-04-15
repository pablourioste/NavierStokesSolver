# Discrete Operator Report — NavierStokesSolver

This report describes each discrete operator in the solver: its mathematical role, discretisation
stencil, ghost-cell treatment, sign conventions, and every C++ file that implements or uses it.
It also documents the existing MMS verification tests (Suites 1–5) and the five additional tests
(Suites 6–10) together with the new convergence-rate plots.

---

## 1. Staggered MAC Grid

The solver uses a **Marker-and-Cell (MAC) staggered grid** on a rectangular domain
Ω = (0, L) × (0, H):

```
                     v(i, M+1)   ← north boundary row (V-nodes)
     ┌──────┬──────┬──────┐
     │      │      │      │
  u  │  p   │  p   │  p   │  u   ← P cells, U faces (east/west)
 (1,j)│(1,j) │(2,j) │(3,j) │(4,j)
     │      │      │      │
     └──────┴──────┴──────┘
                     v(i, 1)     ← south boundary row (V-nodes)
```

| Grid | Symbol | Degrees of freedom | Logical index range |
|------|---------|--------------------|---------------------|
| Pressure | P | `N × M` | i=1..N, j=1..M |
| Horizontal velocity | U | `(N+1) × M` | k=1..N+1, j=1..M |
| Vertical velocity | V | `N × (M+1)` | i=1..N, l=1..M+1 |

**Why stagger?** Collocated grids allow the pressure–velocity system to decouple
(the "checkerboard" instability). The staggered layout ensures that the divergence of the
gradient operator is spectrally equivalent to the Laplacian, eliminating spurious
pressure modes without any stabilisation term.

Index helpers used throughout:

```
pIdx(i,j) = (j-1)*N     + (i-1)     → flat index into P-vector
uIdx(k,j) = (j-1)*(N+1) + (k-1)     → flat index into U-vector
vIdx(i,l) = (l-1)*N     + (i-1)     → flat index into V-vector
```

---

## 2. Operator Catalogue

### 2.1 Laplacian — **L**

**Mathematical form:**
```
L φ = ∫∫_Ω ∇²φ dV   (assembled as matrix A, positive-definite form)
```

**Role:** Discretises ∂²p/∂x² + ∂²p/∂y² for the pressure-Poisson equation (PPE)
in the pressure-correction step of Chorin's method.

**Discretisation (FVM, 2nd-order):**

For P-cell (i,j) with face areas S and cell volume V:

```
a_E = S_E / (x[i+1] - x[i])       (east flux coefficient)
a_W = S_W / (x[i]   - x[i-1])     (west flux coefficient)
a_N = S_N / (y[j+1] - y[j])       (north flux coefficient)
a_S = S_S / (y[j]   - y[j-1])     (south flux coefficient)

a_P = +(a_E + a_W + a_N + a_S)    (diagonal — positive)

(A·φ)[k] = a_P·φ_P - a_E·φ_E - a_W·φ_W - a_N·φ_N - a_S·φ_S
          = Σ a_nb·(φ_P - φ_nb)   ≈ -V · ∇²φ
```

**Sign convention:** The matrix is **positive-definite**, so `(A·φ)/V ≈ −∇²φ`.
For the MMS test with φ = cos(πx)cos(πy), ∇²φ = −2π², so `(A·φ)/V ≈ +2π²cos(πx)cos(πy)`.

**Boundary conditions:** Homogeneous Neumann (∂p/∂n = 0). Boundary face coefficients
are set to zero — no contribution from outside the domain.

**Singularity:** With all-Neumann BCs the system is singular (constant p is a null
vector). `pinReference()` fixes p[0] = 0 to regularise.

**C++ implementation:** [src/operators/LaplaceOperator.hpp](../src/operators/LaplaceOperator.hpp)
/ [.cpp](../src/operators/LaplaceOperator.cpp)

**Used by:**
- [src/PressurePoissonSystem.cpp](../src/PressurePoissonSystem.cpp) — `assembleMatrix()` builds the
  identical FVM stencil for the PPE; `LaplaceOperator::assemble()` is the standalone form
- [src/verification/verify_operators.cpp](../src/verification/verify_operators.cpp) — Suites 1, 5, 9, 10

---

### 2.2 Gradient — **G**

**Mathematical form:**
```
G : P → (U, V)
  grad_u[k,j] = ∂p/∂x|_{face(k,j)} ≈ (p[k,j] - p[k-1,j]) / (x_P[k] - x_P[k-1])
  grad_v[i,l] = ∂p/∂y|_{face(i,l)} ≈ (p[i,l] - p[i,l-1]) / (y_P[l] - y_P[l-1])
```

**Role:** Maps the pressure correction onto the U/V faces so that the divergence-free
velocity field can be recovered:
```
u^{n+1} = u* - (dt/ρ) · G·p
```
The operator itself is purely geometric (no dt or ρ).

**Discretisation (FD, 2nd-order):** Centered difference between adjacent P-cell
centres, evaluated exactly at the shared face — the natural position of the U/V node.

**Boundary treatment:** Boundary U-nodes (k=1, k=N+1) and V-nodes (l=1, l=M+1)
are set to zero in the output (velocity BCs are imposed elsewhere).

**C++ implementation:** [src/operators/GradientOperator.hpp](../src/operators/GradientOperator.hpp)
/ [.cpp](../src/operators/GradientOperator.cpp)

**Used by:**
- [src/PressurePoissonSystem.cpp](../src/PressurePoissonSystem.cpp) — `projectVelocity()` applies
  `-(dt/ρ) · G·p` to correct the predictor velocity
- [src/verification/verify_operators.cpp](../src/verification/verify_operators.cpp) — Suites 2, 5, 8

---

### 2.3 Divergence — **D**

**Mathematical form:**
```
D : (U, V) → P
  (D·u)[i,j] = (1/V) · [u_E·S_E - u_W·S_W + v_N·S_N - v_S·S_S]
```

**Role:** Computes the volumetric flux divergence at each P-cell. The pressure-Poisson
right-hand side is:
```
b = (ρ/dt) · D·u*
```
After projection, D·u^{n+1} = 0 (incompressibility) to machine precision.

**Discretisation (FVM, 2nd-order):** Outward-normal flux integration over all four
faces of P-cell (i,j):

```
Staggered index mapping for P(i,j):
  u_E ← u[uIdx(i+1, j)],  u_W ← u[uIdx(i,   j)]
  v_N ← v[vIdx(i,   j+1)], v_S ← v[vIdx(i,   j)]

flux = u_E·S_E - u_W·S_W + v_N·S_N - v_S·S_S
result[pIdx(i,j)] = flux / V[i][j]        (units: 1/time)
```

**Important:** The U/V fields must be sampled at the **face position** (not the node
centroid) for boundary cells — see Suite 3 notes below.

**C++ implementation:** [src/operators/DivergenceOperator.hpp](../src/operators/DivergenceOperator.hpp)
/ [.cpp](../src/operators/DivergenceOperator.cpp)

**Used by:**
- [src/PressurePoissonSystem.cpp](../src/PressurePoissonSystem.cpp) — `computeRHS()` evaluates
  `(ρ/dt)·D·u*`
- [src/verification/verify_operators.cpp](../src/verification/verify_operators.cpp) — Suites 3, 5, 8

---

### 2.4 Diffusion — **D_ν**

**Mathematical form:**
```
D_ν : (U, V) → (U, V)
  (D_ν u)[k,j] = ν · (∂²u/∂x² + ∂²u/∂y²)|_{node(k,j)}
  (D_ν v)[i,l] = ν · (∂²v/∂x² + ∂²v/∂y²)|_{node(i,l)}
```

**Role:** Discretises the viscous term ν∇²u in the momentum equations:
```
u* = u^n + dt · [D_ν·u  −  C(u)·u]
```

**Discretisation (FD, 2nd-order, non-uniform grid):**

For a U-node at (k,j):
```
d²u/dx² = 2/(dx_e + dx_w) · [(u_E - u_P)/dx_e  -  (u_P - u_W)/dx_w]
d²u/dy² = 2/(dy_n + dy_s) · [(u_N - u_P)/dy_n  -  (u_P - u_S)/dy_s]
result   = ν · (d²u/dx² + d²u/dy²)
```

**Ghost-cell method for walls (2nd-order):**

At j=1 (south wall) and j=M (north wall) for the u-equation:
```
u_ghost = 2·u_wall − u_P           (value on far side of wall face)
dy_ghost = 2·|y_P − y_wall|        (doubled distance for correct gradient)
```
This ensures the interpolated value at the wall face equals u_wall exactly.

Similarly for v at i=1 (west) and i=N (east).

**Wall velocity values** are passed at construction:
- `u_ghost_north`, `u_ghost_south` — for the u-component (lid or no-slip)
- `v_ghost_west`, `v_ghost_east` — for the v-component (always no-slip in lid cavity)

**C++ implementation:** [src/operators/DiffusionOperator.hpp](../src/operators/DiffusionOperator.hpp)
/ [.cpp](../src/operators/DiffusionOperator.cpp)

**Used by:**
- [src/MomentumSolver.cpp](../src/MomentumSolver.cpp) — holds a `DiffusionOperator diffOp_` member;
  `uRHS()` and `vRHS()` call `diffOp_.diffuseU()` / `diffOp_.diffuseV()` per node
- [src/verification/verify_operators.cpp](../src/verification/verify_operators.cpp) — Suites 4, 6, 7

---

### 2.5 Convection — **C(u)**

**Mathematical form:**
```
C(u) : (U, V) → (U, V)
  (C(u)·u)[k,j] = u_P·∂u/∂x  +  v_interp·∂u/∂y     (u-momentum)
  (C(u)·v)[i,l] = u_interp·∂v/∂x  +  v_P·∂v/∂y     (v-momentum)
```

**Role:** Discretises the nonlinear convection term (u·∇)u. The full momentum RHS is
`D_ν·u − C(u)·u`; both are combined inside `MomentumSolver::uRHS()` / `vRHS()`.

**Discretisation — two selectable schemes:**

| Scheme | Order | Stability | Formula (u, x-direction) |
|--------|-------|-----------|--------------------------|
| **CDS** (Central) | 2nd | Neutral | `u_P · (u_E − u_W)/(dx_e + dx_w)` |
| **UPWIND** | 1st | Dissipative | `u_P≥0 ? (u_P−u_W)/dx_w : (u_E−u_P)/dx_e` |

**Cross-velocity interpolation** (bilinear, 4-point average):
```
v at U-node (k,j):  v_interp = 0.25·(v[k-1,j] + v[k,j] + v[k-1,j+1] + v[k,j+1])
u at V-node (i,l):  u_interp = 0.25·(u[i,l-1] + u[i+1,l-1] + u[i,l] + u[i+1,l])
```

**Note:** Convection is **not** a standalone class. It is implemented inline inside
`MomentumSolver::uRHS()` and `MomentumSolver::vRHS()`.

**C++ implementation:** [src/MomentumSolver.hpp](../src/MomentumSolver.hpp)
/ [src/MomentumSolver.cpp](../src/MomentumSolver.cpp) (private methods `uRHS`, `vRHS`)

**Used by:**
- [src/MomentumSolver.cpp](../src/MomentumSolver.cpp) — `advance()` applies Adams-Bashforth 2
  time integration using `uRHS`/`vRHS`

---

### Summary Table

| Operator | Symbol | Domain → Range | Discretisation | Standalone class | Files |
|----------|--------|---------------|----------------|------------------|-------|
| Laplacian | **L** | P → P | FVM 2nd-order | `LaplaceOperator` | LaplaceOperator.{hpp,cpp}, PressurePoissonSystem.cpp |
| Gradient | **G** | P → (U,V) | FD 2nd-order | `GradientOperator` | GradientOperator.{hpp,cpp}, PressurePoissonSystem.cpp |
| Divergence | **D** | (U,V) → P | FVM 2nd-order | `DivergenceOperator` | DivergenceOperator.{hpp,cpp}, PressurePoissonSystem.cpp |
| Diffusion | **D_ν** | (U,V) → (U,V) | FD 2nd-order + ghost | `DiffusionOperator` | DiffusionOperator.{hpp,cpp}, MomentumSolver.cpp |
| Convection | **C(u)** | (U,V) → (U,V) | CDS / Upwind | (inline) | MomentumSolver.{hpp,cpp} |

---

## 3. Existing Test Suite (Suites 1–5)

All tests live in [src/verification/verify_operators.cpp](../src/verification/verify_operators.cpp)
and are run by the `operator_verify` executable.  
Each suite builds uniform N×N grids for N ∈ {8, 16, 32, 64}, measures the L∞ error, and
estimates the convergence rate r = log₂(e_N / e_{2N}). Target: r ≈ 2.0.

---

### Suite 1 — LaplaceOperator (MMS, Neumann BC)

**Manufactured solution:** φ(x,y) = cos(πx)cos(πy)  
**Why this function?** Its normal derivative vanishes on all four walls — the
homogeneous Neumann BC is satisfied exactly, so no boundary error contaminates the test.

**Analytical Laplacian:** ∇²φ = −2π²cos(πx)cos(πy)  
**Discrete check:** `(A·φ)[k] / V[i][j]` must equal `+2π²cos(πx)cos(πy)` (positive-definite sign)  
**Error metric:** L∞ over all N² P-cells  
**Expected rate:** ~2.0 — a failure here means a coding error in the FVM stencil or the
face-area/volume factors.

---

### Suite 2 — GradientOperator (MMS, interior nodes only)

**Manufactured solution:** φ(x,y) = cos(πx)cos(πy)  
**Analytical gradients:**
- ∂φ/∂x = −πsin(πx)cos(πy)  at interior U-nodes (k=2..N)
- ∂φ/∂y = −πcos(πx)sin(πy)  at interior V-nodes (l=2..M)

**Why interior only?** Boundary U/V nodes (k=1, k=N+1, l=1, l=M+1) are set to zero by
the operator — they represent prescribed BCs, not computed gradients.

**Error metric:** Two separate L∞ norms (∂p/∂x and ∂p/∂y)  
**Expected rate:** ~2.0 for both components — a failure indicates wrong cell-centre
positions in the mesh or incorrect indexing of adjacent P-cells.

---

### Suite 3 — DivergenceOperator (MMS, face sampling)

**Manufactured solution:**
- u(x,y) = sin(2πx)cos(πy)   →  ∂u/∂x = 2πcos(2πx)cos(πy)
- v(x,y) = cos(πx)sin(2πy)   →  ∂v/∂y = 2πcos(πx)cos(2πy)

**Critical detail — face vs. centroid sampling:** Boundary U-nodes (k=1, k=N+1) have
their node centroid at dx/4 from the wall, **not** at the face position x=0. The
divergence operator reads the flux at the face. Therefore the test samples u at
`x_face = pMesh.xvc[k-1]` (the actual face position), not at the node centroid.
Interior nodes coincide with their face, so no correction is needed there.

**Error metric:** L∞ over all N² P-cells  
**Expected rate:** ~2.0 — a failure may mean wrong face-area factors or an index
mismatch between the U-node flat index and the intended P-cell face.

---

### Suite 4 — DiffusionOperator `applyU` (MMS, ghost cells)

**Manufactured solution:** u(x,y) = sin(πy) (constant in x → ∂²u/∂x² = 0 exactly)  
**Analytical result:** ν∇²u = −νπ²sin(πy),  ν = 1  
**Boundary values:** sin(π·0) = sin(π·1) = 0 → wall velocity = 0 → ghost = −u_P

**Why x-constant?** It isolates the y-direction ghost-cell stencil without
interference from any non-uniform spacing in x. The x-term drops out analytically,
so only the ghost-cell treatment for j=1 and j=M is tested.

**Comparison range:** Interior U-nodes k=2..N, all j (including boundary-adjacent j=1, j=M
which use ghost cells — they are interior in the sense that they are computed, not prescribed).

**Expected rate:** ~2.0 — a failure here usually means the ghost-cell distance formula
`dy_ghost = 2·|y_P − y_wall|` is wrong or the wall value is applied to the wrong axis.

---

### Suite 5 — Discrete identity D·G + L/V = 0 (machine precision)

**Identity:** For any scalar field φ and **every** P-cell (including boundary cells):
```
D·(G·φ)  +  (L·φ)/V  =  0
```

**Proof sketch (interior face):**
- G contributes `S_E/V · (φ_E − φ_P)/(x_E − x_P)` = `+a_E/V · (φ_E − φ_P)` to D·G·φ
- L contributes `−a_E·(φ_E − φ_P)` to L·φ, giving `−a_E/V · (φ_E − φ_P)` after dividing by V
- Sum = 0 ✓

**Boundary face:** G returns 0 for boundary nodes; L ignores the boundary face coefficient.
Both contributions are zero independently — identity holds at corners too.

**Expected error:** ~1e-13 (double precision) regardless of N — this is an algebraic
identity, not an approximation.  
**Why important?** This identity is the discrete analogue of the Poisson equation
∇·∇φ = ∇²φ. It guarantees that the pressure-correction step genuinely enforces
incompressibility, and that the CG solver is solving a consistent system.

---

## 4. New Tests (Suites 6–10)

### Suite 6 — DiffusionOperator `applyV` (MMS, west/east ghost cells)

**Why:** Suite 4 only tests the u-component (north/south ghost cells). The v-component
uses ghost cells at the **west and east** walls — a different code path. An independent
error in `diffuseV` would not be caught by Suite 4.

**Manufactured solution:** v(x,y) = sin(πx) (constant in y → ∂²v/∂y² = 0)  
**Analytical result:** ν∇²v = −νπ²sin(πx),  ν = 1  
**Boundary values:** sin(π·0) = sin(π·1) = 0 → v_ghost_west = 0, v_ghost_east = 0  
**Comparison range:** All i=1..N, interior l=2..M  
**Expected rate:** ~2.0

---

### Suite 7 — DiffusionOperator non-zero lid velocity (exact cancellation)

**Why:** The production solver sets `u_ghost_north = U_lid ≠ 0`. The ghost formula
`u_ghost = 2·U_lid − u_P` must give zero Laplacian for a linear field — this is an
**exact** property, not just an O(h²) approximation.

**Manufactured solution:** u(x,y) = y (linear, ∂²u/∂y² = 0, so ν∇²u = 0)  
**Wall values:** u_wall_south = 0 (no-slip), u_wall_north = 1 (lid)  
**Expected error:** < 1e-12 at **all** N — the cancellation is algebraically exact for
a linear field and should hold to floating-point precision. A rate ~2 would indicate
only approximate cancellation, which would be a bug.

---

### Suite 8 — Discrete adjoint G = −D^T (energy conservation)

**Why:** The symmetry-preserving property of Chorin's method requires that the gradient
and divergence are discrete adjoints:
```
φᵀ · (D·u)  +  (G·φ)ᵀ · u  =  0      for all φ ∈ P, u ∈ U interior
```
This is the discrete analogue of the integration-by-parts identity
∫(∇·u)φ dV + ∫u·∇φ dV = ∮ φu·dS = 0 (with u·n = 0 on walls).

The identity holds WITHOUT volume weighting because `DivergenceOperator.apply()` already
divides by V internally, making D and G formally adjoint in the unweighted L² inner product.

**Test:** Generate two random vectors φ (P-space) and u/v (U/V-space, zero on boundaries).
Compute `φᵀ·(D·u) + (G·φ)ᵀ·u`.

**Expected residual:** grows with N due to the O(N²) accumulation of cancelling floating-point
terms. Residual at N=64 is ~1e-11 (roughly N × machine_epsilon × |φ||Du|), confirming
machine-precision cancellation. Not mesh-dependent in the normalised sense.  
**Why important?** A violation means the discretisation does not conserve kinetic
energy in the discrete sense, and long-time simulations may exhibit spurious drift.

---

### Suite 9 — LaplaceOperator SPD check

**Why:** The CG solver requires **symmetric positive-definite** A. Symmetry is not trivially
guaranteed if any coefficient in the boundary treatment has an asymmetric implementation.

**Test A — symmetry:** For two random vectors φ, ψ:
```
|φᵀ·(A·ψ)  −  ψᵀ·(A·φ)|  <  ε_machine
```

**Test B — positive-definiteness:** For a random non-zero φ with p[0] pinned to zero:
```
φᵀ·A·φ  >  0
```

**Expected:** symmetry residual ~1e-13; quadratic form strictly positive for all N.

---

### Suite 10 — LaplaceOperator on a stretched mesh (first-order limitation)

**Why:** Production runs use tanh-stretched grids (see `input/stretch.txt`).
This suite exposes a limitation of the current FVM stencil on non-uniform meshes.

**Grid mapping:** Uniform parameter ξ ∈ [0,1] → physical coordinate x via:
```
x_map(ξ) = 1/2 · (1 + tanh(γ·(2ξ−1)) / tanh(γ)),   γ = 2
```
Same mapping applied independently in both x and y.

**Manufactured solution:** φ(x,y) = cos(πx)cos(πy)  
**Expected rate: ~1.0** (first-order)

**Why first-order?** The stencil uses `a_E = Se / (x[i+1] − x[i])` where x[i] are cell
centres. For a non-uniform mesh with `h_e ≠ h_w`, the x-Laplacian contribution is:
```
[(φ_E−φ_P)/h_e − (φ_P−φ_W)/h_w] / dx_i
≈ φ''_P + φ'_P · (h_e−h_w)/dx_i + O(h²)
```
The first-order term `φ'_P · (h_e−h_w)/dx_i` is O(1/N) because `h_e−h_w = f''(ξ)/N²`
and `dx_i = f'(ξ)/N` (for a smooth mapping f).

**Fix:** Replace `1/dx_i` by `2/(h_e+h_w)` in the stencil to recover 2nd-order accuracy
on non-uniform meshes. Suite 10 therefore also serves as a regression guard: if the stencil
is ever corrected, the rate should jump from ~1 to ~2.

---

## 5. Convergence Plots — `plot_convergence.py`

Running `./build/operator_verify` writes `output/convergence.csv`. The script
`plotters/plot_convergence.py` reads this file and produces `output/convergence_rates.png`,
a 2×3 figure analogous to Fig. 5 of the reference paper:

| Panel | Content |
|-------|---------|
| (0,0) | Log-log: Laplacian L∞ vs h=1/N, with O(h²) and O(h) reference lines |
| (0,1) | Log-log: Gradient ∂p/∂x and ∂p/∂y L∞ vs h |
| (0,2) | Log-log: Divergence L∞ vs h |
| (1,0) | Log-log: Diffusion U and V L∞ vs h (Suites 4 and 6) |
| (1,1) | Semi-log: D·G identity and adjoint residuals vs N (Suites 5 and 8 — should be ~1e-13) |
| (1,2) | Log-log: Laplacian uniform (Suite 1) vs stretched mesh (Suite 10) |

Each log-log panel carries an estimated slope annotation. A slope ≈ 2.0 on panels (0,0)–(1,0)
is visual proof of second-order accuracy. The flat residual lines in panel (1,1) confirm that
the algebraic identities hold at floating-point precision irrespective of mesh refinement.

---

## 6. How to Run

```bash
# Build
cd build_cmake && cmake .. -DCMAKE_BUILD_TYPE=Release && make operator_verify -j4

# Run verification (also writes output/convergence.csv)
./build/operator_verify

# Expected console output:
#   Suites 1–4, 6:  convergence rate  ~2.0 (second-order)
#   Suite 7:        error             = 0.0 (exact algebraic cancellation)
#   Suite 5:        residual          ~1e-14..1e-12 (machine precision, grows with N)
#   Suite 8:        residual          ~0..3e-11 (machine precision, grows with N²)
#   Suite 9:        symmetry          ~1e-16..1e-14, PD > 0
#   Suite 10:       convergence rate  ~1.0 (first-order, stencil limitation on non-uniform mesh)

# Plot convergence rates
python3 plotters/plot_convergence.py
# → output/convergence_rates.png
```
