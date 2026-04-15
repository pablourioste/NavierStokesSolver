# NavierStokesSolver — Documentación del Algoritmo

Simulador 2D de flujo incompresible (cavidad con tapa deslizante) mediante el **método de proyección de Chorin** sobre una malla MAC escalonada (staggered grid).

---

## Tabla de Contenidos

1. [Visión general](#1-visión-general)
2. [Estructura de archivos](#2-estructura-de-archivos)
3. [La malla MAC escalonada](#3-la-malla-mac-escalonada)
4. [Estructuras de datos principales](#4-estructuras-de-datos-principales)
5. [Algoritmo global: método de proyección de Chorin](#5-algoritmo-global-método-de-proyección-de-chorin)
6. [Generación de malla — MeshGenerator](#6-generación-de-malla--meshgenerator)
7. [Operadores discretos — src/operators/](#7-operadores-discretos--srcoperators)
8. [Paso predictor — MomentumSolver](#8-paso-predictor--momentumsolver)
9. [Paso corrector — PressurePoissonSystem](#9-paso-corrector--pressurepoissonssystem)
10. [Solucionadores lineales — solver_sparse](#10-solucionadores-lineales--solver_sparse)
11. [Bucle principal — main.cpp](#11-bucle-principal--maincpp)
12. [Concatenación de funciones (grafo de llamadas)](#12-concatenación-de-funciones-grafo-de-llamadas)
13. [Verificación de operadores](#13-verificación-de-operadores)
14. [Postproceso y validación](#14-postproceso-y-validación)
15. [Parámetros numéricos y estabilidad](#15-parámetros-numéricos-y-estabilidad)

---

## 1. Visión general

El solver resuelve las ecuaciones de Navier-Stokes para flujo incompresible en 2D:

```
∂u/∂t + (u·∇)u = -∇p/ρ + ν∇²u
∇·u = 0
```

**Caso de prueba:** Cavidad cuadrada con tapa deslizante (*lid-driven cavity*).
- Dominio: L = H = 1.0 m
- Tapa superior (norte): velocidad prescrita u = U_lid
- Resto de paredes: no-slip (u = v = 0)
- Validación: datos de referencia de Ghia et al. (1982), Re ∈ {100, 400, 1000, …, 10000}

**Método numérico clave:**
- Discretización espacial: diferencias finitas / volúmenes finitos de 2.º orden sobre malla MAC escalonada
- Integración temporal: Adams-Bashforth 2 explícito (AB2), con primer paso en Euler
- Acoplamiento presión-velocidad: proyección de Chorin (fraccionamiento de pasos)
- Sistema lineal de presión: gradiente conjugado (CG) con arranque en caliente

---

## 2. Estructura de archivos

```
NavierStokesSolver/
├── CMakeLists.txt                        # Sistema de construcción CMake (C++17)
├── run.sh                                # Script de lanzamiento rápido
│
├── input/
│   ├── data.txt                          # Geometría del dominio y tamaño de malla
│   ├── stretch.txt                       # Parámetros de estiramiento (beta)
│   ├── sim_params.txt                    # Re, U_lid, dt, rho, max_steps, ss_tol
│   ├── boundaries_lid_cavity.txt         # Condiciones de contorno (WALL, INLET…)
│   └── ghia_reference.txt                # Datos tabulados de Ghia et al.
│
├── output/                               # Resultados de la simulación
│   ├── mesh_P/U/V.vtk                    # Mallas exportadas (ParaView)
│   ├── pressure.vtk                      # Campo de presión final
│   ├── u_velocity.vtk / v_velocity.vtk   # Campos de velocidad finales
│   ├── centerline_u.csv / _v.csv         # Perfiles de línea central (validación)
│   ├── convergence.csv                   # Errores MMS por suite (verify_operators)
│   ├── validation.png                    # Figura de validación vs Ghia
│   └── convergence_rates.png             # Figura de tasas de convergencia
│
├── src/
│   ├── main.cpp                          # Bucle principal de simulación
│   ├── mainUtils.hpp                     # Inspector interactivo de malla (debug)
│   ├── MomentumSolver.hpp/.cpp           # Paso predictor (AB2, difusión + convección)
│   ├── PressurePoissonSystem.hpp/.cpp    # Paso corrector (PPE, proyección)
│   │
│   ├── mesh/
│   │   ├── Mesh_config.h                 # Estructuras GridData, InternalCell, BoundaryFace
│   │   └── MeshGenerator.cpp             # Generación y conectividad de mallas P/U/V
│   │
│   ├── operators/                        # Operadores discretos independientes
│   │   ├── LaplaceOperator.hpp/.cpp      # FVM Laplaciano ∇²p  (P → P)
│   │   ├── GradientOperator.hpp/.cpp     # Gradiente discreto ∇p  (P → U,V)
│   │   ├── DivergenceOperator.hpp/.cpp   # Divergencia discreta ∇·u  (U,V → P)
│   │   └── DiffusionOperator.hpp/.cpp    # Difusión viscosa ν∇²u  (U,V → U,V)
│   │
│   ├── solver/sparse/
│   │   ├── solver_class_sparse.hpp       # Formato CSR, alias Vector, declaraciones
│   │   ├── solver_sparse.cpp             # Jacobi, GS, SOR, CG, PCG
│   │   └── main_solver_sparse.cpp        # Utilidades auxiliares del solucionador
│   │
│   └── verification/
│       └── verify_operators.cpp          # Suite MMS de 10 tests (ejecutable: operator_verify)
│
├── plotters/
│   ├── plot_results.py                   # Figura de validación vs Ghia et al.
│   └── plot_convergence.py               # Figura de tasas de convergencia (MMS)
│
└── docs/
    ├── ALGORITMO.md                      # Este documento
    ├── OPERATORS_REPORT.md               # Informe detallado de operadores y tests
    └── Lid_driven_structured.cpp         # Implementación de referencia
```

**Ejecutables generados en `build/`:**
- `mesh_gen` — solver principal (fuente: `src/main.cpp` + todo `src/`)
- `operator_verify` — verificación independiente de operadores (fuente: `src/verification/verify_operators.cpp`)

---

## 3. La malla MAC escalonada

En una malla MAC (*Marker-and-Cell*), cada variable física vive en un lugar distinto de la celda:

```
      v(i,j+1)
        |
u(i,j)--P(i,j)--u(i+1,j)
        |
      v(i,j)
```

| Variable | Posición | Índices | Fórmula plana |
|----------|----------|---------|---------------|
| Presión `p` | Centro de celda | i=1..N, j=1..M | `pIdx(i,j) = (j-1)·N + (i-1)` |
| Velocidad `u` | Centro de cara vertical | k=1..N+1, j=1..M | `uIdx(k,j) = (j-1)·(N+1) + (k-1)` |
| Velocidad `v` | Centro de cara horizontal | i=1..N, l=1..M+1 | `vIdx(i,l) = (l-1)·N + (i-1)` |

Esta disposición garantiza el acoplamiento natural entre la divergencia de velocidad y el gradiente de presión, **eliminando las oscilaciones en tablero de ajedrez** (checkerboard) que aparecen en mallas colocadas.

Los nodos de frontera llevan los valores prescritos y **no se actualizan** durante el avance temporal:
- U en k=1 (pared Oeste) y k=N+1 (pared Este)
- V en l=1 (pared Sur) y l=M+1 (pared Norte)

---

## 4. Estructuras de datos principales

Definidas en [src/mesh/Mesh_config.h](../src/mesh/Mesh_config.h).

### `InternalCell` — Celda interna del dominio fluido
```cpp
struct InternalCell {
    int global_id;              // Índice lineal global (base 0)
    double x, y;               // Coordenadas del centroide
    double Vol;                 // Volumen (área 2D) de la celda
    int n_east, n_west,        // IDs de celdas vecinas (-1 si frontera)
        n_north, n_south;
    bool bc_east, bc_west,     // ¿La cara da a la frontera?
         bc_north, bc_south;
};
```

### `BoundaryFace` — Cara de frontera
```cpp
struct BoundaryFace {
    double x, y;               // Centro de la cara
    double Area;               // Área de la cara
    int neighbor_id;           // Celda interna adyacente
    BCType type;               // WALL, WALL_FIXED_VALUE, INLET, OUTLET, SYMMETRY
    double value;              // Valor prescrito (velocidad, presión, etc.)
    std::string boco;          // Nombre del grupo de condición de contorno
};
```

### `GridData` — Malla completa (P, U o V)
```cpp
struct GridData {
    std::vector<double> xvc, yvc;                    // Posiciones de caras (vértices)
    std::vector<double> x, y;                        // Centros de celda (incluye fantasmas)
    std::vector<std::vector<double>> Se, Sw, Sn, Ss; // Áreas de caras E/O/N/S
    std::vector<std::vector<double>> V;              // Volúmenes de celda
    std::vector<int> id_to_i, id_to_j;              // Índice global → (i,j)
    int num_active_nodes;                            // Total de celdas internas
    int N_cells_x, M_cells_y;
    std::vector<InternalCell> cells;
    std::vector<BoundaryFace> bound_west, bound_east,
                               bound_south, bound_north;
};
```

Hay **tres instancias de `GridData`**: `pMesh` (presión), `uMesh` (velocidad u), `vMesh` (velocidad v).

---

## 5. Algoritmo global: método de proyección de Chorin

Por cada paso de tiempo `n → n+1`:

```
┌─────────────────────────────────────────────────────────┐
│  Paso 1 — PREDICTOR (MomentumSolver::advance)           │
│                                                          │
│  u* = u^n + dt · [ν∇²u - (u·∇)u]   (AB2 o Euler)       │
│  v* = v^n + dt · [ν∇²v - (u·∇)v]                       │
│                                                          │
│  La velocidad predicha u* NO satisface ∇·u* = 0         │
└───────────────────────┬─────────────────────────────────┘
                        │
                        ▼
┌─────────────────────────────────────────────────────────┐
│  Paso 2a — RHS DE PRESIÓN (PressurePoissonSystem)        │
│                                                          │
│  b = (ρ/dt) · ∇·u*                                      │
└───────────────────────┬─────────────────────────────────┘
                        │
                        ▼
┌─────────────────────────────────────────────────────────┐
│  Paso 2b — RESOLVER PPE (Gradiente Conjugado)           │
│                                                          │
│  ∇²p^{n+1} = b   →   A · p = b                         │
└───────────────────────┬─────────────────────────────────┘
                        │
                        ▼
┌─────────────────────────────────────────────────────────┐
│  Paso 2c — CORRECCIÓN DE VELOCIDAD (projectVelocity)    │
│                                                          │
│  u^{n+1} = u* - (dt/ρ) · ∂p/∂x                         │
│  v^{n+1} = v* - (dt/ρ) · ∂p/∂y                         │
│                                                          │
│  Ahora ∇·u^{n+1} ≈ 0  (solenoidal)                     │
└─────────────────────────────────────────────────────────┘
```

La clave del método es que **desacopla** la ecuación de momentum del campo de presión: primero avanza la velocidad ignorando la presión, y luego corrige para que la divergencia sea cero.

---

## 6. Generación de malla — MeshGenerator

Archivo: [src/mesh/MeshGenerator.cpp](../src/mesh/MeshGenerator.cpp)

### `generate_mesh()` — Malla de presión (P)
Genera la malla base de celdas de presión de N×M celdas.

1. Lee los parámetros de dominio desde `input/data.txt` y `input/stretch.txt`.
2. Genera coordenadas 1D estiradas con `generate_stretched_coords()` en x e y.
3. Calcula áreas de caras `Se, Sw, Sn, Ss` y volúmenes `V` para cada celda.
4. Crea los índices de celda fantasma en las posiciones 0 y N+1/M+1.

### `generate_u_mesh(pMesh)` — Malla de velocidad U
Crea la malla escalonada para la componente horizontal de velocidad.
- Dimensiones: (N+1) × M nodos
- `xvc_u[k] = pMesh.x[k]` para k=1..N (U cae entre centros P consecutivos en x)
- Los bordes k=1 y k=N+1 coinciden con las paredes del dominio

### `generate_v_mesh(pMesh)` — Malla de velocidad V
Análogo al anterior para la componente vertical.
- Dimensiones: N × (M+1) nodos
- `yvc_v[l] = pMesh.y[l]` para l=1..M

### `generate_stretched_coords(vc, num_cells, length, beta, offset)` — Coordenadas estiradas
Genera una distribución 1D de nodos con refinamiento opcional cerca de las paredes mediante transformación `sinh`.

### `build_connectivity(GridData&)` — Conectividad de vecinos
Para cada celda, identifica sus 4 vecinas (E, O, N, S), levanta la bandera `bc_*` si la vecina es una frontera y rellena `id_to_i`, `id_to_j`.

### `define_boundaries(GridData&, filename)` — Condiciones de contorno
Lee el archivo de condiciones de contorno con formato:
```
NOMBRE   CARA    INICIO  FIN   TIPO              VALOR
lid      NORTH   0       9999  WALL_FIXED_VALUE   1.0
south    SOUTH   0       9999  WALL               0.0
```
Rellena los vectores `bound_west`, `bound_east`, `bound_south`, `bound_north` de `GridData`.

### `export_mesh_to_VTK(GridData, filename)` — Exportación para ParaView
Escribe la malla en formato VTK rectilinear grid. Genera `output/mesh_P.vtk`, `mesh_U.vtk`, `mesh_V.vtk`.

---

## 7. Operadores discretos — src/operators/

Los cuatro operadores que intervienen en el algoritmo se han extraído a **clases independientes** en `src/operators/`. Esto permite:
- Reutilizarlos en el solucionador principal y en la suite de verificación sin duplicar código.
- Testearlos de forma aislada con soluciones manufacturadas (Sección 13).
- Sustituir una implementación sin tocar el resto del solver.

Todos toman referencias `const GridData&` y trabajan con el alias `Vector = std::vector<double>`.

---

### 7.1 LaplaceOperator — ∇²p (P → P)

Archivo: [src/operators/LaplaceOperator.hpp](../src/operators/LaplaceOperator.hpp)

Discretiza el operador de Laplace sobre la malla de presión mediante FVM de 2.º orden.

**Coeficientes por celda (i,j):**
```
a_E = Se[i][j] / (x[i+1] - x[i])    a_W = Sw[i][j] / (x[i] - x[i-1])
a_N = Sn[i][j] / (y[j+1] - y[j])    a_S = Ss[i][j] / (y[j] - y[j-1])
a_P = +(a_E + a_W + a_N + a_S)      (diagonal positiva)
```

Convención positivo-definida: `(A·φ)[k] = Σ a_nb·(φ_P − φ_nb) ≈ −V·∇²φ`

Por tanto `(A·φ)/V ≈ −∇²φ`.

**Condición de contorno:** Neumann homogénea (∂p/∂n = 0) — caras de frontera ignoradas.

**Métodos clave:**
- `apply(p)` → `A·p` (producto matriz-vector, útil en tests MMS)
- `assemble()` → `CSRMatrix` (construcción completa en formato CSR)
- `assembleToTriplets()` → `vector<Triplet>` (para pasar a `pinReference` o a `CSRMatrix::fromTriplets`)
- `pinReference(triplets, b)` → fija p[0]=0 para eliminar la singularidad todo-Neumann

**Usado por:** `PressurePoissonSystem::assembleMatrix()`, `verify_operators.cpp` (Suites 1, 5, 9, 10)

---

### 7.2 GradientOperator — ∇p (P → U, V)

Archivo: [src/operators/GradientOperator.hpp](../src/operators/GradientOperator.hpp)

Discretiza el gradiente de presión de la malla P a las mallas escalonadas U/V mediante diferencias finitas centradas de 2.º orden.

**Nodos U interiores** (k=2..N, j=1..M):
```
(∂p/∂x)|_{k,j} = (p[k,j] - p[k-1,j]) / (x_P[k] - x_P[k-1])
```

**Nodos V interiores** (i=1..N, l=2..M):
```
(∂p/∂y)|_{i,l} = (p[i,l] - p[i,l-1]) / (y_P[l] - y_P[l-1])
```

Los nodos de frontera (k=1, k=N+1, l=1, l=M+1) se ponen a cero — la corrección de velocidad en la pared se impone por condición de contorno.

**Operador puramente geométrico** — no incluye factores dt/ρ; estos se aplican en `projectVelocity`.

**Método:** `apply(p, grad_u, grad_v)`

**Usado por:** `PressurePoissonSystem::projectVelocity()`, `verify_operators.cpp` (Suites 2, 5, 8)

**Propiedad algebraica:** G^T = −D en la norma no ponderada (identidad de integración por partes discreta, verificada en Suite 8).

---

### 7.3 DivergenceOperator — ∇·u (U, V → P)

Archivo: [src/operators/DivergenceOperator.hpp](../src/operators/DivergenceOperator.hpp)

Discretiza la divergencia de velocidad integrando el flujo neto sobre cada celda P mediante FVM.

**Para celda (i,j):**
```
flux  = u_E·Se - u_W·Sw + v_N·Sn - v_S·Ss

donde:  u_E = u[uIdx(i+1, j)],   u_W = u[uIdx(i, j)]
        v_N = v[vIdx(i, j+1)],   v_S = v[vIdx(i, j)]

resultado[pIdx(i,j)] = flux / V[i][j]     (unidades: 1/tiempo)
```

El resultado está dividido por el volumen para poder comparar directamente con ∇·u analítico.

**Operador puramente geométrico** — no incluye factores ρ/dt.

**Método:** `apply(u, v)` → `Vector` de tamaño `pMesh.num_active_nodes`

**Nota de muestreo en frontera:** Los nodos U/V de frontera deben muestrearse en la posición de **cara** (`xvc[k-1]`), no en el centroide, para que la divergencia sea exacta (véase Suite 3).

**Usado por:** `PressurePoissonSystem::computeRHS()`, `verify_operators.cpp` (Suites 3, 5, 8)

---

### 7.4 DiffusionOperator — ν∇²u (U, V → U, V)

Archivo: [src/operators/DiffusionOperator.hpp](../src/operators/DiffusionOperator.hpp)

Discretiza el término viscoso ν∇²u y ν∇²v sobre las mallas escalonadas mediante diferencias finitas de 2.º orden en malla no uniforme.

**Fórmula en nodo U interior (k,j):**
```
d²u/dx² = 2/(dx_e + dx_w) · [(u_E - u_P)/dx_e - (u_P - u_W)/dx_w]
d²u/dy² = 2/(dy_n + dy_s) · [(u_N - u_P)/dy_n - (u_P - u_S)/dy_s]
resultado = ν · (d²u/dx² + d²u/dy²)
```

**Celdas fantasma en paredes (método imagen):**
Para nodos adyacentes a las paredes N/S (componente u) y O/E (componente v):
```
u_ghost = 2·u_wall − u_P      →  (u_ghost + u_P)/2 = u_wall  ✓
dy_ghost = 2·(y_P − y_wall)   →  distancia correcta a la cara de frontera
```
Los valores de pared se pasan al constructor: `u_ghost_north`, `u_ghost_south`, `v_ghost_west`, `v_ghost_east`.

**Constructor:**
```cpp
DiffusionOperator(pMesh, uMesh, vMesh, nu,
                  u_ghost_north, u_ghost_south,
                  v_ghost_west,  v_ghost_east);
```

**Métodos clave:**
- `applyU(u, result)` → ν∇²u en todos los nodos U (k=2..N, j=1..M)
- `applyV(v, result)` → ν∇²v en todos los nodos V (i=1..N, l=2..M)
- `diffuseU(u, k, j)` / `diffuseV(v, i, l)` → núcleos por nodo (usados por MomentumSolver)

**Usado por:** `MomentumSolver` (como miembro `diffOp_` en `uRHS`/`vRHS`), `verify_operators.cpp` (Suites 4, 6, 7)

---

### Resumen de operadores

| Operador | Dominio → Rango | Discretización | Clase | Consumidores |
|----------|----------------|----------------|-------|--------------|
| **L** Laplaciano | P → P | FVM 2.º orden | `LaplaceOperator` | `PressurePoissonSystem`, `verify_operators` |
| **G** Gradiente | P → (U,V) | FD 2.º orden | `GradientOperator` | `PressurePoissonSystem`, `verify_operators` |
| **D** Divergencia | (U,V) → P | FVM 2.º orden | `DivergenceOperator` | `PressurePoissonSystem`, `verify_operators` |
| **D_ν** Difusión | (U,V) → (U,V) | FD 2.º orden + fantasma | `DiffusionOperator` | `MomentumSolver`, `verify_operators` |
| **C(u)** Convección | (U,V) → (U,V) | CDS / Upwind | (inline en MomentumSolver) | `MomentumSolver` |

**Identidades algebraicas** verificadas a precisión máquina (Suites 5, 8 de `verify_operators`):
```
D·(G·φ) + (L·φ)/V = 0    (consistencia del laplaciano con D y G)
φᵀ·(D·u) + (G·φ)ᵀ·u = 0  (adjunto discreto G^T = −D)
```

---

## 8. Paso predictor — MomentumSolver

Archivos: [src/MomentumSolver.hpp](../src/MomentumSolver.hpp) / [src/MomentumSolver.cpp](../src/MomentumSolver.cpp)

Implementa el **paso predictor** de Chorin: avanza las velocidades con difusión y convección, ignorando el gradiente de presión.

### Constructor `MomentumSolver(pMesh, uMesh, vMesh, nu, scheme)`

1. Lee los valores de contorno desde `pMesh.bound_*` mediante `ghostValue()`.
2. Construye el miembro `DiffusionOperator diffOp_` con los valores de pared extraídos:
   - `u_ghost_north_` (lid o no-slip), `u_ghost_south_` (no-slip)
   - `v_ghost_west_`, `v_ghost_east_` (no-slip)
3. El esquema de convección se selecciona en construcción: `CDS` (2.º orden) o `UPWIND` (1.º orden).

---

### `applyBoundaryConditions(u, v)` — Imponer condiciones de contorno

Asigna valores prescritos a los nodos de frontera de los vectores de velocidad:
- Nodos U en k=1 y k=N+1: velocidad de entrada o cero
- Nodos V en l=1 y l=M+1: cero (impermeabilidad)

---

### `uRHS(u, v, k, j)` — Lado derecho de la ecuación de momentum x

Calcula `ν∇²u − (u·∇)u` en el nodo U interior (k,j):

**Difusión** (delegada a `DiffusionOperator::diffuseU`):
```
diffuseU(u, k, j)  →  ν·(d²u/dx² + d²u/dy²)
```

**Convección (CDS):**
```
u_P · (u_E − u_W)/(dx_e + dx_w)   +   v_interp · (u_N − u_S)/(dy_n + dy_s)
```

**Convección (UPWIND):**
```
u_P≥0 ? (u_P − u_W)/dx_w : (u_E − u_P)/dx_e     +   análogo para v
```

**Interpolación de v al nodo U** (bilineal, 4 puntos):
```
v_interp = 0.25·(v[vIdx(k-1,j)] + v[vIdx(k,j)] + v[vIdx(k-1,j+1)] + v[vIdx(k,j+1)])
```

---

### `vRHS(u, v, i, l)` — Lado derecho de la ecuación de momentum y

Simétrico a `uRHS` para nodos V interiores (i,l). Difusión delegada a `diffOp_.diffuseV`. Interpolación de u al nodo V usando los 4 nodos U circundantes.

---

### `advance(u, v, dt)` — Avance temporal Adams-Bashforth 2

A partir del segundo paso, usa AB2 para 2.º orden en tiempo:

```
RHS^n  = [uRHS(k,j)  para todos los nodos interiores]
RHS^{n-1} = guardado del paso anterior en rhs_u_old_

u* = u^n + dt · [(3/2)·RHS^n − (1/2)·RHS^{n-1}]   (AB2)
```

El primer paso utiliza Euler explícito (`rhs_u_old_` vacío). Los incrementos se calculan **todos antes de actualizar** usando siempre los valores de `u^n`.

---

## 9. Paso corrector — PressurePoissonSystem

Archivos: [src/PressurePoissonSystem.hpp](../src/PressurePoissonSystem.hpp) / [src/PressurePoissonSystem.cpp](../src/PressurePoissonSystem.cpp)

Resuelve la **ecuación de Poisson para la presión (PPE)** y proyecta la velocidad al subespacio libre de divergencia.

La clase usa internamente los mismos coeficientes FVM que `LaplaceOperator` y las mismas fórmulas de gradiente/divergencia que `GradientOperator`/`DivergenceOperator`, pero las implementa directamente para mantener la clase autónoma.

---

### `assembleMatrix()` → `vector<Triplet>`

Construye la matriz A de la PPE con los mismos coeficientes FVM que `LaplaceOperator::assembleToTriplets()`:

```
a_E = Se[i][j] / (x[i+1] - x[i]),   a_W = Sw[i][j] / (x[i] - x[i-1])
a_N = Sn[i][j] / (y[j+1] - y[j]),   a_S = Ss[i][j] / (y[j] - y[j-1])
a_P = +(a_E + a_W + a_N + a_S)
```

Caras de frontera → coeficiente 0 (Neumann homogéneo). La matriz resultante es **SPD** → apta para CG.

---

### `pinReferencePressure(triplets, b)` — Fijación de la presión de referencia

Elimina la singularidad del sistema todo-Neumann fijando `p[0] = 0` mediante sustitución simétrica de fila/columna 0. También impone `b[0] = 0`.

---

### `setVelocityFields(u_star, v_star)` — Almacenar velocidad predicha

Guarda internamente los campos `u*` y `v*` para que `computeRHS` y `projectVelocity` los usen.

---

### `computeRHS(dt, rho)` → `Vector`

Para cada celda de presión (i,j):
```
b[k] = (ρ/dt) · [u*_E·Se - u*_W·Sw + v*_N·Sn - v*_S·Ss]
```

Mapeo de índices escalonados:
- Cara Este de P(i,j) ↔ U-nodo (k=i+1, j)
- Cara Oeste de P(i,j) ↔ U-nodo (k=i, j)
- Cara Norte de P(i,j) ↔ V-nodo (i, l=j+1)
- Cara Sur de P(i,j)  ↔ V-nodo (i, l=j)

El término `b[0] = 0` se aplica explícitamente en `main.cpp` para mantener la presión de referencia.

---

### `projectVelocity(p_sol, dt, rho, u_next, v_next)` — Corrección de velocidad

Resta el gradiente de presión discreto a la velocidad predicha:

**Nodos U interiores** (k=2..N, j=1..M):
```
u_next[k,j] = u*[k,j] - (dt/ρ) · (p[k,j] - p[k-1,j]) / (x_P[k] - x_P[k-1])
```

**Nodos V interiores** (i=1..N, l=2..M):
```
v_next[i,l] = v*[i,l] - (dt/ρ) · (p[i,l] - p[i,l-1]) / (y_P[l] - y_P[l-1])
```

Los nodos de frontera no se corrigen.

---

## 10. Solucionadores lineales — solver_sparse

Archivos: [src/solver/sparse/solver_class_sparse.hpp](../src/solver/sparse/solver_class_sparse.hpp) / [solver_sparse.cpp](../src/solver/sparse/solver_sparse.cpp)

### Tipo `Vector` y formato CSR
```cpp
using Vector = std::vector<double>;   // alias global

struct CSRMatrix {
    int n;
    std::vector<double> values;       // valores no nulos
    std::vector<int>    col_indices;  // columna de cada valor
    std::vector<int>    row_ptr;      // puntero de inicio de cada fila (tamaño n+1)
    std::vector<double> diagonal;     // diagonal extraída para O(1) de precondicionador
};
```
Construido desde triplets con `CSRMatrix::fromTriplets(n, triplets)`.

---

### Solucionadores disponibles

| Clase | Método | Notas |
|-------|--------|-------|
| `Jacobi` | Iteración de Jacobi | Lento, paralelizable |
| `GS` | Gauss-Seidel | Convergencia ~2× Jacobi |
| `SOR` | GS + relajación (ω=1.5) | Más rápido que GS puro |
| `CG` | Gradiente Conjugado | **Solucionador por defecto** para PPE SPD |
| `PCG` | CG precondicionado (Jacobi diagonal) | Útil en mallas muy no uniformes |

**CG con arranque en caliente** (`warm-start`):
```cpp
Vector solve(const CSRMatrix& A, const Vector& b, const Vector& x0);
```
La solución del paso anterior (`p^n`) se usa como estimación inicial, reduciendo las iteraciones a ~5–20 en flujo quasi-estacionario.

---

## 11. Bucle principal — main.cpp

Archivo: [src/main.cpp](../src/main.cpp)

### Parámetros leídos de archivos de entrada

| Archivo | Parámetros |
|---------|-----------|
| `input/data.txt` | `L_domain`, `H_domain`, `W_depth`, coordenadas de nodos |
| `input/stretch.txt` | Factores beta de estiramiento en x e y |
| `input/sim_params.txt` | `Re`, `U_lid`, `dt`, `rho`, `max_steps`, `ss_tol` |

La viscosidad cinemática se deriva: `nu = U_lid · L_domain / Re`.

### Inicialización (se ejecuta una sola vez)

```
1. MeshConfig::get_data(), get_stretch(), get_sim_params()
2. generate_mesh()          → pMesh
   generate_u_mesh(pMesh)   → uMesh
   generate_v_mesh(pMesh)   → vMesh
3. define_boundaries(pMesh, "boundaries_lid_cavity.txt")
4. build_connectivity(pMesh), build_connectivity(uMesh), build_connectivity(vMesh)
5. export_mesh_to_VTK()  ×3  → output/mesh_P/U/V.vtk

6. PressurePoissonSystem::assembleMatrix()   → triplets
   pinReferencePressure(triplets, dummy_b)
   CSRMatrix::fromTriplets()                 → A (CSR, ensamblada una sola vez)

7. Inicializar u=v=p=0
8. MomentumSolver::applyBoundaryConditions(u, v)
```

### Bucle temporal

```
Para step = 0, 1, 2, ... hasta max_steps o convergencia:

  A) PREDICTOR
     u_old = u_star;  v_old = v_star
     momentum.advance(u_star, v_star, dt)           → u*, v* (AB2)

  B) RHS DE PRESIÓN
     poisson.setVelocityFields(u_star, v_star)
     b_rhs = poisson.computeRHS(dt, rho)            → b = (ρ/dt)·∇·u*
     b_rhs[0] = 0                                   → referencia de presión

  C) RESOLVER PPE
     p_sol = CG_solver.solve(A, b_rhs, p_sol)       → arranque en caliente

  D) CORRECCIÓN DE VELOCIDAD
     poisson.projectVelocity(p_sol, dt, rho, u_star, v_star)
                                                    → u^{n+1} = u* - (dt/ρ)·∇p

  E) REIMPOSICIÓN DE CONTORNO
     momentum.applyBoundaryConditions(u_star, v_star)

  F) CHEQUEO DE ESTADO ESTACIONARIO (cada 100 pasos)
     si ||u^{n+1} - u^n||∞ / dt < ss_tol → CONVERGE, salir
```

### Parámetros de simulación por defecto

| Parámetro | Valor típico | Descripción |
|-----------|-------------|-------------|
| `Re` | 400 | Número de Reynolds |
| `U_lid` | 1.0 | Velocidad de la tapa (m/s) |
| `nu` | U·L/Re | Viscosidad cinemática |
| `dt` | 1e-3 | Paso temporal (s) |
| `rho` | 1.0 | Densidad (kg/m³) |
| `max_steps` | 50 000 | Máximo de pasos temporales |
| `ss_tol` | 1e-6 | Tolerancia de estado estacionario |

---

## 12. Concatenación de funciones (grafo de llamadas)

```
main()
│
├── Malla
│   ├── get_data() / get_stretch() / get_sim_params()
│   ├── generate_mesh() → pMesh
│   │   └── generate_stretched_coords() ×2
│   ├── generate_u_mesh(pMesh) → uMesh
│   ├── generate_v_mesh(pMesh) → vMesh
│   ├── define_boundaries(pMesh/uMesh/vMesh, ...)
│   ├── build_connectivity(pMesh/uMesh/vMesh)
│   └── export_mesh_to_VTK() ×3
│
├── Sistema lineal de presión (una vez)
│   ├── PressurePoissonSystem::assembleMatrix()   → triplets
│   │   └── [mismos coef. FVM que LaplaceOperator]
│   ├── pinReferencePressure(triplets, dummy_b)
│   └── CSRMatrix::fromTriplets()                → A (CSR)
│
├── MomentumSolver::MomentumSolver(pMesh, uMesh, vMesh, nu)
│   └── DiffusionOperator::DiffusionOperator(...)  ← miembro diffOp_
│
├── MomentumSolver::applyBoundaryConditions(u, v)   ← init
│
└── BUCLE TEMPORAL ────────────────────────────────────────────
    │
    ├── MomentumSolver::advance(u, v, dt)
    │   ├── [j=1..M, k=2..N]  uRHS(u, v, k, j)
    │   │   ├── diffOp_.diffuseU(u, k, j)        ← DiffusionOperator
    │   │   │   ├── d²u/dx²  (FD 2.º orden)
    │   │   │   └── d²u/dy²  (celda fantasma en j=1 y j=M)
    │   │   └── Convección CDS/Upwind  +  v interpolado bilineal
    │   ├── [l=2..M, i=1..N]  vRHS(u, v, i, l)
    │   │   ├── diffOp_.diffuseV(v, i, l)        ← DiffusionOperator
    │   │   └── Convección  +  u interpolado bilineal
    │   └── u += dt·(3/2·RHS^n − 1/2·RHS^{n-1});  v += … (AB2)
    │
    ├── PressurePoissonSystem::setVelocityFields(u*, v*)
    │
    ├── PressurePoissonSystem::computeRHS(dt, rho)
    │   └── [k=0..N_P-1]  b[k] = (ρ/dt)·(u*_E·Se - u*_W·Sw + v*_N·Sn - v*_S·Ss)
    │
    ├── b_rhs[0] = 0
    │
    ├── CG::solve(A, b_rhs, p_anterior)           ← arranque en caliente
    │   ├── r = b - A·p_ant
    │   └── [iter] α, x, r, β, p  (convergencia ||r||<1e-6·||b||)
    │
    ├── PressurePoissonSystem::projectVelocity(p_sol, dt, rho, u*, v*)
    │   ├── [k=2..N, j=1..M]  u*[k,j] -= (dt/ρ)·(p[k,j]-p[k-1,j])/(x[k]-x[k-1])
    │   └── [i=1..N, l=2..M]  v*[i,l] -= (dt/ρ)·(p[i,l]-p[i,l-1])/(y[l]-y[l-1])
    │
    ├── MomentumSolver::applyBoundaryConditions(u*, v*)
    │
    └── [cada 100 pasos]  ||u*−u_old||∞/dt < ss_tol → BREAK
```

---

## 13. Verificación de operadores

Ejecutable: `./build/operator_verify`  
Fuente: [src/verification/verify_operators.cpp](../src/verification/verify_operators.cpp)  
Documentación detallada: [docs/OPERATORS_REPORT.md](OPERATORS_REPORT.md)

El programa de verificación usa el **Método de Soluciones Manufacturadas (MMS)**: se elige una función analítica suave, se muestrea en la malla, se aplica el operador discreto y se compara con el resultado analítico. Se repite en mallas uniformes N ∈ {8, 16, 32, 64} y se estima la tasa de convergencia `r = log₂(e_N / e_{2N})`.

Al terminar escribe `output/convergence.csv` con todos los errores por suite y N.

### Suites implementadas

| Suite | Operador | Función manufacturada | Resultado esperado |
|-------|----------|-----------------------|--------------------|
| **1** | `LaplaceOperator` | φ=cos(πx)cos(πy), Neumann | tasa ~2.0 |
| **2** | `GradientOperator` | φ=cos(πx)cos(πy), nodos interiores | tasa ~2.0 (∂p/∂x y ∂p/∂y) |
| **3** | `DivergenceOperator` | u=sin(2πx)cos(πy), muestreo en cara | tasa ~2.0 |
| **4** | `DiffusionOperator::applyU` | u=sin(πy) (cte. en x), ν=1 | tasa ~2.0 |
| **5** | Identidad D·G + L/V = 0 | φ=cos(πx)cos(πy) | error ~1e-14..1e-12 |
| **6** | `DiffusionOperator::applyV` | v=sin(πx) (cte. en y), ν=1 | tasa ~2.0 |
| **7** | Difusión con BC no nula | u=y (lineal), u_N=1, u_S=0 | error = 0 exacto |
| **8** | Adjunto G^T = −D | vectores aleatorios interiores | residuo ~1e-11 (prec. máquina) |
| **9** | SPD de LaplaceOperator | vectores aleatorios | simetría ~1e-14, φᵀAφ > 0 |
| **10** | Laplaciano en malla tanh (γ=2) | φ=cos(πx)cos(πy) | tasa ~1.0 (*) |

(*) La Suite 10 expone una **limitación del estencil actual**: en mallas no uniformes (`h_e ≠ h_w`), el término cruzado `φ'·(h_e−h_w)/dx_i = O(1/N)` degrada la precisión a primer orden. El fix consiste en sustituir `1/dx_i` por `2/(h_e+h_w)` en los coeficientes FVM.

### Figura de tasas de convergencia

```bash
./build/operator_verify          # genera output/convergence.csv
python3 plotters/plot_convergence.py   # genera output/convergence_rates.png
```

La figura muestra 6 paneles:
- (0,0) Laplaciano uniforme — pendiente ~2.0
- (0,1) Gradiente ∂p/∂x y ∂p/∂y — pendiente ~2.0
- (0,2) Divergencia — pendiente ~2.0
- (1,0) Difusión U y V — pendiente ~2.0
- (1,1) Residuos de identidades algebraicas (Suites 5, 7–9) — a precisión máquina
- (1,2) Laplaciano uniforme vs estirado — pendiente 2.0 vs 1.0

---

## 14. Postproceso y validación

### Exportación de campos a VTK

Al finalizar el bucle temporal, `main.cpp` exporta tres ficheros mediante la función lambda `writeFieldVTK`:

```
output/pressure.vtk     — Campo escalar de presión (malla P)
output/u_velocity.vtk   — Componente u (malla U)
output/v_velocity.vtk   — Componente v (malla V)
```

Formato: VTK Rectilinear Grid ASCII, visualizable en **ParaView**.

### Extracción de perfiles de línea central

```
output/centerline_u.csv  — u(y) en x = L/2  (columna U más cercana al centro)
output/centerline_v.csv  — v(x) en y = H/2  (fila V más cercana al centro)
```

Columnas: `y, y_norm, u, u_norm` y `x, x_norm, v, v_norm`.

### Figura de validación

```bash
python3 plotters/plot_results.py   # genera output/validation.png
```

La figura (2×3 paneles) muestra:
- Campos de presión, u y v en 2D (contornos de color)
- Comparación cuantitativa de u(y) y v(x) frente a Ghia et al. (1982)
- Líneas de corriente del flujo

El script lee `input/sim_params.txt` para seleccionar automáticamente la columna correcta de la tabla de Ghia según el Re simulado (disponibles: 100, 400, 1000, 3200, 5000, 7500, 10000).

**Precisión esperada:**
- Malla 15×15: error ~2–5% en perfiles de velocidad
- Malla 65×65: error < 1%

---

## 15. Parámetros numéricos y estabilidad

### Condición CFL (estabilidad temporal para difusión explícita)
```
dt ≤ C · dx² / ν
```
Con dx ≈ 1/15 ≈ 0.067 m y ν = U·L/Re:
- Re=100: ν=0.01 → dt_máx ≈ 3×10⁻⁵ s (difusión pura)
- dt = 1×10⁻³ s usado → dentro del margen estable con CFL ≈ 0.9

### Convergencia del sistema lineal (PPE)
- Tolerancia relativa: `||r|| < 1e-6 · ||b||`
- CG con arranque en caliente: ~5–20 iteraciones por paso
- Número de condición: κ(A) ~ O(N²) — bien condicionado para mallas pequeñas

### Convergencia a estado estacionario
- Criterio: `||u^{n+1} - u^n||∞ / dt < ss_tol` (por defecto 1e-6)
- Re=100, malla 15×15: ~500–1000 pasos temporales
- Re=400, malla 65×65: ~10 000–30 000 pasos

### Integración temporal AB2 vs Euler
- El primer paso usa Euler (sin historia previa) — ligero transitorio inicial
- A partir del segundo paso, AB2 proporciona 2.º orden en tiempo
- El estado `rhs_u_old_` / `rhs_v_old_` en `MomentumSolver` guarda el RHS del paso anterior

---

*Documentación generada para el proyecto NavierStokesSolver.*
