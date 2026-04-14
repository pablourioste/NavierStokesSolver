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
7. [Paso predictor — MomentumSolver](#7-paso-predictor--momentumsolver)
8. [Paso corrector — PressurePoissonSystem](#8-paso-corrector--pressurepoissonssystem)
9. [Solucionadores lineales — solver_sparse](#9-solucionadores-lineales--solver_sparse)
10. [Bucle principal — main.cpp](#10-bucle-principal--maincpp)
11. [Concatenación de funciones (grafo de llamadas)](#11-concatenación-de-funciones-grafo-de-llamadas)
12. [Postproceso y validación](#12-postproceso-y-validación)
13. [Parámetros numéricos y estabilidad](#13-parámetros-numéricos-y-estabilidad)

---

## 1. Visión general

El solver resuelve las ecuaciones de Navier-Stokes para flujo incompresible en 2D:

```
∂u/∂t + (u·∇)u = -∇p/ρ + ν∇²u
∇·u = 0
```

**Caso de prueba:** Cavidad cuadrada con tapa deslizante (*lid-driven cavity*), Re = 100.
- Dominio: L = H = 1.0 m
- Tapa superior (norte): velocidad prescrita u = 1.0 m/s
- Resto de paredes: no-slip (u = v = 0)
- Validación: datos de referencia de Ghia et al. (1982)

**Método numérico clave:**
- Discretización espacial: diferencias finitas de 2.º orden sobre malla MAC escalonada
- Integración temporal: Euler explícito
- Acoplamiento presión-velocidad: proyección de Chorin (dos pasos por instante de tiempo)
- Sistema lineal de presión: gradiente conjugado (CG) con arranque en caliente

---

## 2. Estructura de archivos

```
NavierStokesSolver/
└── Mesh_creation/src/
    ├── main.cpp                        # Bucle principal de simulación
    ├── mainUtils.hpp                   # Utilidad interactiva de inspección de malla
    ├── MomentumSolver.hpp/.cpp         # Paso predictor (momentum)
    ├── PressurePoissonSystem.hpp/.cpp  # Paso corrector (presión + proyección)
    ├── mesh/
    │   ├── Mesh_config.h               # Estructuras de datos (celdas, caras, malla)
    │   └── MeshGenerator.cpp           # Generación de malla con estiramiento
    └── solver/sparse/
        ├── solver_class_sparse.hpp     # Formato CSR y declaraciones de solucionadores
        └── solver_sparse.cpp           # Jacobi, Gauss-Seidel, SOR, CG, PCG
```

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

| Variable | Posición | Por qué |
|----------|----------|---------|
| Presión `p` | Centro de celda | Punto de equilibrio del balance de fuerzas |
| Velocidad `u` (horizontal) | Centro de cara vertical | Flujo neto en x a través de la cara |
| Velocidad `v` (vertical) | Centro de cara horizontal | Flujo neto en y a través de la cara |

Esta disposición garantiza el acoplamiento natural entre la divergencia de velocidad y el gradiente de presión, **eliminando las oscilaciones en tablero de ajedrez** (checkerboard) que aparecen en mallas colocadas.

---

## 4. Estructuras de datos principales

Definidas en [Mesh_config.h](Mesh_creation/src/mesh/Mesh_config.h).

### `InternalCell` — Celda interna del dominio fluido
```cpp
struct InternalCell {
    int global_id;              // Índice lineal global (base 0)
    double x, y;               // Coordenadas del centroide
    double Vol;                 // Volumen (área 2D) de la celda
    int n_east, n_west,        // IDs de celdas vecinas
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
    double value;              // Valor prescrito (velocidad, temperatura, etc.)
    std::string boco;          // Nombre del grupo de condición de contorno
};
```

### `GridData` — Malla completa (P, U o V)
```cpp
struct GridData {
    std::vector<double> xvc, yvc;          // Posiciones de caras (vértices)
    std::vector<double> x, y;              // Centros de celda
    std::vector<std::vector<double>> Se, Sw, Sn, Ss;  // Áreas de caras E/O/N/S
    std::vector<std::vector<double>> V;    // Volúmenes de celda
    std::vector<std::vector<int>> connectivity;        // [id][dir] → id vecino
    std::vector<int> id_to_i, id_to_j;    // Índice global → (i,j)
    int num_active_nodes;                  // Total de celdas internas
    std::vector<InternalCell> cells;
    std::vector<BoundaryFace> bound_west, bound_east,
                               bound_south, bound_north;
    int N_cells_x, M_cells_y;
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
│  u* = u^n + dt · [ν∇²u - (u·∇)u]                       │
│  v* = v^n + dt · [ν∇²v - (u·∇)v]                       │
│                                                          │
│  La velocidad predicha u* NO satisface ∇·u* = 0         │
└───────────────────────┬─────────────────────────────────┘
                        │
                        ▼
┌─────────────────────────────────────────────────────────┐
│  Paso 2a — RHS DE PRESIÓN (PressurePoissonSystem)        │
│                                                          │
│  b = -(ρ/dt) · ∇·u*                                     │
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

Archivo: [MeshGenerator.cpp](Mesh_creation/src/mesh/MeshGenerator.cpp)

### `generate_mesh()` — Malla de presión (P)
Genera la malla base de celdas de presión de N×M celdas.

1. Lee los parámetros de dominio desde `input/data.txt` y `input/stretch.txt`.
2. Genera coordenadas 1D estiradas con `generate_stretched_coords()` en x e y.
3. Calcula áreas de caras `Se, Sw, Sn, Ss` para cada celda (necesarias en la integración FVM).
4. Añade celdas fantasma para gestión de condiciones de contorno.

### `generate_u_mesh(pMesh)` — Malla de velocidad U
Crea la malla escalonada para la componente horizontal de velocidad.
- Dimensiones: (N+1) × M nodos
- Los nodos U se sitúan **entre los centros de celdas P** en la dirección x
- Los bordes oeste (k=1) y este (k=N+1) coinciden con las paredes del dominio

### `generate_v_mesh(pMesh)` — Malla de velocidad V
Análogo al anterior para la componente vertical.
- Dimensiones: N × (M+1) nodos
- Los nodos V se sitúan entre centros de celdas P en la dirección y

### `generate_stretched_coords(vc, num_cells, length, beta, offset)` — Coordenadas estiradas
Genera una distribución 1D de nodos con refinamiento opcional cerca de las paredes.

- **`beta = 0`:** Malla uniforme (espaciado regular)
- **`beta > 0`:** Refinamiento en x = L (extremo derecho/superior)
- **`beta < 0`:** Refinamiento en x = 0 (extremo izquierdo/inferior)

Transformación: ξ ∈ [-1, 1] → x físico mediante `sinh(β·(ξ - ξ_focus))`.

### `build_connectivity(GridData&)` — Conectividad de vecinos
Para cada celda (i,j), identifica sus 4 vecinas (E, O, N, S) y levanta la bandera `bc_*` si la vecina es una frontera. Rellena los mapas inversos `id_to_i`, `id_to_j`.

### `define_boundaries(GridData&, filename)` — Condiciones de contorno
Lee el archivo `input/boundaries_lid_cavity.txt` con formato:
```
BOCO_NAME   FACE    START   END    TYPE              VALUE
lid         NORTH   0       9999   WALL_FIXED_VALUE  1.0
wall_south  SOUTH   0       9999   WALL              0.0
```
Rellena los vectores `bound_west`, `bound_east`, `bound_south`, `bound_north` de `GridData`.

### `export_mesh_to_VTK(GridData, filename)` — Exportación para ParaView
Escribe la malla en formato VTK rectilinear grid con los IDs de celda como dato escalar. Permite visualizar la estructura de la malla en ParaView.

---

## 7. Paso predictor — MomentumSolver

Archivos: [MomentumSolver.hpp](Mesh_creation/src/MomentumSolver.hpp) / [MomentumSolver.cpp](Mesh_creation/src/MomentumSolver.cpp)

Implementa el **paso predictor** de Chorin: avanza las velocidades con difusión y convección, ignorando el gradiente de presión.

### Indexación lógica (base 1)

| Campo | Índices | Fórmula de índice lineal |
|-------|---------|--------------------------|
| U-nodos | k ∈ [1, N+1], j ∈ [1, M] | `uIdx(k,j) = (j-1)·(N+1) + (k-1)` |
| V-nodos | i ∈ [1, N], l ∈ [1, M+1] | `vIdx(i,l) = (l-1)·N + (i-1)` |

Los nodos de frontera (k=1, k=N+1 para U; l=1, l=M+1 para V) llevan los valores prescritos de contorno y **no se actualizan** durante el paso temporal.

---

### Constructor `MomentumSolver(pMesh, uMesh, vMesh, nu, scheme)`
Extrae los valores de las celdas fantasma desde las caras de frontera de `pMesh.bound_*`:
- `u_ghost_north_`, `u_ghost_south_`: velocidad fantasma para las paredes N/S (usadas en difusión de u)
- `v_ghost_west_`, `v_ghost_east_`: velocidad fantasma para las paredes E/O (usadas en difusión de v)

---

### `applyBoundaryConditions(u, v)` — Imponer condiciones de contorno
Asigna valores prescritos a los nodos de frontera de los vectores de velocidad:
- Nodos U en k=1 y k=N+1: velocidad de entrada o cero (no-slip)
- Nodos V en l=1 y l=M+1: siempre cero (impermeabilidad)

---

### `uRHS(u, v, k, j)` — Lado derecho de la ecuación de momentum x

Calcula `ν∇²u - (u·∇)u` en el nodo U interior (k,j) para k ∈ [2, N].

**Término de difusión** (diferencias centradas en malla no uniforme):
```
∂²u/∂x² = 2/(dx_e + dx_w) · [(u_E - u_P)/dx_e - (u_P - u_W)/dx_w]
∂²u/∂y² = 2/(dy_n + dy_s) · [(u_N - u_P)/dy_n - (u_P - u_S)/dy_s]
```

**Celdas fantasma en la frontera** (método imagen — garantiza el valor correcto en la pared):
```
u_ghost = 2·u_wall - u_P     →     (u_ghost + u_P)/2 = u_wall  ✓
dy      = 2·(y_P - y_wall)   →     distancia correcta a la cara de frontera
```

**Término convectivo (interpolación de v al nodo U):**
La velocidad v se interpola bilinearmente desde los 4 nodos V circundantes:
```
v_interp = 0.25 · (v[vIdx(k-1,j)] + v[vIdx(k,j)] + v[vIdx(k-1,j+1)] + v[vIdx(k,j+1)])
```

**Esquemas de convección disponibles:**
- `CDS` (diferencias centradas): 2.º orden, menor difusión numérica, menos estable
- `UPWIND` (primer orden): más difusivo, más estable con celdas gruesas

---

### `vRHS(u, v, i, l)` — Lado derecho de la ecuación de momentum y

Simétrico a `uRHS` para los nodos V interiores (i,l) con l ∈ [2, M].
La interpolación de u al nodo V usa los 4 nodos U circundantes.

---

### `advance(u, v, dt)` — Avance temporal explícito

```
1. Calcular du[k,j] = uRHS(u, v, k, j)  para todos los nodos U interiores
2. Calcular dv[i,l] = vRHS(u, v, i, l)  para todos los nodos V interiores
3. u = u + dt·du
4. v = v + dt·dv
   (Los nodos de frontera no se modifican, su delta = 0)
```

Los incrementos se calculan **todos antes de actualizar**, usando siempre los valores del instante `n`. Esto evita la propagación de errores dentro del mismo paso (coherencia con Euler explícito).

---

## 8. Paso corrector — PressurePoissonSystem

Archivos: [PressurePoissonSystem.hpp](Mesh_creation/src/PressurePoissonSystem.hpp) / [PressurePoissonSystem.cpp](Mesh_creation/src/PressurePoissonSystem.cpp)

Resuelve la **ecuación de Poisson para la presión (PPE)** y luego proyecta la velocidad al subespacio libre de divergencia.

---

### `assembleMatrix()` — Ensamblado de la matriz Laplaciana

Construye la matriz A de la forma `A·p = b` mediante el método de volúmenes finitos (FVM).

Para cada celda de presión k con índices lógicos (i,j):

```
Coeficiente este:   a_E = Se[i][j] / (x[i+1] - x[i])
Coeficiente oeste:  a_W = Sw[i][j] / (x[i]   - x[i-1])
Coeficiente norte:  a_N = Sn[i][j] / (y[j+1] - y[j])
Coeficiente sur:    a_S = Ss[i][j] / (y[j]   - y[j-1])
Diagonal:           a_P = a_E + a_W + a_N + a_S   (suma de vecinos activos)
```

- Las caras de **frontera** aportan coeficiente 0 (Neumann homogéneo: ∂p/∂n = 0).
- Solo los **vecinos internos** contribuyen a la matriz.
- El resultado es una matriz **simétrica definida positiva** (SPD) — idónea para CG.

Devuelve un vector de `Triplet(fila, columna, valor)` que luego se convierte a formato CSR.

---

### `pinReferencePressure(triplets, b)` — Fijación de la presión de referencia

Con condiciones de contorno puramente de Neumann (∂p/∂n = 0 en todas las paredes), la PPE tiene infinitas soluciones (la presión está definida salvo una constante). Para eliminar la singularidad:

- **Fijar p[0] = 0** mediante eliminación simétrica de la fila/columna 0:
  - Fila 0: todos los fuera-diagonal a cero, diagonal = 1
  - Columna 0: todos los fuera-diagonal a cero
  - b[0] = 0.0

Esto hace el sistema singular → no singular sin alterar la solución física (solo fija el nivel de referencia de presión).

---

### `setVelocityFields(u_star, v_star)` — Almacenar velocidad predicha

Guarda internamente los campos `u*` y `v*` procedentes del paso predictor.
Valida que los tamaños coincidan con `uMesh` y `vMesh`.

---

### `computeRHS(dt, rho)` — Vector del lado derecho (divergencia de u*)

Para cada celda de presión k con índices (i,j), calcula la integral de volumen de ∇·u*:

```
∫∇·u* dV ≈ u*_E · Se - u*_W · Sw + v*_N · Sn - v*_S · Ss
```

donde los índices de velocidad corresponden a las caras escalonadas:
- Cara este de P(i,j)  ↔  U-nodo (k=i+1, j)
- Cara oeste de P(i,j) ↔  U-nodo (k=i,   j)
- Cara norte de P(i,j) ↔  V-nodo (i,   l=j+1)
- Cara sur de P(i,j)   ↔  V-nodo (i,   l=j)

El RHS de la PPE es:
```
b[k] = -(ρ/dt) · div(u*)[k]
```
El signo negativo es coherente con la convención Laplaciana positivo-definida.

---

### `projectVelocity(p_sol, dt, rho, u_next, v_next)` — Corrección de velocidad

Resta el gradiente de presión discreto a la velocidad predicha:

**Nodos U interiores** (i ∈ [2, N], j ∈ [1, M]):
```
u_next[i,j] = u*[i,j] - (dt/ρ) · (p[i,j] - p[i-1,j]) / (x[i] - x[i-1])
```

**Nodos V interiores** (l ∈ [2, M], i ∈ [1, N]):
```
v_next[i,l] = v*[i,l] - (dt/ρ) · (p[i,l] - p[i,l-1]) / (y[l] - y[l-1])
```

Los nodos de frontera **no se corrigen** (retienen el valor prescrito de la condición de contorno).

---

## 9. Solucionadores lineales — solver_sparse

Archivos: [solver_class_sparse.hpp](Mesh_creation/src/solver/sparse/solver_class_sparse.hpp) / [solver_sparse.cpp](Mesh_creation/src/solver/sparse/solver_sparse.cpp)

### Formato de almacenamiento: CSR (Compressed Sparse Row)
```cpp
struct CSRMatrix {
    int n;                        // Tamaño del sistema
    std::vector<double> values;   // Valores no nulos
    std::vector<int> col_indices; // Columnas de cada valor no nulo
    std::vector<int> row_ptr;     // Punteros de inicio de cada fila (tamaño n+1)
    std::vector<double> diagonal; // Diagonal extraída para acceso O(1)
};
```
Construido desde triplets con `CSRMatrix::fromTriplets(n, triplets)`.

---

### Solucionadores implementados

#### Jacobi
Iteración: `x_i^{n+1} = (b_i - Σ_{j≠i} A_ij · x_j^n) / A_ii`
- Lento pero trivialmente paralelizable.

#### Gauss-Seidel (GS)
Iteración: usa los valores más recientes disponibles en cada fila.
```
x_i^{n+1} = (b_i - Σ_{j<i} A_ij · x_j^{n+1} - Σ_{j>i} A_ij · x_j^n) / A_ii
```
- Convergencia más rápida que Jacobi.
- Criterio: `||b - Ax|| < 1e-6 · ||b||`, máx. 10n iteraciones.

#### SOR (Successive Over-Relaxation)
GS con relajación: `x_i^{n+1} = (1-ω)·x_i^n + (ω/A_ii)·(b_i - resto)`
- Factor de relajación: ω = 1.25 (fijo).
- Más rápido que GS en sistemas con autovalores sesgados.

#### Gradiente Conjugado (CG) — **solucionador principal**
Algoritmo estándar para matrices SPD. Incluye sobrecarga con **arranque en caliente**:
```cpp
Vector solve(const CSRMatrix& A, const Vector& b, const Vector& x0)
```
El arranque en caliente (`x0 = p^n`) reutiliza la solución del paso temporal anterior como estimación inicial, reduciendo el número de iteraciones en flujos quasi-estacionarios.

Pasos del algoritmo:
```
r = b - A·x0
p = r
loop:
  Ap = A·p
  α  = (r·r) / (p·Ap)         ← longitud de paso óptima
  x  = x + α·p                ← actualizar solución
  r  = r - α·Ap               ← actualizar residuo
  β  = (r_new·r_new)/(r·r)    ← dirección de búsqueda nueva
  p  = r + β·p
  converge si ||r|| < 1e-6·||b||
```

#### PCG (Gradiente Conjugado Precondicionado)
CG con precondicionador diagonal (Jacobi): M = diag(√A_ii).
Reduce el número de condición efectivo, útil en mallas muy no uniformes.

---

## 10. Bucle principal — main.cpp

Archivo: [main.cpp](Mesh_creation/src/main.cpp)

### Inicialización (se ejecuta una sola vez)

```
1. Leer parámetros de dominio  →  MeshConfig::get_data("input/data.txt")
2. Leer parámetros de malla    →  MeshConfig::get_stretch("input/stretch.txt")
3. Generar malla P             →  generate_mesh()
4. Generar malla U             →  generate_u_mesh(pMesh)
5. Generar malla V             →  generate_v_mesh(pMesh)
6. Leer condiciones de contorno → define_boundaries(...)
7. Construir conectividad      →  build_connectivity()  [×3, una por malla]
8. Exportar mallas a VTK       →  export_mesh_to_VTK()  [×3]

9. Ensamblar matriz Laplaciana →  PressurePoissonSystem::assembleMatrix()
10. Fijar presión de referencia →  pinReferencePressure(triplets, b)
11. Construir matriz CSR        →  CSRMatrix::fromTriplets(n, triplets)

12. Inicializar u = 0, v = 0
13. Imponer condiciones de contorno iniciales
```

### Bucle temporal (hasta 50 000 pasos o convergencia)

```
Para cada paso n = 0, 1, 2, ...:

  A) PREDICTOR
     MomentumSolver::advance(u, v, dt)
       → calcula u*, v*

  B) RHS DE PRESIÓN
     PressurePoissonSystem::setVelocityFields(u*, v*)
     PressurePoissonSystem::computeRHS(dt, rho)
       → calcula b = -(ρ/dt)·∇·u*
     b[0] = 0   ← imponer presión de referencia

  C) RESOLVER PPE
     CG::solve(A, b, p_anterior)  ← arranque en caliente
       → calcula p^{n+1}

  D) CORRECCIÓN DE VELOCIDAD
     PressurePoissonSystem::projectVelocity(p^{n+1}, dt, rho, u^{n+1}, v^{n+1})
       → u^{n+1} = u* - (dt/ρ)·∇p
       → v^{n+1} = v* - (dt/ρ)·∇p

  E) REIMPOSICIÓN DE CONTORNO
     MomentumSolver::applyBoundaryConditions(u^{n+1}, v^{n+1})

  F) CHEQUEO DE ESTADO ESTACIONARIO (cada 100 pasos)
     Si ||u^{n+1} - u^n||∞ / dt < 1e-6  →  CONVERGE, salir del bucle
```

### Parámetros de simulación

| Parámetro | Valor | Descripción |
|-----------|-------|-------------|
| `dt` | 1e-3 | Paso temporal |
| `rho` | 1.0 | Densidad del fluido |
| `nu` | 0.01 | Viscosidad cinemática (Re = U·L/ν = 100) |
| `max_steps` | 50 000 | Máximo de iteraciones temporales |
| `ss_tol` | 1e-6 | Tolerancia de estado estacionario |
| Malla | 15×15 | Celdas uniformes |

---

## 11. Concatenación de funciones (grafo de llamadas)

```
main()
│
├── Inicialización de malla
│   ├── get_data()
│   ├── get_stretch()
│   ├── generate_mesh()               → pMesh
│   │   └── generate_stretched_coords()
│   ├── generate_u_mesh(pMesh)        → uMesh
│   ├── generate_v_mesh(pMesh)        → vMesh
│   ├── define_boundaries(pMesh, ...) + define_boundaries(uMesh, ...) + define_boundaries(vMesh, ...)
│   ├── build_connectivity(pMesh)     + build_connectivity(uMesh) + build_connectivity(vMesh)
│   └── export_mesh_to_VTK(...)       ×3
│
├── Ensamblado del sistema lineal de presión
│   ├── PressurePoissonSystem::assembleMatrix()   → triplets
│   ├── PressurePoissonSystem::pinReferencePressure(triplets, b)
│   └── CSRMatrix::fromTriplets(n, triplets)      → A (CSR)
│
├── MomentumSolver::MomentumSolver(pMesh, uMesh, vMesh, nu, CDS)
│   └── [extrae celdas fantasma de pMesh.bound_*]
│
├── MomentumSolver::applyBoundaryConditions(u, v)  ← inicialización
│
└── BUCLE TEMPORAL ─────────────────────────────────────────────────
    │
    ├── MomentumSolver::advance(u, v, dt)
    │   ├── [j=1..M, k=2..N]  uRHS(u, v, k, j)
    │   │   ├── Difusión x:  ∂²u/∂x²  (FD 2.º orden, malla no uniforme)
    │   │   ├── Difusión y:  ∂²u/∂y²  (con celda fantasma si j=1 o j=M)
    │   │   ├── Convección:  u·∂u/∂x  (CDS o Upwind)
    │   │   └── Convección:  v·∂u/∂y  (v interpolado bilinealmente)
    │   ├── [l=2..M, i=1..N]  vRHS(u, v, i, l)
    │   │   └── [análogo, con u interpolado]
    │   └── u += dt·du;  v += dt·dv
    │
    ├── PressurePoissonSystem::setVelocityFields(u*, v*)
    │
    ├── PressurePoissonSystem::computeRHS(dt, rho)
    │   └── [k=0..N_P-1]  b[k] = -(ρ/dt)·(u*_E·Se - u*_W·Sw + v*_N·Sn - v*_S·Ss)
    │
    ├── b[0] = 0.0
    │
    ├── CG::solve(A, b, p_anterior)
    │   ├── r = b - A·p_anterior
    │   ├── [iter=1..10n]
    │   │   ├── multiply(A, p) → Ap
    │   │   ├── α = dot(r,r)/dot(p,Ap)
    │   │   ├── x += α·p;  r -= α·Ap
    │   │   ├── β = dot(r_new,r_new)/dot(r,r)
    │   │   └── p = r + β·p
    │   └── return p_sol
    │
    ├── PressurePoissonSystem::projectVelocity(p_sol, dt, rho, u_next, v_next)
    │   ├── [j=1..M, i=2..N]  u_next[i,j] -= (dt/ρ)·(p[i,j]-p[i-1,j])/(x[i]-x[i-1])
    │   └── [l=2..M, i=1..N]  v_next[i,l] -= (dt/ρ)·(p[i,l]-p[i,l-1])/(y[l]-y[l-1])
    │
    ├── MomentumSolver::applyBoundaryConditions(u_next, v_next)
    │
    └── [cada 100 pasos]  ||u_next - u_old||∞ / dt < 1e-6  →  BREAK
```

---

## 12. Postproceso y validación

### Exportación de campos a VTK
Al finalizar el bucle, se exportan tres ficheros:
- `pressure.vtk` — Campo escalar de presión sobre la malla P
- `u_velocity.vtk` — Componente u sobre la malla U
- `v_velocity.vtk` — Componente v sobre la malla V

Visualizables directamente en **ParaView**.

### Validación con Ghia et al. (1982)
El código extrae perfiles de velocidad a lo largo de las líneas centrales:

- **Perfil u(y)** a lo largo de x = L/2 (columna de U-nodos más cercana al centro)
- **Perfil v(x)** a lo largo de y = H/2 (fila de V-nodos más cercana al centro)

Escribe `centerline_u.csv` y `centerline_v.csv` e imprime los datos de referencia de Ghia para comparación directa.

Precisión esperada:
- Malla 15×15: error ~2-5% en velocidades
- Malla 65×65: error < 1%

---

## 13. Parámetros numéricos y estabilidad

### Condición CFL (estabilidad temporal)
Para el paso difusivo explícito:
```
dt ≤ C · dx² / ν
```
Con dx ≈ 1/15 ≈ 0.067 m y ν = 0.01:
- dt_máx ≈ 3×10⁻⁵ s (difusión)
- dt = 1×10⁻³ s usado → CFL ≈ 0.9 (dentro del margen estable)

### Convergencia del sistema lineal
- Tolerancia relativa: `||r|| < 1e-6 · ||b||`
- CG con arranque en caliente: ~5-20 iteraciones por paso en estado quasi-estacionario
- Número de condición de A (malla 15×15): κ ≈ N² ≈ 225 → bien condicionado

### Convergencia a estado estacionario
- Criterio: `||u^{n+1} - u^n||∞ / dt < 1e-6`
- Para Re = 100 en malla 15×15: convergencia típica en ~500-1000 pasos temporales
- Máximo permitido: 50 000 pasos (margen de seguridad amplio)

---

*Documentación generada para el proyecto NavierStokesSolver — Cavidad con tapa deslizante, Re = 100.*
