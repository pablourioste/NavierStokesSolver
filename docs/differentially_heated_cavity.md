# Differentially Heated Cavity (DHC) — Documentación de la implementación

Esta guía explica, de principio a fin, cómo está implementada la **cavidad
diferencialmente calentada** (DHC, *Differentially Heated Cavity*) en el solver,
qué ecuaciones resuelve, y **cómo viaja la información** desde los archivos de
entrada → la malla y su discretización (incluido el *stretching*) → los
operadores → los solvers de momento, presión y energía.

El código reutiliza por completo el solver de *Lid-Driven Cavity* (cavidad con
tapa móvil). La DHC es un **modo** que se activa con un único flag y que añade
dos piezas modulares: la ecuación de la energía ([EnergySolver.cpp](../src/EnergySolver.cpp))
y el término de flotabilidad de Boussinesq ([Buoyancy.hpp](../src/Buoyancy.hpp)).
**Ningún archivo del núcleo (momento, presión, operadores) se modificó** para
soportar la DHC.

---

## 1. Las ecuaciones físicas

La DHC resuelve las ecuaciones de Navier–Stokes incompresibles **acopladas con
la ecuación de la energía** bajo la **aproximación de Boussinesq**, en la
adimensionalización **Ra–Pr** (de Vahl Davis, 1983).

### 1.1 Conservación de masa (continuidad)

$$
\nabla \cdot \mathbf{u} = 0
$$

Se impone de forma exacta mediante el **método de proyección** (paso de presión).

### 1.2 Momento (con flotabilidad de Boussinesq)

$$
\frac{\partial \mathbf{u}}{\partial t} + (\mathbf{u}\cdot\nabla)\mathbf{u}
= -\nabla p + \nu\,\nabla^2 \mathbf{u} + \mathbf{f}_b
$$

En el escalado Ra–Pr, la **difusión de momento** vale `nu = Pr` y la fuerza de
flotabilidad sólo actúa en la dirección vertical:

$$
\mathbf{f}_b = (0,\; Ra\,Pr\,\theta)
$$

donde `θ` es la temperatura adimensional. Esto se ve directamente en el código:

- En [main.cpp:196](../src/main.cpp#L196): `const double nu = energy_on ? Pr : ...`
- En [Buoyancy.hpp:34](../src/Buoyancy.hpp#L34): `const double coeff = dt * Ra * Pr;`

### 1.3 Energía (advección–difusión de la temperatura)

$$
\frac{\partial \theta}{\partial t} + (\mathbf{u}\cdot\nabla)\theta
= \alpha\,\nabla^2 \theta
$$

Con la difusividad térmica `α = 1` en el escalado Ra–Pr
([main.cpp:221](../src/main.cpp#L221), se pasa `alpha=1.0` al `EnergySolver`).
La ecuación implementada literalmente en el encabezado de
[EnergySolver.hpp:14](../src/EnergySolver.hpp#L14):

```
d(theta)/dt + (u . grad) theta = alpha * laplacian(theta)
```

### 1.4 Números adimensionales

| Símbolo | Significado | Dónde entra |
|---------|-------------|-------------|
| `Ra` (Rayleigh) | intensidad de la convección natural | barrido en `rayleigh_cases.txt` → `coeff = dt·Ra·Pr` |
| `Pr` (Prandtl)  | relación difusividad de momento / térmica | `sim_params.txt` → `nu = Pr` |
| `α = 1`         | difusividad térmica adimensional | fijada en `main.cpp` |

---

## 2. Selección del problema: cómo se activa la DHC

El conmutador entre **Lid-Driven Cavity** y **Differentially Heated Cavity** es
un único entero en [input/sim_params.txt](../input/sim_params.txt):

```text
solve_energy = 1     # 0 -> Lid cavity (momento + presión)
                     # 1 -> Differentially Heated Cavity (+ energía + flotabilidad)
Pr = 0.71            # Número de Prandtl (sólo se usa si solve_energy = 1)
```

Recorrido del flag:

1. [`MeshConfig::get_sim_params`](../src/mesh/MeshGenerator.cpp#L686) lee el
   archivo y `update_parameter` asigna `solve_energy` y `Pr` a los miembros de
   la clase ([MeshGenerator.cpp:732](../src/mesh/MeshGenerator.cpp#L732)).
2. En [main.cpp:129](../src/main.cpp#L129):
   `const bool energy_on = (myConfig.solve_energy != 0);`
3. Si `energy_on`:
   - Se leen los Rayleigh a barrer desde `rayleigh_cases.txt`
     ([main.cpp:135](../src/main.cpp#L135)).
   - Se cargan las condiciones de contorno **de temperatura** desde
     `boundaries_T.txt` ([main.cpp:136](../src/main.cpp#L136)).
   - El bucle de casos itera sobre **Rayleigh** en vez de Reynolds.

---

## 3. Los archivos de entrada (input) y qué controla cada uno

Todo en `input/`. El `main()` los lee al arrancar en este orden
([main.cpp:82-93](../src/main.cpp#L82-L93)):

| Archivo | Lo lee | Controla |
|---------|--------|----------|
| [data.txt](../input/data.txt) | `get_data` | Geometría: `L_domain`, `H_domain`, profundidad `W_depth`, y los nodos NW/NE/SW/SE que dividen el dominio en zonas |
| [stretch.txt](../input/stretch.txt) | `get_stretch` | Número de celdas por zona y los **β de estiramiento** (refinamiento) |
| [sim_params.txt](../input/sim_params.txt) | `get_sim_params` | `dt`, `rho`, `U_lid`, `nu`, `max_steps`, `ss_tol`, **`solve_energy`**, **`Pr`** |
| [boundaries.txt](../input/boundaries.txt) | `define_boundaries` | Condiciones de contorno de **velocidad** |
| [boundaries_T.txt](../input/boundaries_T.txt) | `define_temperature_boundaries` | Condiciones de contorno de **temperatura** (sólo DHC) |
| [rayleigh_cases.txt](../input/rayleigh_cases.txt) | `read_rayleigh_cases` | Lista de números de Rayleigh a barrer |

### 3.1 `data.txt` — geometría

```text
L_domain = 1.0    # ancho del dominio
H_domain = 1.0    # alto del dominio
W_depth  = 0.1    # profundidad (2D real, pero las áreas se multiplican por W_depth)
# nodos NW/NE/SW/SE = 1.0  -> dominio de una sola zona (sin obstáculo)
```

Para la cavidad cuadrada estándar, los cuatro nodos se colocan de modo que la
**zona "west/south" cubre todo el dominio** y las zonas "center"/"east" tienen
longitud cero (con `N_cells = 0` no se copian a la malla).

### 3.2 `stretch.txt` — discretización y refinamiento

```text
beta_x_west   = 0.8     # estiramiento en X
beta_y_south  = 0.8     # estiramiento en Y
beta_x_center = 0.0     # 0.0 -> malla uniforme en esa zona
N_cells_west  = 50      # 50 celdas en X
M_cells_west  = 50      # 50 celdas en Y
N_cells_center = 0      # zonas vacías
```

El número total de celdas es la suma de las tres zonas:
`N_cells = N_cells_west + N_cells_center + N_cells_east`
([MeshGenerator.cpp:158](../src/mesh/MeshGenerator.cpp#L158)).

### 3.3 `boundaries_T.txt` — contornos térmicos (clave de la DHC)

```text
# BOCO_NAME  FACE   IDX_START IDX_END  TYPE             VALUE
hot_west     WEST   0   9999   WALL_ISOTHERMAL   1.0    # pared caliente θ = 1
cold_east    EAST   0   9999   WALL_ISOTHERMAL   0.0    # pared fría     θ = 0
top          NORTH  0   9999   WALL_ADIABATIC    0.0    # flujo nulo
bottom       SOUTH  0   9999   WALL_ADIABATIC    0.0    # flujo nulo
```

Esta es la configuración física que **define** la cavidad diferencialmente
calentada: pared oeste caliente, pared este fría, tapas superior e inferior
adiabáticas. El gradiente horizontal de temperatura genera la flotabilidad que
mueve el fluido.

---

## 4. La malla: generación, *staggered grid* y stretching

### 4.1 Malla escalonada (MAC / staggered)

Se generan **tres mallas** a partir de la malla de presión
([main.cpp:86-88](../src/main.cpp#L86-L88)):

```
pMesh = generate_mesh();           // presión y temperatura θ (centro de celda)
uMesh = generate_u_mesh(pMesh);    // velocidad u  (caras verticales,  N+1 en X)
vMesh = generate_v_mesh(pMesh);    // velocidad v  (caras horizontales, M+1 en Y)
```

- **p y θ** viven en el **centro** de cada celda.
- **u** vive en las caras verticales (la malla U tiene `N+1` columnas y comparte
  las caras Y de P, [MeshGenerator.cpp:343](../src/mesh/MeshGenerator.cpp#L343)).
- **v** vive en las caras horizontales (la malla V tiene `M+1` filas y comparte
  las caras X de P, [MeshGenerator.cpp:400](../src/mesh/MeshGenerator.cpp#L400)).

Esta disposición es la razón por la que `θ` se almacena igual que la presión: las
velocidades escalonadas **son directamente** las velocidades normales a las caras
de cada celda P, sin necesidad de interpolar para el flujo convectivo (ver §6.2).

### 4.2 El stretching (refinamiento)

El refinamiento se hace con una función `sinh` en
[`generate_stretched_coords`](../src/mesh/MeshGenerator.cpp#L79). Idea:

1. Si `|β| < 1e-6` → **malla uniforme** (`Δx = L/num_cells`).
2. Si `β > 0` → se refina **al final** del segmento (`x = L`).
3. Si `β < 0` → se refina **al principio** (`x = 0`).

El algoritmo:

```
xi      = -1 + 2*i/num_cells          # coordenada lógica uniforme [-1,1]
g_xi    = (sinh(β(xi - xi_focus)) - val_min) / val_range * 2 - 1   # estirada
x       = (L/2)*(1 + g_xi) + offset   # mapeo a coordenada física
```

Cuanto mayor `|β|`, más concentradas quedan las celdas cerca del foco. Para la
DHC esto permite **refinar cerca de las paredes caliente/fría**, donde las capas
límite térmicas son finas y donde se mide el número de Nusselt.

El `offset` posiciona cada zona (west/center/east) en su lugar dentro del
dominio global, y luego los segmentos se **fusionan** en `data.xvc` / `data.yvc`
([MeshGenerator.cpp:208-237](../src/mesh/MeshGenerator.cpp#L208-L237)).

### 4.3 De coordenadas a la geometría que usa el solver

Tras generar las caras (`xvc`, `yvc`), `generate_mesh` calcula y guarda en
`GridData` todo lo que los operadores necesitan
([MeshGenerator.cpp:239-261](../src/mesh/MeshGenerator.cpp#L239-L261)):

- `x[i]`, `y[j]` — **centroides** de celda (incluye celdas fantasma 0 y N+1).
- `Se, Sw, Sn, Ss` — **áreas de las caras** (= `Δy·W_depth` o `Δx·W_depth`).
- `V` — **volumen** de celda (`Δx·Δy·W_depth`).

Estas estructuras son **el contrato** entre la malla y los solvers: el
`EnergySolver` lee `pMesh.Se/Sw/Sn/Ss/V` y los centroides directamente para
construir sus flujos y su laplaciano (ver §6).

### 4.4 Conectividad y contornos de temperatura

- [`build_connectivity`](../src/mesh/MeshGenerator.cpp#L448) crea la tabla de
  vecinos (E/W/N/S) y marca qué celdas tocan una frontera.
- [`define_temperature_boundaries`](../src/mesh/MeshGenerator.cpp#L993) crea las
  listas paralelas `bound_*_T` (oeste/este/norte/sur). **Copia la geometría** de
  las listas de velocidad `bound_*` y luego sobreescribe sólo el `type` y el
  `value` con lo que diga `boundaries_T.txt`. Las caras no mencionadas quedan
  **adiabáticas por defecto**. Esto deja intactas las BC de velocidad.

---

## 5. El bucle temporal: cómo se alimenta el solver paso a paso

El corazón está en [main.cpp:225-258](../src/main.cpp#L225-L258). Cada `caseVal`
del barrido es un Rayleigh; dentro, el bucle de tiempo es:

```cpp
for (int step = 0; step < max_steps; ++step) {
    u_old = u_star;

    // (1) PREDICTOR de momento (AB2): u* = u^n + dt(3/2 R^n - 1/2 R^{n-1})
    momentum.advance(u_star, v_star, dt);

    // (2) FLOTABILIDAD: inyecta dt·Ra·Pr·θ en v* (sólo DHC)
    if (energy_on)
        addBuoyancy(T, v_star, dt, Ra, Pr, pMesh, vMesh);

    // (3) PRESIÓN: monta el RHS de Poisson a partir de la divergencia de u*
    poisson.setVelocityFields(u_star, v_star);
    b_rhs = poisson.computeRHS(dt, rho);
    b_rhs[0] = 0.0;                       // fija la presión de referencia (pinning)

    // (4) Resuelve el sistema de Poisson con Gradiente Conjugado (warm-start)
    p_sol = CG_solver.solve(A, b_rhs, p_sol);

    // (5) PROYECCIÓN: corrige u* -> campo de velocidad sin divergencia
    poisson.projectVelocity(p_sol, dt, rho, u_star, v_star);
    momentum.applyBoundaryConditions(u_star, v_star);

    // (6) ENERGÍA: transporta θ con el campo ya corregido (sólo DHC)
    if (energy_on)
        energy.advance(T, u_star, v_star, dt);

    // (7) Chequeo de estado estacionario cada 100 pasos
}
```

**El orden importa.** La flotabilidad se inyecta en `v*` **después** del
predictor de momento y **antes** de la proyección de presión: así la proyección
hace que el campo final sea incompresible **incluyendo** el forzado boyante. Y la
energía se avanza **al final**, usando el campo de velocidad ya corregido y libre
de divergencia (`u_star`, `v_star`), lo cual mantiene consistente el transporte
de θ. Esta es precisamente la razón por la que **no hubo que tocar** el
`MomentumSolver` ni el `PressurePoissonSystem`: la DHC se acopla por inyección
externa.

### 5.1 Condición inicial

- `u_star, v_star, p_sol = 0` ([main.cpp:211-213](../src/main.cpp#L211-L213)).
- `T = 0.5` en todas las celdas (valor intermedio entre la pared caliente θ=1 y
  la fría θ=0) ([main.cpp:223](../src/main.cpp#L223)).

---

## 6. El EnergySolver en detalle

Archivo: [src/EnergySolver.cpp](../src/EnergySolver.cpp). Resuelve
`dθ/dt + (u·∇)θ = α∇²θ` sobre la **malla de presión**. Su estructura imita a la
del `MomentumSolver`.

### 6.1 Construcción: lectura de las BC de temperatura

En el constructor ([EnergySolver.cpp:17-34](../src/EnergySolver.cpp#L17-L34)),
`extractWallT` lee de `pMesh.bound_*_T` el tipo y valor de cada pared
(`bc_west_`, `Tw_west_`, …). Por defecto, lista vacía → adiabática.

### 6.2 El término convectivo (volúmenes finitos conservativos)

En [`tRHS`](../src/EnergySolver.cpp#L48), las velocidades normales a las caras de
la celda P son **directamente** las velocidades escalonadas (sin interpolar):

```cpp
const double u_e = u[uIdx(i+1, j)];  // cara este
const double u_w = u[uIdx(i,   j)];  // cara oeste
const double v_n = v[vIdx(i, j+1)];  // cara norte
const double v_s = v[vIdx(i, j  )];  // cara sur
```

Los flujos se calculan con **UPWIND** de primer orden por defecto (recomendado
para estabilidad a Rayleigh alto) o **CDS** (diferencias centradas). Con UPWIND,
el valor en cada cara se toma de la celda *aguas arriba*
([EnergySolver.cpp:100-106](../src/EnergySolver.cpp#L100-L106)):

```cpp
Fe = u_e * ((u_e >= 0.0) ? T_P : T_E) * Se;   // y análogamente Fw, Fn, Fs
conv = (Fe - Fw + Fn - Fs) / Vol;             // divergencia por unidad de volumen
```

Como la velocidad **normal a las paredes es cero** (no penetración), los flujos
en las caras de contorno se anulan y la convección no necesita ningún ghost de
temperatura.

### 6.3 El término difusivo (laplaciano con ghost-cell image)

Diferencias finitas de 2º orden sobre la malla no uniforme. En las paredes se usa
el **método de la celda imagen** ([EnergySolver.cpp:59-78](../src/EnergySolver.cpp#L59-L78)):

- **Isoterma (Dirichlet):** `θ_ghost = 2·θ_pared − θ_P`
- **Adiabática (Neumann):** `θ_ghost = θ_P` (gradiente normal nulo)

```cpp
d2T_dx2 = 2/(dx_e+dx_w) * ((T_E - T_P)/dx_e - (T_P - T_W)/dx_w);
d2T_dy2 = 2/(dy_n+dy_s) * ((T_N - T_P)/dy_n - (T_P - T_S)/dy_s);
diff    = alpha * (d2T_dx2 + d2T_dy2);
```

El RHS total es `return diff - conv;`.

### 6.4 Integración temporal (Adams-Bashforth 2)

En [`advance`](../src/EnergySolver.cpp#L118): se calcula **todo** el RHS desde
`θ^n` en un buffer (estilo Euler diferido, para que actualizar una celda no
contamine a sus vecinas dentro del paso), y luego:

- **Primer paso:** Euler — `θ += dt·rhs`.
- **Resto:** AB2 — `θ += dt·(1.5·rhs − 0.5·rhs_old)`.

> **Nota — ¿por qué la temperatura usa los mismos pesos que el momento, pero la
> presión no tiene pesos?**
>
> La integración temporal de `θ` es **idéntica, línea por línea**, a la del
> `MomentumSolver`: mismos pesos AB2 (`3/2`, `−1/2`), mismo arranque con Euler y
> el mismo buffer diferido. El `EnergySolver` se diseñó a propósito como un
> espejo del `MomentumSolver` para el escalar `θ`. Compárese
> [EnergySolver.cpp:128-134](../src/EnergySolver.cpp#L128-L134) con
> [MomentumSolver.cpp:289-300](../src/MomentumSolver.cpp#L289-L300).
>
> La **presión es un caso distinto**: no es una ecuación de evolución (no tiene
> `∂p/∂t`), sino una **restricción elíptica**. En cada paso se resuelve la
> ecuación de Poisson `∇²p = (ρ/dt)·∇·u*` con el Gradiente Conjugado
> ([main.cpp:236](../src/main.cpp#L236)) para imponer `∇·u = 0`. Por eso **no
> tiene pesos AB2 ni `rhs_old`**: es un *solve* implícito que da la presión
> instantánea, no una marcha temporal.
>
> | Campo | Tipo de ecuación | Avance temporal |
> |-------|------------------|-----------------|
> | `u, v` (momento) | evolución (parabólica) | AB2 explícito `3/2, −1/2` |
> | `θ` (energía)    | evolución (parabólica) | AB2 explícito `3/2, −1/2` ← igual que momento |
> | `p` (presión)    | restricción (elíptica) | no se integra; Poisson con CG cada paso |

---

## 7. El acoplamiento de flotabilidad (Buoyancy.hpp)

[`addBuoyancy`](../src/Buoyancy.hpp#L24) es header-only y se aplica como un
incremento explícito externo sobre `v*`:

```cpp
const double coeff = dt * Ra * Pr;
for (int l = 2; l <= M; ++l)          // sólo nodos V interiores
    for (int i = 1; i <= N; ++i) {
        // θ vive en celdas P; el nodo V (i,l) está entre P(i,l-1) y P(i,l)
        const double theta_v = 0.5 * (T[pIdx(i, l-1)] + T[pIdx(i, l)]);
        v[vIdx(i, l)] += coeff * theta_v;
    }
```

Detalles importantes:

- Como `θ` está en P y `v` en las caras horizontales, **se interpola θ
  linealmente** al nodo V (promedio de las dos celdas P adyacentes).
- Sólo se tocan los nodos V **interiores** (`l = 2..M`); los nodos de pared
  (no penetración) se dejan intactos.
- Sólo añade fuerza **vertical** (es la columna `v`), coherente con
  `f_b = (0, Ra·Pr·θ)`.

---

## 8. Validación y salida

Cuando `energy_on`, al final de cada caso ([main.cpp:406-461](../src/main.cpp#L406-L461))
se calculan las métricas del **benchmark de de Vahl Davis (1983)**:

- `u_max` en el plano vertical medio (`x = L/2`).
- `v_max` en el plano horizontal medio (`y = H/2`).
- **Número de Nusselt medio** en la pared caliente (oeste):

$$
\overline{Nu} = \frac{1}{H}\sum_j \frac{\theta_{pared} - \theta_P}{x_P - x_{pared}}\,\Delta y_j
$$

```cpp
Nu_avg += (Twall - T[k]) / dist * dy;   // flujo local * altura de celda
Nu_avg /= H_domain;
```

Salidas por caso (en `output/Ra_<valor>/`):

- `temperature.vtk`, `u_velocity.vtk`, `v_velocity.vtk`, `pressure.vtk` (campos en
  formato VTK rectilíneo, [main.cpp:266-270](../src/main.cpp#L266-L270)).
- `dhc_metrics.csv` con `Ra, Pr, u_max, v_max, Nu_avg`.

---

## 9. Aviso CRÍTICO de estabilidad

Tanto la energía como el momento son **totalmente explícitos**, así que la DHC
necesita un `dt` **mucho menor** que la cavidad con tapa. El límite de difusión
explícita es:

$$
dt \lesssim \frac{\Delta x^2}{2\,\max(\alpha,\nu)}, \qquad \alpha = 1
$$

En una malla de producción 150² esto implica `dt ~ 1e-5` o menos. Usar
`dt = 0.001` (el default del lid cavity) **diverge a NaN** (el CG falla en
converger primero). El `sim_params.txt` actual ya está puesto en `dt = 0.00001`
para la DHC.

> Verificado a Ra=1e4 en malla 40² con `dt=2e-5`: `u_max≈16.7`, `v_max≈20.1`,
> `Nu≈2.31` — coincide con de Vahl Davis (1983) dentro de un pequeño porcentaje.

---

## 10. Resumen del flujo de información

```
data.txt ─────────┐
stretch.txt ──────┤→ MeshConfig ──→ generate_mesh ──→ pMesh (xvc,yvc,x,y,Se..Ss,V)
sim_params.txt ───┘                      │                 ├─→ generate_u_mesh → uMesh
   (solve_energy, Pr, dt, rho)           │                 └─→ generate_v_mesh → vMesh
                                         │
boundaries.txt ───→ define_boundaries ──→ bound_*  (BC velocidad, ghosts)
boundaries_T.txt ─→ define_temperature_boundaries → bound_*_T (BC temperatura)
rayleigh_cases.txt → read_rayleigh_cases → lista de Ra

                         ┌──────────────── BUCLE TEMPORAL (por cada Ra) ───────────────┐
 pMesh,uMesh,vMesh,nu=Pr │ 1. momentum.advance(u*,v*)          (ν∇²u − (u·∇)u, AB2)    │
 ───────────────────────→│ 2. addBuoyancy → v* += dt·Ra·Pr·θ                            │
 T (=θ), Ra, Pr          │ 3-5. Poisson(CG) + projectVelocity  (∇·u = 0)               │
 ───────────────────────→│ 6. energy.advance(θ, u*, v*)        (α∇²θ − (u·∇)θ, AB2)    │
                         └──────────────────────────────────────────────────────────────┘
                                         │
                                         └─→ output/Ra_*/{temperature,u,v,pressure}.vtk
                                                        dhc_metrics.csv (u_max,v_max,Nu)
```

---

### Archivos clave

| Pieza | Archivo |
|-------|---------|
| Bucle principal y selección de modo | [src/main.cpp](../src/main.cpp) |
| Generación de malla + stretching | [src/mesh/MeshGenerator.cpp](../src/mesh/MeshGenerator.cpp) |
| Estructuras de malla y config | [src/mesh/Mesh_config.h](../src/mesh/Mesh_config.h) |
| Ecuación de la energía (θ) | [src/EnergySolver.cpp](../src/EnergySolver.cpp) · [.hpp](../src/EnergySolver.hpp) |
| Flotabilidad de Boussinesq | [src/Buoyancy.hpp](../src/Buoyancy.hpp) |
| Predictor de momento | [src/MomentumSolver.cpp](../src/MomentumSolver.cpp) |
| Presión / proyección | [src/PressurePoissonSystem.cpp](../src/PressurePoissonSystem.cpp) |
</content>
</invoke>
