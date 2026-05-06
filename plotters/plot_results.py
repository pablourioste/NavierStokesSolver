import os
import re
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from scipy.interpolate import RegularGridInterpolator

BASE_DIR = os.path.dirname(os.path.abspath(__file__))
PROJECT_ROOT = os.path.dirname(BASE_DIR)
OUTPUT_DIR = os.path.join(PROJECT_ROOT, "output")

# =============================================================================
# Datos de referencia — Ghia et al. (1982), todas las Re
# =============================================================================
# TABLE I: col 0 = y/H, cols 1-7 = u/U para Re 100,400,1000,3200,5000,7500,10000
ghia_u_table = np.array([
    [1.0000,  1.00000,  1.00000,  1.00000,  1.00000,  1.00000,  1.00000,  1.00000],
    [0.9766,  0.84123,  0.75837,  0.65928,  0.53236,  0.48223,  0.47244,  0.47221],
    [0.9688,  0.78871,  0.68439,  0.57492,  0.48296,  0.46120,  0.47048,  0.47783],
    [0.9609,  0.73722,  0.61756,  0.51117,  0.46547,  0.45992,  0.47323,  0.48070],
    [0.9531,  0.68717,  0.55892,  0.46604,  0.45863,  0.46036,  0.47048,  0.47804],
    [0.8516,  0.23151,  0.29093,  0.33304,  0.33556,  0.33754,  0.34228,  0.34635],
    [0.7344,  0.00332,  0.16256,  0.18719,  0.20087,  0.20580,  0.20591,  0.20673],
    [0.6172, -0.13641,  0.02135,  0.05702,  0.08344,  0.08183,  0.08342,  0.08344],
    [0.5000, -0.20581, -0.11477, -0.06080, -0.03111, -0.03039, -0.03432, -0.03581],
    [0.4531, -0.21090, -0.17119, -0.10648, -0.07540, -0.07404, -0.07503, -0.07496],
    [0.2813, -0.15662, -0.32726, -0.27805, -0.23186, -0.22855, -0.23176, -0.23186],
    [0.1719, -0.10150, -0.24299, -0.38289, -0.32709, -0.33050, -0.32393, -0.32709],
    [0.1016, -0.06434, -0.14612, -0.29730, -0.38000, -0.37633, -0.38324, -0.38000],
    [0.0703, -0.04775, -0.10338, -0.22220, -0.41657, -0.41810, -0.43643, -0.41657],
    [0.0625, -0.04192, -0.09266, -0.20196, -0.42768, -0.43590, -0.43048, -0.42768],
    [0.0547, -0.03717, -0.08186, -0.18109, -0.43025, -0.43025, -0.41487, -0.43025],
    [0.0000,  0.00000,  0.00000,  0.00000,  0.00000,  0.00000,  0.00000,  0.00000],
])
# TABLE II: col 0 = x/L, cols 1-7 = v/U para Re 100,400,1000,3200,5000,7500,10000
ghia_v_table = np.array([
    [1.0000,  0.00000,  0.00000,  0.00000,  0.00000,  0.00000,  0.00000,  0.00000],
    [0.9688, -0.05906, -0.12146, -0.21388, -0.39017, -0.49774, -0.53858, -0.54302],
    [0.9609, -0.07391, -0.15663, -0.27669, -0.47425, -0.55069, -0.55216, -0.52987],
    [0.9531, -0.08864, -0.19254, -0.33714, -0.52357, -0.55408, -0.52347, -0.49099],
    [0.9453, -0.10313, -0.22847, -0.39188, -0.54053, -0.52876, -0.48590, -0.45863],
    [0.9063, -0.16914, -0.23827, -0.51500, -0.44307, -0.41442, -0.41050, -0.41496],
    [0.8594, -0.22445, -0.44993, -0.42665, -0.37401, -0.36214, -0.36213, -0.36737],
    [0.8047, -0.24533, -0.38598, -0.31966, -0.31184, -0.30018, -0.30448, -0.30719],
    [0.5000,  0.05454,  0.05186,  0.02526,  0.00999,  0.00945,  0.00824,  0.00831],
    [0.2344,  0.17527,  0.30174,  0.32235,  0.28188,  0.27280,  0.27348,  0.27224],
    [0.2266,  0.17507,  0.30203,  0.33075,  0.29030,  0.28066,  0.28117,  0.28003],
    [0.1563,  0.16077,  0.28124,  0.37095,  0.37119,  0.35368,  0.35060,  0.35070],
    [0.0938,  0.12317,  0.22965,  0.32627,  0.42768,  0.42951,  0.41824,  0.41487],
    [0.0781,  0.10890,  0.20920,  0.30353,  0.41906,  0.43648,  0.43564,  0.43124],
    [0.0703,  0.10091,  0.19713,  0.29012,  0.40917,  0.43329,  0.44030,  0.43733],
    [0.0625,  0.09233,  0.18360,  0.27485,  0.39560,  0.42447,  0.43979,  0.43983],
    [0.0000,  0.00000,  0.00000,  0.00000,  0.00000,  0.00000,  0.00000,  0.00000],
])
GHIA_RE_LIST = [100, 400, 1000, 3200, 5000, 7500, 10000]


def ghia_column_for_re(re_value):
    return min(range(len(GHIA_RE_LIST)), key=lambda i: abs(GHIA_RE_LIST[i] - re_value))


# =============================================================================
# Leer VTK rectilinear grid (formato ASCII simple generado por el solver)
# =============================================================================
def read_vtk_rectilinear(path):
    """Devuelve (x_coords, y_coords, field_2d) del VTK RECTILINEAR_GRID."""
    if not os.path.exists(path):
        raise FileNotFoundError(f"❌ No se encontró: {path}")
        
    with open(path) as f:
        lines = f.readlines()

    x_coords = y_coords = None
    values = []
    i = 0
    while i < len(lines):
        line = lines[i].strip()
        if line.startswith("X_COORDINATES"):
            n = int(line.split()[1])
            data = []
            i += 1
            while len(data) < n:
                data += list(map(float, lines[i].split()))
                i += 1
            x_coords = np.array(data)
            continue
        if line.startswith("Y_COORDINATES"):
            n = int(line.split()[1])
            data = []
            i += 1
            while len(data) < n:
                data += list(map(float, lines[i].split()))
                i += 1
            y_coords = np.array(data)
            continue
        if line.startswith("LOOKUP_TABLE"):
            i += 1
            while i < len(lines) and not lines[i].strip().startswith(("SCALARS", "VECTORS")):
                tok = lines[i].split()
                values += list(map(float, tok))
                i += 1
            continue
        i += 1

    N = len(x_coords) - 1
    M = len(y_coords) - 1
    
    # Manejo seguro del reshape
    if len(values) < (N * M):
        print(f"Advertencia en {os.path.basename(path)}: Faltan valores para llenar la malla.")
        field = np.zeros((M, N)) # Evita que el script pete
    elif len(values) > (N * M) and len(values) % (N * M) == 0:
        # Probablemente exportó vectores de 3 componentes (u, v, w)
        num_components = len(values) // (N * M)
        field_3d = np.array(values).reshape(M, N, num_components)
        field = field_3d[:, :, 0] # Toma la magnitud principal
    else:
        # Formato normal escalar
        field = np.array(values[:N * M]).reshape(M, N)   # (j, i)
        
    # Cell centres
    xc = 0.5 * (x_coords[:-1] + x_coords[1:])
    yc = 0.5 * (y_coords[:-1] + y_coords[1:])
    return xc, yc, field


def list_reynolds_cases(output_dir):
    cases = []
    if not os.path.isdir(output_dir):
        return cases
    for name in os.listdir(output_dir):
        m = re.fullmatch(r"Re_(\d+)", name)
        if m and os.path.isdir(os.path.join(output_dir, name)):
            cases.append(int(m.group(1)))
    return sorted(cases)


def load_case(case_dir):
    df_u = pd.read_csv(os.path.join(case_dir, "centerline_u.csv"))
    df_v = pd.read_csv(os.path.join(case_dir, "centerline_v.csv"))
    xc_p, yc_p, p_field = read_vtk_rectilinear(os.path.join(case_dir, "pressure.vtk"))
    xc_u, yc_u, u_field = read_vtk_rectilinear(os.path.join(case_dir, "u_velocity.vtk"))
    xc_v, yc_v, v_field = read_vtk_rectilinear(os.path.join(case_dir, "v_velocity.vtk"))
    return df_u, df_v, xc_p, yc_p, p_field, xc_u, yc_u, u_field, xc_v, yc_v, v_field


def plot_case(re_value, case_dir):
    df_u, df_v, xc_p, yc_p, p_field, xc_u, yc_u, u_field, xc_v, yc_v, v_field = load_case(case_dir)

    ghia_col = ghia_column_for_re(re_value)
    ghia_u = np.column_stack([ghia_u_table[:, 0], ghia_u_table[:, 1 + ghia_col]])
    ghia_v = np.column_stack([ghia_v_table[:, 0], ghia_v_table[:, 1 + ghia_col]])

    fig = plt.figure(figsize=(16, 10))
    fig.suptitle(f"Lid-Driven Cavity  —  Re = {re_value}", fontsize=14, fontweight="bold")
    gs = gridspec.GridSpec(2, 3, figure=fig, hspace=0.38, wspace=0.32)

    ax_p = fig.add_subplot(gs[0, 0])
    Xp, Yp = np.meshgrid(xc_p, yc_p)
    cp = ax_p.contourf(Xp, Yp, p_field, levels=40, cmap="RdBu_r")
    fig.colorbar(cp, ax=ax_p, label="p")
    ax_p.set_title("Presion")
    ax_p.set_xlabel("x"); ax_p.set_ylabel("y"); ax_p.set_aspect("equal")

    ax_u2d = fig.add_subplot(gs[0, 1])
    Xu, Yu = np.meshgrid(xc_u, yc_u)
    cu = ax_u2d.contourf(Xu, Yu, u_field, levels=40, cmap="coolwarm")
    fig.colorbar(cu, ax=ax_u2d, label="u")
    ax_u2d.set_title("Velocidad u")
    ax_u2d.set_xlabel("x"); ax_u2d.set_ylabel("y"); ax_u2d.set_aspect("equal")

    ax_v2d = fig.add_subplot(gs[0, 2])
    Xv, Yv = np.meshgrid(xc_v, yc_v)
    cv = ax_v2d.contourf(Xv, Yv, v_field, levels=40, cmap="coolwarm")
    fig.colorbar(cv, ax=ax_v2d, label="v")
    ax_v2d.set_title("Velocidad v")
    ax_v2d.set_xlabel("x"); ax_v2d.set_ylabel("y"); ax_v2d.set_aspect("equal")

    ax_uc = fig.add_subplot(gs[1, 0])
    ax_uc.plot(df_u["u_norm"], df_u["y_norm"], color="steelblue", linewidth=2, label="Solver")
    ax_uc.scatter(ghia_u[:, 1], ghia_u[:, 0], color="red", zorder=5, s=40, label=f"Ghia Re={GHIA_RE_LIST[ghia_col]}")
    ax_uc.set_xlabel("u / U_lid"); ax_uc.set_ylabel("y / H")
    ax_uc.set_title("u(y)  x=L/2"); ax_uc.grid(True, alpha=0.3); ax_uc.legend(fontsize=8)

    ax_vc = fig.add_subplot(gs[1, 1])
    ax_vc.plot(df_v["x_norm"], df_v["v_norm"], color="steelblue", linewidth=2, label="Solver")
    ax_vc.scatter(ghia_v[:, 0], ghia_v[:, 1], color="red", zorder=5, s=40, label=f"Ghia Re={GHIA_RE_LIST[ghia_col]}")
    ax_vc.set_xlabel("x / L"); ax_vc.set_ylabel("v / U_lid")
    ax_vc.set_title("v(x)  y=H/2"); ax_vc.grid(True, alpha=0.3); ax_vc.legend(fontsize=8)

    ax_str = fig.add_subplot(gs[1, 2])
    ny, nx = len(yc_p), len(xc_p)
    u_on_p = 0.5 * (u_field[:, :-1] + u_field[:, 1:]) if u_field.shape[1] > nx else u_field[:ny, :nx]
    v_on_p = 0.5 * (v_field[:-1, :] + v_field[1:, :]) if v_field.shape[0] > ny else v_field[:ny, :nx]
    u_on_p = u_on_p[:ny, :nx]
    v_on_p = v_on_p[:ny, :nx]
    x_uni = np.linspace(xc_p[0], xc_p[-1], nx)
    y_uni = np.linspace(yc_p[0], yc_p[-1], ny)
    interp_u = RegularGridInterpolator((yc_p, xc_p), u_on_p, method="linear")
    interp_v = RegularGridInterpolator((yc_p, xc_p), v_on_p, method="linear")
    Xuu, Yuu = np.meshgrid(x_uni, y_uni)
    pts = np.stack([Yuu.ravel(), Xuu.ravel()], axis=1)
    u_uni = interp_u(pts).reshape(ny, nx)
    v_uni = interp_v(pts).reshape(ny, nx)
    speed = np.sqrt(u_uni ** 2 + v_uni ** 2)
    ax_str.streamplot(x_uni, y_uni, u_uni, v_uni, color=speed, cmap="inferno", density=1.5, linewidth=0.8)
    ax_str.set_title("Lineas de corriente")
    ax_str.set_xlabel("x"); ax_str.set_ylabel("y"); ax_str.set_aspect("equal")

    save_path = os.path.join(case_dir, "validation.png")
    plt.savefig(save_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"Figura guardada: {save_path}")
    return df_u


def plot_multi_re_u(output_dir, re_cases, u_profiles):
    plt.figure(figsize=(7, 6))
    for re_value, df_u in zip(re_cases, u_profiles):
        plt.plot(df_u["u_norm"], df_u["y_norm"], linewidth=1.8, label=f"Solver Re={re_value}")
        ghia_col = ghia_column_for_re(re_value)
        ghia_u = np.column_stack([ghia_u_table[:, 0], ghia_u_table[:, 1 + ghia_col]])
        plt.scatter(ghia_u[:, 1], ghia_u[:, 0], s=14, alpha=0.8)
    plt.xlabel("u / U_lid")
    plt.ylabel("y / H")
    plt.title("Comparacion de perfiles u(x=L/2,y)")
    plt.grid(True, alpha=0.3)
    plt.legend(fontsize=8)
    save_path = os.path.join(output_dir, "validation_all_re_u.png")
    plt.tight_layout()
    plt.savefig(save_path, dpi=150)
    plt.close()
    print(f"Figura global guardada: {save_path}")


def main():
    print(f"Buscando resultados en: {OUTPUT_DIR}")
    re_cases = list_reynolds_cases(OUTPUT_DIR)
    if not re_cases:
        raise RuntimeError("No se encontraron carpetas output/Re_<N>.")

    print(f"Casos detectados: {re_cases}")
    u_profiles = []
    for re_value in re_cases:
        case_dir = os.path.join(OUTPUT_DIR, f"Re_{re_value}")
        u_profiles.append(plot_case(re_value, case_dir))

    plot_multi_re_u(OUTPUT_DIR, re_cases, u_profiles)


if __name__ == "__main__":
    main()