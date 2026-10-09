"""
plot_dhc.py — visualise Differentially Heated Cavity results.

Scans output/Ra_<value>/ folders produced by the solver (solve_energy = 1) and,
for each case, draws:
  - temperature field with isotherms
  - streamlines coloured by speed
and annotates the de Vahl Davis metrics from dhc_metrics.csv.

Reuses read_vtk_rectilinear() from plot_results.py (same ASCII VTK format).
Run from anywhere:  python3 plotters/plot_dhc.py
"""
import os
import re
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from scipy.interpolate import RegularGridInterpolator

from plot_results import read_vtk_rectilinear  # reuse the VTK reader

BASE_DIR = os.path.dirname(os.path.abspath(__file__))
PROJECT_ROOT = os.path.dirname(BASE_DIR)
OUTPUT_DIR = os.path.join(PROJECT_ROOT, "output")

# de Vahl Davis (1983) benchmark — average Nusselt on the hot wall
DVD_NU = {1000: 1.118, 10000: 2.243, 100000: 4.519, 1000000: 8.800}


def list_rayleigh_cases(output_dir):
    cases = []
    if not os.path.isdir(output_dir):
        return cases
    for name in os.listdir(output_dir):
        m = re.fullmatch(r"Ra_(\d+)", name)
        if m and os.path.isdir(os.path.join(output_dir, name)):
            cases.append(int(m.group(1)))
    return sorted(cases)


def interp_to_grid(xc, yc, field, x_uni, y_uni):
    """Interpolate a (len(yc), len(xc)) field onto a uniform grid for streamplot."""
    interp = RegularGridInterpolator((yc, xc), field,
                                     bounds_error=False, fill_value=None)
    X, Y = np.meshgrid(x_uni, y_uni)
    pts = np.stack([Y.ravel(), X.ravel()], axis=1)
    return interp(pts).reshape(len(y_uni), len(x_uni))


def plot_case(ra_value, case_dir):
    xc_t, yc_t, T = read_vtk_rectilinear(os.path.join(case_dir, "temperature.vtk"))
    xc_u, yc_u, U = read_vtk_rectilinear(os.path.join(case_dir, "u_velocity.vtk"))
    xc_v, yc_v, V = read_vtk_rectilinear(os.path.join(case_dir, "v_velocity.vtk"))

    fig = plt.figure(figsize=(13, 5.5))
    fig.suptitle(f"Differentially Heated Cavity  —  Ra = {ra_value:g}",
                 fontsize=14, fontweight="bold")
    gs = gridspec.GridSpec(1, 2, figure=fig, wspace=0.28)

    # --- Temperature field + isotherms ---
    ax_T = fig.add_subplot(gs[0, 0])
    XT, YT = np.meshgrid(xc_t, yc_t)
    cf = ax_T.contourf(XT, YT, T, levels=30, cmap="coolwarm")
    ax_T.contour(XT, YT, T, levels=15, colors="k", linewidths=0.4)
    fig.colorbar(cf, ax=ax_T, label=r"$\theta$")
    ax_T.set_title("Temperature (hot left, cold right)")
    ax_T.set_xlabel("x"); ax_T.set_ylabel("y"); ax_T.set_aspect("equal")

    # --- Streamlines on a uniform grid (u, v interpolated onto P-grid extent) ---
    ax_S = fig.add_subplot(gs[0, 1])
    nx, ny = len(xc_t), len(yc_t)
    x_uni = np.linspace(xc_t[0], xc_t[-1], nx)
    y_uni = np.linspace(yc_t[0], yc_t[-1], ny)
    u_uni = interp_to_grid(xc_u, yc_u, U, x_uni, y_uni)
    v_uni = interp_to_grid(xc_v, yc_v, V, x_uni, y_uni)
    speed = np.sqrt(u_uni ** 2 + v_uni ** 2)
    strm = ax_S.streamplot(x_uni, y_uni, u_uni, v_uni, color=speed,
                           cmap="inferno", density=1.6, linewidth=0.8)
    fig.colorbar(strm.lines, ax=ax_S, label="|u|")
    ax_S.set_title("Streamlines")
    ax_S.set_xlabel("x"); ax_S.set_ylabel("y"); ax_S.set_aspect("equal")

    # --- Annotate metrics ---
    metrics_path = os.path.join(case_dir, "dhc_metrics.csv")
    if os.path.exists(metrics_path):
        m = pd.read_csv(metrics_path).iloc[0]
        ref = DVD_NU.get(int(round(ra_value)))
        txt = (f"u_max = {m['u_max']:.3f}\n"
               f"v_max = {m['v_max']:.3f}\n"
               f"Nu_avg = {m['Nu_avg']:.4f}")
        if ref is not None:
            txt += f"\nNu (deVahlDavis) = {ref:.3f}"
        ax_S.text(0.02, 0.98, txt, transform=ax_S.transAxes, va="top", ha="left",
                  fontsize=8, family="monospace",
                  bbox=dict(boxstyle="round", fc="white", alpha=0.8))

    save_path = os.path.join(case_dir, "dhc_plot.png")
    plt.savefig(save_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"Figure saved: {save_path}")


def main():
    print(f"Looking for results in: {OUTPUT_DIR}")
    cases = list_rayleigh_cases(OUTPUT_DIR)
    if not cases:
        raise RuntimeError("No output/Ra_<value> folders found. "
                           "Run the solver with solve_energy = 1 first.")
    print(f"Detected cases: {cases}")
    for ra in cases:
        plot_case(ra, os.path.join(OUTPUT_DIR, f"Ra_{ra}"))


if __name__ == "__main__":
    main()
