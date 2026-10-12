"""
plot_square_cylinder.py — visualise Square Cylinder results.

Scans output/square_cylinder/Re_<value>/ folders produced by the solver
(input/scenario.txt = square_cylinder) and, for each case, draws:
  - vorticity field (finite-difference curl) with the obstacle outline overlaid
  - Cd(t) / Cl(t) time series (from forces.csv)
and produces one summary plot across the Re sweep:
  - Strouhal number vs Re, with the literature reference band (~0.13-0.15 at
    Re ~ 100-200) annotated for visual reference only — NOT a claim of
    agreement. See docs/scenarios/square_cylinder.md for validation status.

Reuses read_vtk_rectilinear() from plot_results.py (same ASCII VTK format).
Run from anywhere:  python3 plotters/plot_square_cylinder.py
"""
import os
import re
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.patches as patches
from scipy.interpolate import RegularGridInterpolator

from plot_results import read_vtk_rectilinear  # reuse the VTK reader

BASE_DIR = os.path.dirname(os.path.abspath(__file__))
PROJECT_ROOT = os.path.dirname(BASE_DIR)
OUTPUT_DIR = os.path.join(PROJECT_ROOT, "output", "square_cylinder")

# Literature reference band for a square cylinder, Re ~ 100-200 (for visual
# reference only — this solver's results are NOT yet validated against it).
ST_LITERATURE_BAND = (0.13, 0.15)
ST_LITERATURE_RE_RANGE = (100, 200)


def list_reynolds_cases(output_dir):
    cases = []
    if not os.path.isdir(output_dir):
        return cases
    for name in os.listdir(output_dir):
        m = re.fullmatch(r"Re_(\d+)", name)
        if m and os.path.isdir(os.path.join(output_dir, name)):
            cases.append(int(m.group(1)))
    return sorted(cases)


def read_obstacle_mask(output_root):
    mask_path = os.path.join(output_root, "obstacle_mask.vtk")
    if not os.path.exists(mask_path):
        return None
    return read_vtk_rectilinear(mask_path)


def obstacle_outline_from_mask(xc, yc, mask):
    """Returns (x0, x1, y0, y1) bounding box of the solid region, or None."""
    solid = mask > 0.5
    if not solid.any():
        return None
    rows = np.any(solid, axis=1)
    cols = np.any(solid, axis=0)
    y_idx = np.where(rows)[0]
    x_idx = np.where(cols)[0]
    return (xc[x_idx[0]], xc[x_idx[-1]], yc[y_idx[0]], yc[y_idx[-1]])


def interp_to_grid(xc, yc, field, x_uni, y_uni):
    interp = RegularGridInterpolator((yc, xc), field,
                                     bounds_error=False, fill_value=None)
    X, Y = np.meshgrid(x_uni, y_uni)
    pts = np.stack([Y.ravel(), X.ravel()], axis=1)
    return interp(pts).reshape(len(y_uni), len(x_uni))


def plot_case(re_value, case_dir, obstacle_box):
    xc_u, yc_u, U = read_vtk_rectilinear(os.path.join(case_dir, "u_velocity.vtk"))
    xc_v, yc_v, V = read_vtk_rectilinear(os.path.join(case_dir, "v_velocity.vtk"))

    fig = plt.figure(figsize=(13, 8))
    fig.suptitle(f"Square Cylinder  —  Re = {re_value}  (validation pending)",
                 fontsize=14, fontweight="bold")
    gs = gridspec.GridSpec(2, 1, figure=fig, height_ratios=[1.3, 1], hspace=0.35)

    # --- Vorticity field (finite-difference curl on a common uniform grid) ---
    ax_w = fig.add_subplot(gs[0, 0])
    nx, ny = max(len(xc_u), len(xc_v)), max(len(yc_u), len(yc_v))
    x_uni = np.linspace(min(xc_u[0], xc_v[0]), max(xc_u[-1], xc_v[-1]), nx)
    y_uni = np.linspace(min(yc_u[0], yc_v[0]), max(yc_u[-1], yc_v[-1]), ny)
    u_uni = interp_to_grid(xc_u, yc_u, U, x_uni, y_uni)
    v_uni = interp_to_grid(xc_v, yc_v, V, x_uni, y_uni)
    dvdx = np.gradient(v_uni, x_uni, axis=1)
    dudy = np.gradient(u_uni, y_uni, axis=0)
    vorticity = dvdx - dudy

    X, Y = np.meshgrid(x_uni, y_uni)
    vmax = np.nanpercentile(np.abs(vorticity), 95) or 1.0
    cf = ax_w.contourf(X, Y, vorticity, levels=40, cmap="RdBu_r", vmin=-vmax, vmax=vmax)
    fig.colorbar(cf, ax=ax_w, label=r"$\omega = \partial v/\partial x - \partial u/\partial y$")
    if obstacle_box is not None:
        x0, x1, y0, y1 = obstacle_box
        ax_w.add_patch(patches.Rectangle((x0, y0), x1 - x0, y1 - y0,
                                          facecolor="black", edgecolor="white", linewidth=1.0))
    ax_w.set_title("Vorticity field")
    ax_w.set_xlabel("x"); ax_w.set_ylabel("y"); ax_w.set_aspect("equal")

    # --- Cd(t) / Cl(t) time series ---
    # The first few steps carry a large impulsive-start transient (the field
    # begins at rest) that dwarfs the later physical range on a linear axis —
    # skip a short warm-up window so the plot is legible.
    ax_f = fig.add_subplot(gs[1, 0])
    forces_path = os.path.join(case_dir, "forces.csv")
    if os.path.exists(forces_path):
        df = pd.read_csv(forces_path)
        skip = min(10, max(0, len(df) // 20))
        df_plot = df.iloc[skip:]
        ax_f.plot(df_plot["time"], df_plot["Cd"], label="$C_d$", color="tab:blue", linewidth=0.8)
        ax_f.plot(df_plot["time"], df_plot["Cl"], label="$C_l$", color="tab:orange", linewidth=0.8)
        ax_f.set_xlabel("time"); ax_f.set_ylabel("coefficient")
        ax_f.legend(loc="upper right", fontsize=8)
        ax_f.grid(True, alpha=0.3)
        if skip > 0:
            ax_f.set_title(f"Drag / lift coefficients (first {skip} steps' impulsive-start transient omitted)")
        else:
            ax_f.set_title("Drag / lift coefficients")
    else:
        ax_f.set_title("Drag / lift coefficients")

    strouhal_path = os.path.join(case_dir, "strouhal.csv")
    if os.path.exists(strouhal_path):
        st = pd.read_csv(strouhal_path).iloc[0]
        ax_f.text(0.02, 0.98, f"St (this run, unvalidated) = {st['St']:.4f}",
                  transform=ax_f.transAxes, va="top", ha="left", fontsize=8,
                  family="monospace", bbox=dict(boxstyle="round", fc="white", alpha=0.8))

    save_path = os.path.join(case_dir, "square_cylinder_plot.png")
    plt.savefig(save_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"Figure saved: {save_path}")


def plot_strouhal_summary(cases, output_dir):
    re_vals, st_vals = [], []
    for re in cases:
        path = os.path.join(output_dir, f"Re_{re}", "strouhal.csv")
        if not os.path.exists(path):
            continue
        st = pd.read_csv(path).iloc[0]["St"]
        re_vals.append(re)
        st_vals.append(st)
    if not re_vals:
        return

    fig, ax = plt.subplots(figsize=(6, 5))
    ax.scatter(re_vals, st_vals, color="tab:blue", label="This solver (unvalidated)", zorder=3)
    ax.axhspan(*ST_LITERATURE_BAND, xmin=0, xmax=1, color="tab:green", alpha=0.15,
               label=f"Literature band {ST_LITERATURE_BAND} (Re {ST_LITERATURE_RE_RANGE})")
    ax.set_xlabel("Re"); ax.set_ylabel("Strouhal number St")
    ax.set_title("Strouhal number vs Re — validation pending, NOT yet calibrated")
    ax.legend(fontsize=8)
    ax.grid(True, alpha=0.3)
    save_path = os.path.join(output_dir, "strouhal_summary.png")
    plt.savefig(save_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"Figure saved: {save_path}")


def main():
    print(f"Looking for results in: {OUTPUT_DIR}")
    cases = list_reynolds_cases(OUTPUT_DIR)
    if not cases:
        raise RuntimeError("No output/square_cylinder/Re_<value> folders found. "
                           "Run the solver with input/scenario.txt = square_cylinder first.")
    print(f"Detected cases: {cases}")

    mask_result = read_obstacle_mask(OUTPUT_DIR)
    obstacle_box = None
    if mask_result is not None:
        xc_m, yc_m, mask = mask_result
        obstacle_box = obstacle_outline_from_mask(xc_m, yc_m, mask)

    for re in cases:
        plot_case(re, os.path.join(OUTPUT_DIR, f"Re_{re}"), obstacle_box)

    plot_strouhal_summary(cases, OUTPUT_DIR)


if __name__ == "__main__":
    main()
