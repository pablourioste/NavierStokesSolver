"""
plot_convergence.py — convergence of the discrete operators (Method of
Manufactured Solutions).

Reads  output/verification/data/convergence.csv   (written by ./build/operator_verify)
Writes output/verification/figures/convergence_rates.{pdf,png}
         (a) truncation error of each operator on a uniform mesh
         (b) Laplacian on a uniform vs. a stretched (tanh, gamma = 2) mesh
         (c) residuals of the discrete algebraic identities

The observed order p is the least-squares slope of log(error) vs log(h) over
the two finest meshes, and is reported in the legend.

Usage:
    bash run.sh verify
    python3 plotters/plot_convergence.py
"""
import os
import sys

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

import nsplot as ns

VERIFY_DIR = os.path.join(ns.OUTPUT_ROOT, "verification")
CSV_PATH = os.path.join(VERIFY_DIR, "data", "convergence.csv")
MACHINE_EPS = np.finfo(float).eps


def suite(df, name):
    s = df[df["suite"] == name].sort_values("N")
    return 1.0 / s["N"].values, s["error_linf"].values, s["N"].values


def observed_order(h, e):
    return np.polyfit(np.log(h[-2:]), np.log(e[-2:]), 1)[0]


def slope_triangle(ax, h0, h1, e0, p, upper=False):
    """Reference triangle of slope p spanning h0..h1, with its horizontal leg
    at e0 (below the data) or, if `upper`, its horizontal leg on top."""
    e1 = e0 * (h1 / h0) ** p
    if upper:
        xs, ys = [h0, h0, h1, h0], [e0, e1, e1, e0]
        ax.text(h0 / 1.07, np.sqrt(e0 * e1), f"${p}$", ha="right", va="center", fontsize=8)
        ax.text(np.sqrt(h0 * h1), e1 * 1.12, "$1$", ha="center", va="bottom", fontsize=8)
    else:
        xs, ys = [h0, h1, h1, h0], [e0, e0, e1, e0]
        ax.text(h1 * 1.07, np.sqrt(e0 * e1), f"${p}$", ha="left", va="center", fontsize=8)
        ax.text(np.sqrt(h0 * h1), e0 / 1.12, "$1$", ha="center", va="top", fontsize=8)
    ax.plot(xs, ys, color="k", lw=0.5)


def legend_below(ax):
    ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.27), fontsize=7.5,
              handlelength=2.2, borderaxespad=0.0)


def series(ax, h, e, marker, label, ls="-"):
    p = observed_order(h, e)
    ax.loglog(h, e, ls=ls, color="k", lw=0.8, marker=marker, mfc="white", mec="k",
              ms=4, mew=0.7, label=rf"{label} ($p = {p:.2f}$)")


def h_axis(ax):
    ax.set_xlabel(r"$h = 1/N$")
    ax.set_xticks([1 / 64, 1 / 32, 1 / 16, 1 / 8],
                  [r"$2^{-6}$", r"$2^{-5}$", r"$2^{-4}$", r"$2^{-3}$"])
    ax.xaxis.set_minor_locator(plt.NullLocator())
    ax.set_xlim(1 / 64 / 1.4, 1 / 8 * 1.4)


def main():
    if not os.path.exists(CSV_PATH):
        sys.exit(f"{CSV_PATH} not found — run  bash run.sh verify  first.")
    ns.apply_style()
    df = pd.read_csv(CSV_PATH)

    fig, (ax_a, ax_b, ax_c) = plt.subplots(1, 3, figsize=(ns.TEXTWIDTH, 2.3))
    fig.subplots_adjust(wspace=0.42)

    # (a) every operator on a uniform mesh. dp/dy and the v-diffusion are
    # bitwise identical to dp/dx and the u-diffusion (symmetric test fields),
    # so one series each is drawn.
    for name, marker, label in (("Laplace_uniform", "o", r"Laplacian $L$"),
                                ("Gradient_dpx", "s", r"Gradient $G$"),
                                ("Divergence", "^", r"Divergence $D$"),
                                ("Diffusion_U", "D", r"Diffusion $\nu\nabla^2$")):
        h, e, _ = suite(df, name)
        series(ax_a, h, e, marker, label)
    slope_triangle(ax_a, 0.05, 0.11, 4e-4, 2)
    ax_a.set_ylim(1e-4, 1)
    h_axis(ax_a)
    ax_a.set_ylabel(r"$\|\varepsilon\|_\infty$")
    ns.panel_title(ax_a, "a", "Uniform mesh")
    legend_below(ax_a)

    # (b) uniform vs. stretched Laplacian
    h, e, _ = suite(df, "Laplace_uniform")
    series(ax_b, h, e, "o", "Uniform")
    slope_triangle(ax_b, 0.05, 0.11, 5e-3, 2)
    h, e, _ = suite(df, "Laplace_stretched")
    series(ax_b, h, e, "s", "Stretched", ls="--")
    slope_triangle(ax_b, 0.019, 0.045, 3.0, 1, upper=True)
    h_axis(ax_b)
    ns.panel_title(ax_b, "b", "Mesh stretching")
    legend_below(ax_b)
    ax_b.set_ylim(1e-3, 1e2)

    # (c) algebraic identities (machine precision; exact zeros are not drawn)
    for name, marker, label in (("DGL_identity", "o", r"$\|D G\phi + L\phi/V\|_\infty$"),
                                ("Adjoint_GD", "s", r"$|\phi^{\mathsf T} D u + (G\phi)^{\mathsf T} u|$"),
                                ("SPD_symmetry", "^", r"$|\phi^{\mathsf T} A\psi - \psi^{\mathsf T} A\phi|$")):
        _, e, N = suite(df, name)
        keep = e > 0
        ax_c.semilogy(N[keep], e[keep], color="k", lw=0.8, marker=marker, mfc="white",
                      mec="k", ms=4, mew=0.7, label=label)
    ax_c.axhline(MACHINE_EPS, color="0.5", lw=0.6, ls=":")
    ax_c.text(64 * 1.3, MACHINE_EPS / 1.6, r"$\epsilon_{\mathrm{mach}}$", ha="right",
              va="top", fontsize=8, color="0.3")
    ax_c.set_xscale("log", base=2)
    ax_c.set_xticks([8, 16, 32, 64], ["8", "16", "32", "64"])
    ax_c.xaxis.set_minor_locator(plt.NullLocator())
    ax_c.set_xlim(8 / 1.4, 64 * 1.4)
    ax_c.set_ylim(1e-17, 1e-9)
    ax_c.set_xlabel(r"$N$")
    ax_c.set_ylabel("Residual")
    ns.panel_title(ax_c, "c", "Discrete identities")
    legend_below(ax_c)

    ns.save(fig, os.path.join(VERIFY_DIR, "figures"), "convergence_rates")


if __name__ == "__main__":
    main()
