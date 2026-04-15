"""
plot_convergence.py — Convergence-rate plots for discrete operator verification
Reads output/convergence.csv (written by ./build/operator_verify) and produces
output/convergence_rates.png, a 2×3 figure analogous to Fig. 5 of the reference
paper (Trias et al. style: log-log error vs. mesh step h = 1/N).

Usage:
    ./build/operator_verify          # writes output/convergence.csv
    python3 plotters/plot_convergence.py

Requires: numpy, matplotlib, pandas
"""

import os
import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

# ---------------------------------------------------------------------------
# Paths (project root is one level above plotters/)
# ---------------------------------------------------------------------------
PLOTTERS_DIR = os.path.dirname(os.path.abspath(__file__))
BASE_DIR     = os.path.dirname(PLOTTERS_DIR)
OUTPUT_DIR   = os.path.join(BASE_DIR, "output")
CSV_PATH   = os.path.join(OUTPUT_DIR, "convergence.csv")

if not os.path.exists(CSV_PATH):
    print(f"Error: {CSV_PATH} not found.")
    print("Run  ./build/operator_verify  first to generate the data.")
    sys.exit(1)

df = pd.read_csv(CSV_PATH)
df["h"] = 1.0 / df["N"]    # mesh step h = 1/N

# Convenience: get error array for a named suite
def suite(name):
    s = df[df["suite"] == name].sort_values("N")
    return s["h"].values, s["error_linf"].values, s["N"].values

# ---------------------------------------------------------------------------
# Reference lines O(h^p) anchored at first data point
# ---------------------------------------------------------------------------
def ref_line(h, e0, h0, p, label=None, ax=None, **kw):
    e_ref = e0 * (h / h0) ** p
    ax.plot(h, e_ref, **kw)

# Slope annotation in log-log space
def annotate_slope(ax, h, e, color="black"):
    rates = []
    for k in range(1, len(e)):
        if e[k] > 0 and e[k-1] > 0:
            r = np.log2(e[k-1] / e[k])
            rates.append(r)
    if rates:
        mean_r = np.mean(rates)
        # Place annotation near middle of the line
        mid = len(h) // 2
        ax.annotate(f"slope ≈ {mean_r:.2f}", xy=(h[mid], e[mid]),
                    xytext=(h[mid]*1.3, e[mid]*2.0),
                    fontsize=8, color=color,
                    arrowprops=dict(arrowstyle="->", color=color, lw=0.8))

# ---------------------------------------------------------------------------
# Figure
# ---------------------------------------------------------------------------
fig = plt.figure(figsize=(16, 10))
fig.suptitle("Operator Verification — Convergence Rates\n"
             "(L∞ error vs mesh step h = 1/N,   reference lines O(h¹) and O(h²))",
             fontsize=13, fontweight="bold")
gs = gridspec.GridSpec(2, 3, figure=fig, hspace=0.42, wspace=0.35)

COLORS = {
    "Laplace_uniform": "#1f77b4",
    "Gradient_dpx":    "#2ca02c",
    "Gradient_dpy":    "#d62728",
    "Divergence":      "#9467bd",
    "Diffusion_U":     "#8c564b",
    "Diffusion_V":     "#e377c2",
    "Laplace_stretched": "#ff7f0e",
}

# ---------------------------------------------------------------------------
# Panel (0,0) — LaplaceOperator uniform
# ---------------------------------------------------------------------------
ax = fig.add_subplot(gs[0, 0])
h, e, N = suite("Laplace_uniform")
ax.loglog(h, e, "o-", color=COLORS["Laplace_uniform"], lw=2, ms=7, label="Laplace (uniform)")
ref_line(h, e[0], h[0], 2, ax=ax, color="gray", ls="--", lw=1, label="O(h²)")
ref_line(h, e[0], h[0], 1, ax=ax, color="lightgray", ls=":", lw=1, label="O(h¹)")
annotate_slope(ax, h, e, color=COLORS["Laplace_uniform"])
ax.set_xlabel("h = 1/N"); ax.set_ylabel("L∞ error")
ax.set_title("Suite 1 — LaplaceOperator\nφ = cos(πx)cos(πy),  Neumann BC")
ax.legend(fontsize=8); ax.grid(True, alpha=0.3, which="both")

# ---------------------------------------------------------------------------
# Panel (0,1) — GradientOperator
# ---------------------------------------------------------------------------
ax = fig.add_subplot(gs[0, 1])
h_x, e_x, _ = suite("Gradient_dpx")
h_y, e_y, _ = suite("Gradient_dpy")
ax.loglog(h_x, e_x, "s-", color=COLORS["Gradient_dpx"], lw=2, ms=7, label="∂p/∂x")
ax.loglog(h_y, e_y, "^-", color=COLORS["Gradient_dpy"], lw=2, ms=7, label="∂p/∂y")
ref_line(h_x, e_x[0], h_x[0], 2, ax=ax, color="gray", ls="--", lw=1, label="O(h²)")
ref_line(h_x, e_x[0], h_x[0], 1, ax=ax, color="lightgray", ls=":", lw=1, label="O(h¹)")
annotate_slope(ax, h_x, e_x, color=COLORS["Gradient_dpx"])
ax.set_xlabel("h = 1/N"); ax.set_ylabel("L∞ error")
ax.set_title("Suite 2 — GradientOperator\nφ = cos(πx)cos(πy),  interior nodes")
ax.legend(fontsize=8); ax.grid(True, alpha=0.3, which="both")

# ---------------------------------------------------------------------------
# Panel (0,2) — DivergenceOperator
# ---------------------------------------------------------------------------
ax = fig.add_subplot(gs[0, 2])
h, e, _ = suite("Divergence")
ax.loglog(h, e, "D-", color=COLORS["Divergence"], lw=2, ms=7, label="Divergence")
ref_line(h, e[0], h[0], 2, ax=ax, color="gray", ls="--", lw=1, label="O(h²)")
ref_line(h, e[0], h[0], 1, ax=ax, color="lightgray", ls=":", lw=1, label="O(h¹)")
annotate_slope(ax, h, e, color=COLORS["Divergence"])
ax.set_xlabel("h = 1/N"); ax.set_ylabel("L∞ error")
ax.set_title("Suite 3 — DivergenceOperator\nu = sin(2πx)cos(πy),  face sampling")
ax.legend(fontsize=8); ax.grid(True, alpha=0.3, which="both")

# ---------------------------------------------------------------------------
# Panel (1,0) — DiffusionOperator U and V
# ---------------------------------------------------------------------------
ax = fig.add_subplot(gs[1, 0])
h_u, e_u, _ = suite("Diffusion_U")
h_v, e_v, _ = suite("Diffusion_V")
ax.loglog(h_u, e_u, "o-", color=COLORS["Diffusion_U"], lw=2, ms=7, label="applyU  u=sin(πy)")
ax.loglog(h_v, e_v, "s--", color=COLORS["Diffusion_V"], lw=2, ms=7, label="applyV  v=sin(πx)")
ref_line(h_u, e_u[0], h_u[0], 2, ax=ax, color="gray", ls="--", lw=1, label="O(h²)")
ref_line(h_u, e_u[0], h_u[0], 1, ax=ax, color="lightgray", ls=":", lw=1, label="O(h¹)")
annotate_slope(ax, h_u, e_u, color=COLORS["Diffusion_U"])
ax.set_xlabel("h = 1/N"); ax.set_ylabel("L∞ error")
ax.set_title("Suites 4 & 6 — DiffusionOperator\nν∇²u (north/south) and ν∇²v (east/west)")
ax.legend(fontsize=8); ax.grid(True, alpha=0.3, which="both")

# ---------------------------------------------------------------------------
# Panel (1,1) — Algebraic identities (machine-precision residuals)
# ---------------------------------------------------------------------------
ax = fig.add_subplot(gs[1, 1])
# Suite 5: D·G identity
h5, e5, N5 = suite("DGL_identity")
ax.semilogy(N5, e5, "o-", color="#1f77b4", lw=2, ms=7, label="Suite 5: |D·G·φ + L·φ/V|∞")
# Suite 7: linear diffusion (should be exactly 0; show as points at 1e-16 floor)
h7, e7, N7 = suite("Diffusion_linear_BC")
e7_plot = np.maximum(e7, 1e-17)   # avoid log(0)
ax.semilogy(N7, e7_plot, "s-", color="#2ca02c", lw=2, ms=7, label="Suite 7: |ν∇²(y)|∞ = 0 (exact)")
# Suite 8: adjoint
h8, e8, N8 = suite("Adjoint_GD")
e8_plot = np.maximum(e8, 1e-17)
ax.semilogy(N8, e8_plot, "^-", color="#d62728", lw=2, ms=7, label="Suite 8: |φᵀDu + (Gφ)ᵀu|")
# Suite 9: SPD symmetry
h9, e9, N9 = suite("SPD_symmetry")
ax.semilogy(N9, e9, "D-", color="#9467bd", lw=2, ms=7, label="Suite 9: |φᵀAψ − ψᵀAφ|")
# Machine epsilon reference
ax.axhline(2.2e-16, color="black", ls=":", lw=1, label="machine ε ≈ 2.2×10⁻¹⁶")
ax.set_xlabel("N")
ax.set_ylabel("Residual (log scale)")
ax.set_title("Suites 5, 7–9 — Algebraic Identities\n(machine-precision residuals)")
ax.legend(fontsize=7.5); ax.grid(True, alpha=0.3, which="both")

# ---------------------------------------------------------------------------
# Panel (1,2) — Uniform vs. stretched LaplaceOperator
# ---------------------------------------------------------------------------
ax = fig.add_subplot(gs[1, 2])
h1, e1, _ = suite("Laplace_uniform")
h10, e10, _ = suite("Laplace_stretched")
ax.loglog(h1,  e1,  "o-",  color=COLORS["Laplace_uniform"],   lw=2, ms=7, label="Suite 1: uniform mesh")
ax.loglog(h10, e10, "s--", color=COLORS["Laplace_stretched"],  lw=2, ms=7, label="Suite 10: tanh stretch (γ=2)")
ref_line(h1, e1[0], h1[0], 2, ax=ax, color="gray", ls="--", lw=1, label="O(h²)")
ref_line(h10, e10[0], h10[0], 1, ax=ax, color="orange", ls=":", lw=1.5, label="O(h¹)")
annotate_slope(ax, h1,  e1,  color=COLORS["Laplace_uniform"])
annotate_slope(ax, h10, e10, color=COLORS["Laplace_stretched"])
ax.set_xlabel("h = 1/N"); ax.set_ylabel("L∞ error")
ax.set_title("Suites 1 & 10 — Laplacian: uniform vs. stretched\n"
             "(stretched → O(h¹); fix: use 2/(h_e+h_w) factor)")
ax.legend(fontsize=8); ax.grid(True, alpha=0.3, which="both")

# ---------------------------------------------------------------------------
# Save
# ---------------------------------------------------------------------------
save_path = os.path.join(OUTPUT_DIR, "convergence_rates.png")
plt.savefig(save_path, dpi=150, bbox_inches="tight")
print(f"Figure saved to: {save_path}")
plt.show()
