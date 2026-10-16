"""
plot_dhc.py — Differentially Heated Cavity post-processing and validation.

Reads output/differentially_heated_cavity/Ra_<N>/ and writes

  per case   Ra_<N>/figures/flow_fields.{pdf,png}   (a) temperature + isotherms,
                                                     (b) stream function
  summary    figures/fields_overview                isotherms and streamlines vs Ra
             figures/benchmark                      Nu, u_max, v_max vs de Vahl Davis
             validation_summary.csv                 relative errors vs the benchmark

Reference: G. de Vahl Davis, Int. J. Numer. Methods Fluids 3 (1983) 249-264.
Run from anywhere:  python3 plotters/plot_dhc.py
"""
import os

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.ticker import NullFormatter, ScalarFormatter

import nsplot as ns

SCENARIO = "differentially_heated_cavity"

# de Vahl Davis (1983) benchmark: u_max on x = 1/2, v_max on y = 1/2, mean Nu
DVD = {
    1000:    {"u_max": 3.649,  "v_max": 3.697,  "Nu_avg": 1.118},
    10000:   {"u_max": 16.178, "v_max": 19.617, "Nu_avg": 2.243},
    100000:  {"u_max": 34.73,  "v_max": 68.59,  "Nu_avg": 4.519},
    1000000: {"u_max": 64.63,  "v_max": 219.36, "Nu_avg": 8.800},
}
THETA_LEVELS = np.linspace(0, 1, 11)


def padded_temperature(F):
    """theta on cell centres plus wall values: hot west (1), cold east (0),
    adiabatic north/south (zero normal gradient)."""
    x = np.concatenate([[F.xf[0]], F.xc, [F.xf[-1]]])
    y = np.concatenate([[F.yf[0]], F.yc, [F.yf[-1]]])
    T = np.pad(F.T, 1, mode="edge")
    T[:, 0], T[:, -1] = 1.0, 0.0
    return x, y, T


def draw_temperature(fig, ax, F, with_colorbar=True):
    x, y, T = padded_temperature(F)
    X, Y = np.meshgrid(x, y)
    cf = ax.contourf(X, Y, T, levels=np.linspace(0, 1, 21), cmap=ns.CMAP_THERMAL)
    ax.contour(X, Y, T, levels=THETA_LEVELS[1:-1], colors="k", linewidths=0.45)
    frame(ax, F)
    if with_colorbar:
        ns.colorbar(fig, cf, ax, r"$\theta$", ticks=[0, 0.5, 1])
    return cf


def draw_streamfunction(ax, F):
    psi = ns.stream_function(F)
    ext = psi.min() if abs(psi.min()) > abs(psi.max()) else psi.max()
    levels = np.sort(np.linspace(0, ext, 11)[1:])
    X, Y = np.meshgrid(F.xf, F.yf)
    ax.contour(X, Y, psi, levels=levels, colors="k", linewidths=0.55,
               linestyles="solid")
    frame(ax, F)
    return abs(ext)


def frame(ax, F):
    ax.set_xlim(F.xf[0], F.xf[-1])
    ax.set_ylim(F.yf[0], F.yf[-1])
    ns.unit_ticks(ax)
    ns.plain_axes(ax)


def grid_label(F):
    M, N = F.shape
    return rf"${N}\times{M}$"


def plot_flow_fields(ra, F, paths):
    fig, (ax_t, ax_s) = plt.subplots(1, 2, figsize=(0.7 * ns.TEXTWIDTH, 2.1))
    fig.subplots_adjust(wspace=0.5)
    draw_temperature(fig, ax_t, F)
    ns.panel_title(ax_t, "a", r"Temperature $\theta$")
    psi_max = draw_streamfunction(ax_s, F)
    ns.panel_title(ax_s, "b", rf"$\psi$, $|\psi|_{{\max}} = {psi_max:.3g}$")
    for ax in (ax_t, ax_s):
        ax.set_xlabel(r"$x/L$")
    ax_t.set_ylabel(r"$y/L$")
    ax_s.set_yticklabels([])
    ns.heading(fig, [ax_t, ax_s], rf"Differentially heated cavity, $Ra = {ns.sci_tex(ra)}$, "
                                  rf"$Pr = 0.71$, grid {grid_label(F)}")
    ns.save(fig, paths.figures, "flow_fields")


def plot_fields_overview(cases, out_dir):
    n = len(cases)
    w = ns.TEXTWIDTH
    fig, axes = plt.subplots(2, n, figsize=(w * min(1, n / 3), 2 * w / 3 * 1.05),
                             squeeze=False)
    fig.subplots_adjust(wspace=0.12, hspace=0.12)
    cf = None
    for k, c in enumerate(cases):
        cf = draw_temperature(fig, axes[0, k], c["F"], with_colorbar=False)
        draw_streamfunction(axes[1, k], c["F"])
        axes[0, k].set_title(rf"$Ra = {ns.sci_tex(c['ra'])}$")
        axes[0, k].set_xticklabels([])
        axes[1, k].set_xlabel(r"$x/L$")
        if k:
            axes[0, k].set_yticklabels([])
            axes[1, k].set_yticklabels([])
    axes[0, 0].set_ylabel(r"$y/L$")
    axes[1, 0].set_ylabel(r"$y/L$")
    axes[0, 0].text(-0.62, 0.5, "Isotherms", rotation=90, va="center", ha="center",
                    transform=axes[0, 0].transAxes)
    axes[1, 0].text(-0.62, 0.5, "Streamlines", rotation=90, va="center", ha="center",
                    transform=axes[1, 0].transAxes)
    ns.colorbar(fig, cf, axes[0, -1], r"$\theta$", ticks=[0, 0.5, 1])
    ns.save(fig, out_dir, "fields_overview")


def plot_benchmark(df, out_dir):
    ref = df.dropna(subset=["Nu_avg_ref"])
    fig, (ax_n, ax_v) = plt.subplots(1, 2, figsize=(ns.TEXTWIDTH, 2.5))
    fig.subplots_adjust(wspace=0.3)
    blue = ns.CMAP_SEQUENTIAL(1.0)

    def present(colour="k"):
        return dict(color=colour, marker="o", ms=3.0, lw=0.9, zorder=3)

    def reference(colour="k"):
        return dict(ls="none", marker="s", mfc="white", mec=colour, ms=6, mew=0.7, zorder=2)

    ax_n.plot(ref["Ra"], ref["Nu_avg_ref"], **reference())
    ax_n.plot(df["Ra"], df["Nu_avg"], **present())
    ax_n.set_ylabel(r"$\overline{Nu}$")
    ns.panel_title(ax_n, "a", r"Mean Nusselt number, hot wall")

    handles = []
    for key, colour, label in (("u_max", "k", r"$u_{\max}$ on $x = L/2$"),
                               ("v_max", blue, r"$v_{\max}$ on $y = L/2$")):
        ax_v.plot(ref["Ra"], ref[f"{key}_ref"], **reference(colour))
        ax_v.plot(df["Ra"], df[key], **present(colour))
        handles.append(Line2D([], [], color=colour, lw=0.9, label=label))
    ax_v.legend(handles=handles, loc="lower right", fontsize=7.5)
    ax_v.set_ylabel(r"Peak velocity, $u L/\alpha$")
    ns.panel_title(ax_v, "b", r"Centreline velocity maxima")

    for ax in (ax_n, ax_v):
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.yaxis.set_major_formatter(ScalarFormatter())
        ax.yaxis.set_minor_formatter(NullFormatter())
        ax.set_xlabel(r"$Ra$")
    ax_n.set_yticks([1, 2, 3, 4, 5])
    ax_v.set_yticks([3, 10, 30, 100])
    key = [Line2D([], [], **present(), label="Present"),
           Line2D([], [], **reference(), label="de Vahl Davis (1983)")]
    ns.legend_above(fig, [ax_n, ax_v], key, ncol=2)
    ns.save(fig, out_dir, "benchmark")


def main():
    ns.apply_style()
    ra_list = ns.list_cases(SCENARIO, "Ra")
    if not ra_list:
        raise SystemExit(f"No output/{SCENARIO}/Ra_<N>/ folders found — run the solver "
                         "with input/scenario.txt = differentially_heated_cavity first.")
    print(f"Differentially heated cavity, cases: {ra_list}")

    cases, rows = [], []
    for ra in ra_list:
        paths = ns.case_paths(SCENARIO, "Ra", ra)
        F = ns.load_staggered(paths.fields)
        if F.T is None:
            raise SystemExit(f"{paths.fields}: temperature.vtk missing")
        plot_flow_fields(ra, F, paths)
        cases.append({"ra": ra, "F": F})

        m = pd.read_csv(os.path.join(paths.data, "dhc_metrics.csv")).iloc[0]
        M, N = F.shape
        row = {"Ra": ra, "grid": f"{N}x{M}", "u_max": m["u_max"], "v_max": m["v_max"],
               "Nu_avg": m["Nu_avg"]}
        for key, val in DVD.get(ra, {}).items():
            row[f"{key}_ref"] = val
            row[f"{key}_err_pct"] = 100 * abs(m[key] - val) / val
        rows.append(row)

    df = pd.DataFrame(rows)
    for key in ("u_max", "v_max", "Nu_avg"):
        if f"{key}_ref" not in df:
            df[f"{key}_ref"] = np.nan

    out = os.path.join(ns.scenario_dir(SCENARIO), "figures")
    print("Summary figures:")
    plot_fields_overview(cases, out)
    plot_benchmark(df, out)

    csv = os.path.join(ns.scenario_dir(SCENARIO), "validation_summary.csv")
    df.to_csv(csv, index=False, float_format="%.6g")
    print(f"  {os.path.relpath(csv, ns.PROJECT_ROOT)}\n")
    print(df.to_string(index=False, float_format=lambda v: f"{v:.4g}"))


if __name__ == "__main__":
    main()
