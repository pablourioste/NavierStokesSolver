"""
plot_results.py — Lid-Driven Cavity post-processing and validation.

Reads output/lid_driven_cavity/Re_<N>/ and writes

  per case   Re_<N>/figures/flow_fields.{pdf,png}   (a) stream function,
                                                     (b) vorticity, (c) pressure
             Re_<N>/figures/centerlines.{pdf,png}   u(y) at x = L/2, v(x) at y = H/2
                                                     (+ Ghia et al. when available)
  summary    figures/centerline_comparison          all Ghia Reynolds numbers overlaid
             figures/centerline_trend               every computed Re, coloured by Re
             figures/streamfunction_overview        flow topology across Re
             validation_summary.csv                 errors vs. Ghia et al. (1982)

Reference: U. Ghia, K. N. Ghia, C. T. Shin, J. Comput. Phys. 48 (1982) 387-411.
Run from anywhere:  python3 plotters/plot_results.py
"""
import os

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import BoundaryNorm, LogNorm, Normalize
from matplotlib.lines import Line2D
from matplotlib.cm import ScalarMappable

import nsplot as ns

SCENARIO = "lid_driven_cavity"
U_LID = 1.0

# =============================================================================
# Reference data — Ghia, Ghia & Shin (1982)
# =============================================================================
GHIA_RE = [100, 400, 1000, 3200, 5000, 7500, 10000]
# Table I: y/H, then u/U for each Re in GHIA_RE
GHIA_U = np.array([
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
# Table II: x/L, then v/U for each Re in GHIA_RE
GHIA_V = np.array([
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
# Table V: primary vortex — psi_min, x, y
GHIA_PRIMARY = {
    100:   (-0.103423, 0.6172, 0.7344),
    400:   (-0.113909, 0.5547, 0.6055),
    1000:  (-0.117929, 0.5313, 0.5625),
    3200:  (-0.120377, 0.5165, 0.5469),
    5000:  (-0.118966, 0.5117, 0.5352),
    7500:  (-0.119976, 0.5117, 0.5322),
    10000: (-0.119731, 0.5117, 0.5333),
}
# Stream-function contour levels used in Ghia et al. (1982), Table III
PSI_LEVELS_PRIMARY = [-0.1175, -0.115, -0.11, -0.1, -0.09, -0.07, -0.05, -0.03,
                      -0.01, -1e-4, -1e-5, -1e-7, -1e-10]
PSI_LEVELS_SECONDARY = [1e-8, 1e-7, 1e-6, 1e-5, 5e-5, 1e-4, 2.5e-4, 5e-4,
                        1e-3, 1.5e-3, 3e-3]
OMEGA_LEVELS = [-5, -4, -3, -2, -1, -0.5, 0, 0.5, 1, 2, 3, 4, 5]

SECONDARY_COLOUR = "#B2182B"   # counter-rotating (psi > 0) eddies

# Cross-case figures only show runs with at least this many cells per
# direction: under-resolved grids (e.g. 50x50 at Re >= 1000) would otherwise be
# plotted next to converged ones. Per-case figures and the summary CSV still
# list every case together with its grid.
MIN_CELLS_SUMMARY = 100


def ghia(re):
    """(y, u), (x, v) reference arrays, or None if Re is not tabulated."""
    if re not in GHIA_RE:
        return None
    c = 1 + GHIA_RE.index(re)
    return (GHIA_U[:, 0], GHIA_U[:, c]), (GHIA_V[:, 0], GHIA_V[:, c])


# =============================================================================
# Loading
# =============================================================================
def load_profiles(paths):
    """Centreline profiles with the wall values (boundary conditions) appended."""
    du = pd.read_csv(os.path.join(paths.data, "centerline_u.csv"))
    dv = pd.read_csv(os.path.join(paths.data, "centerline_v.csv"))
    y = np.concatenate([[0.0], du["y_norm"].values, [1.0]])
    u = np.concatenate([[0.0], du["u_norm"].values, [U_LID]])
    x = np.concatenate([[0.0], dv["x_norm"].values, [1.0]])
    v = np.concatenate([[0.0], dv["v_norm"].values, [0.0]])
    return (y, u), (x, v)


def grid_label(F):
    M, N = F.shape
    return rf"${N}\times{M}$"


# =============================================================================
# Per-case figures
# =============================================================================
def draw_streamfunction(ax, F, psi, mark_centre=True):
    X, Y = np.meshgrid(F.xf, F.yf)
    ax.contour(X, Y, psi, levels=PSI_LEVELS_PRIMARY, colors="k",
               linewidths=0.55, linestyles="solid")
    ax.contour(X, Y, psi, levels=PSI_LEVELS_SECONDARY, colors=SECONDARY_COLOUR,
               linewidths=0.45, linestyles="solid")
    if mark_centre:
        _, xc, yc = ns.vortex_centre(psi, F.xf, F.yf, sign=-1)
        ax.plot(xc, yc, marker="+", color="k", ms=5, mew=0.8)
    domain_frame(ax, F)


def domain_frame(ax, F):
    ax.set_xlim(F.xf[0], F.xf[-1])
    ax.set_ylim(F.yf[0], F.yf[-1])
    ns.unit_ticks(ax)
    ns.plain_axes(ax)


def plot_flow_fields(re, F, paths):
    psi = ns.stream_function(F)
    omega = ns.vorticity(F, u_north=U_LID)

    fig, axes = plt.subplots(1, 3, figsize=(ns.TEXTWIDTH, 2.0))
    fig.subplots_adjust(wspace=0.55)

    ax = axes[0]
    draw_streamfunction(ax, F, psi)
    ns.panel_title(ax, "a", r"Stream function $\psi$")

    ax = axes[1]
    X, Y = np.meshgrid(F.xf, F.yf)
    norm = BoundaryNorm(OMEGA_LEVELS, ncolors=256, extend="both")
    cf = ax.contourf(X, Y, omega, levels=OMEGA_LEVELS, cmap=ns.CMAP_DIVERGING,
                     norm=norm, extend="both")
    ax.contour(X, Y, omega, levels=OMEGA_LEVELS, colors="k", linewidths=0.3)
    domain_frame(ax, F)
    ns.colorbar(fig, cf, ax, r"$\omega_z$", ticks=[-5, -3, -1, 0, 1, 3, 5])
    ns.panel_title(ax, "b", r"Vorticity $\omega_z$")

    ax = axes[2]
    p = F.p - F.p.mean()
    lim = np.percentile(np.abs(p), 98)
    levels = np.linspace(-lim, lim, 21)
    Xc, Yc = np.meshgrid(F.xc, F.yc)
    cf = ax.contourf(Xc, Yc, np.clip(p, -lim, lim), levels=levels,
                     cmap=ns.CMAP_DIVERGING, extend="both")
    ax.contour(Xc, Yc, p, levels=levels[::2], colors="k", linewidths=0.3)
    domain_frame(ax, F)
    tick = float(f"{lim * 0.9:.1g}")
    ns.colorbar(fig, cf, ax, r"$p - \langle p \rangle$", ticks=[-tick, 0, tick])
    ns.panel_title(ax, "c", r"Pressure")

    for ax in axes:
        ax.set_xlabel(r"$x/L$")
    axes[0].set_ylabel(r"$y/H$")
    for ax in axes[1:]:
        ax.set_yticklabels([])

    ns.heading(fig, axes, rf"Lid-driven cavity, $Re = {re}$, grid {grid_label(F)}")
    ns.save(fig, paths.figures, "flow_fields")


def style_profile_axes(ax_u, ax_v):
    ax_u.axvline(0, color="0.6", lw=0.4, ls=":", zorder=0)
    ax_v.axhline(0, color="0.6", lw=0.4, ls=":", zorder=0)
    ax_u.set_xlabel(r"$u/U$")
    ax_u.set_ylabel(r"$y/H$")
    ax_u.set_ylim(0, 1)
    ax_v.set_xlabel(r"$x/L$")
    ax_v.set_ylabel(r"$v/U$")
    ax_v.set_xlim(0, 1)


def plot_centerlines(re, F, profiles, paths):
    (y, u), (x, v) = profiles
    ref = ghia(re)

    fig, (ax_u, ax_v) = plt.subplots(1, 2, figsize=(ns.TEXTWIDTH, 2.6))
    fig.subplots_adjust(wspace=0.32)
    ax_u.plot(u, y, color="k", label=rf"Present, {grid_label(F)}")
    ax_v.plot(x, v, color="k")
    if ref is not None:
        (gy, gu), (gx, gv) = ref
        kw = dict(ls="none", marker="o", mfc="white", mec="k", ms=4, zorder=3)
        ax_u.plot(gu, gy, label="Ghia et al. (1982)", **kw)
        ax_v.plot(gx, gv, **kw)
    style_profile_axes(ax_u, ax_v)
    ax_u.legend(loc="lower right")
    ns.panel_title(ax_u, "a", r"$u$ along $x = L/2$")
    ns.panel_title(ax_v, "b", r"$v$ along $y = H/2$")
    ns.heading(fig, [ax_u, ax_v], rf"Lid-driven cavity, $Re = {re}$")
    ns.save(fig, paths.figures, "centerlines")


# =============================================================================
# Cross-case figures
# =============================================================================
def re_colours(re_values):
    norm = LogNorm(vmin=min(re_values), vmax=max(re_values)) if len(re_values) > 1 \
        else Normalize(vmin=0, vmax=1)
    return {re: ns.CMAP_SEQUENTIAL(norm(re) if len(re_values) > 1 else 1.0)
            for re in re_values}, norm


def plot_centerline_comparison(cases, out_dir):
    """Ghia-tabulated Re only: line = present, open symbol = Ghia."""
    sel = [c for c in cases if c["re"] in GHIA_RE]
    if not sel:
        return
    colours, _ = re_colours([c["re"] for c in sel])

    fig, (ax_u, ax_v) = plt.subplots(1, 2, figsize=(ns.TEXTWIDTH, 3.0))
    fig.subplots_adjust(wspace=0.32)
    handles = []
    for k, c in enumerate(sel):
        re, col, mk = c["re"], colours[c["re"]], ns.MARKERS[k % len(ns.MARKERS)]
        (y, u), (x, v) = c["profiles"]
        (gy, gu), (gx, gv) = ghia(re)
        ax_u.plot(u, y, color=col, lw=0.9)
        ax_v.plot(x, v, color=col, lw=0.9)
        kw = dict(ls="none", marker=mk, mfc="white", mec=col, ms=3.6, mew=0.7, zorder=3)
        ax_u.plot(gu, gy, **kw)
        ax_v.plot(gx, gv, **kw)
        handles.append(Line2D([], [], color=col, lw=0.9, marker=mk, mfc="white",
                              mec=col, ms=3.6, mew=0.7, label=rf"$Re = {re}$"))
    style_profile_axes(ax_u, ax_v)
    ns.panel_title(ax_u, "a", r"$u$ along $x = L/2$")
    ns.panel_title(ax_v, "b", r"$v$ along $y = H/2$")

    key = [Line2D([], [], color="k", lw=0.9, label="Present"),
           Line2D([], [], ls="none", marker="o", mfc="white", mec="k", ms=3.6,
                  label="Ghia et al. (1982)")]
    ns.legend_above(fig, [ax_u, ax_v], handles + key, ncol=min(4, len(handles) + 2))
    ns.save(fig, out_dir, "centerline_comparison")


def plot_centerline_trend(cases, out_dir):
    """Every computed Re (no reference), coloured continuously by Re."""
    res = [c["re"] for c in cases]
    if len(res) < 2:
        return
    colours, norm = re_colours(res)

    fig, (ax_u, ax_v) = plt.subplots(1, 2, figsize=(ns.TEXTWIDTH, 2.6))
    fig.subplots_adjust(wspace=0.32, right=0.86)
    for c in cases:
        (y, u), (x, v) = c["profiles"]
        ax_u.plot(u, y, color=colours[c["re"]], lw=0.8)
        ax_v.plot(x, v, color=colours[c["re"]], lw=0.8)
    style_profile_axes(ax_u, ax_v)
    ns.panel_title(ax_u, "a", r"$u$ along $x = L/2$")
    ns.panel_title(ax_v, "b", r"$v$ along $y = H/2$")

    cax = fig.add_axes([0.89, 0.11, 0.015, 0.77])
    cb = fig.colorbar(ScalarMappable(norm=norm, cmap=ns.CMAP_SEQUENTIAL), cax=cax)
    cb.set_label(r"$Re$")
    cb.set_ticks([r for r in res if r in GHIA_RE] or res)
    cb.set_ticklabels([str(r) for r in (cb.get_ticks())])
    cb.ax.minorticks_off()
    cb.ax.tick_params(direction="in", width=0.5, length=2.5)
    cb.outline.set_linewidth(0.5)
    ns.save(fig, out_dir, "centerline_trend")


def plot_streamfunction_overview(cases, out_dir):
    sel = [c for c in cases if c["re"] in GHIA_RE] or cases
    sel = sel[:6]
    ncol = min(3, len(sel))
    nrow = int(np.ceil(len(sel) / ncol))
    w = ns.TEXTWIDTH
    fig, axes = plt.subplots(nrow, ncol, figsize=(w, w / ncol * nrow * 1.02),
                             squeeze=False)
    fig.subplots_adjust(wspace=0.12, hspace=0.22)
    for k, ax in enumerate(axes.flat):
        if k >= len(sel):
            ax.set_visible(False)
            continue
        c = sel[k]
        draw_streamfunction(ax, c["F"], ns.stream_function(c["F"]))
        ax.set_title(rf"$Re = {c['re']}$")
        r, q = divmod(k, ncol)
        if q == 0:
            ax.set_ylabel(r"$y/H$")
        else:
            ax.set_yticklabels([])
        if r == nrow - 1 or k + ncol >= len(sel):
            ax.set_xlabel(r"$x/L$")
        else:
            ax.set_xticklabels([])
    key = [Line2D([], [], color="k", lw=0.6, label=r"Primary vortex, $\psi < 0$"),
           Line2D([], [], color=SECONDARY_COLOUR, lw=0.6,
                  label=r"Secondary eddies, $\psi > 0$"),
           Line2D([], [], ls="none", marker="+", color="k", ms=5, mew=0.8,
                  label=r"Vortex centre, $\psi_{\min}$")]
    ns.legend_above(fig, axes, key, ncol=3)
    ns.save(fig, out_dir, "streamfunction_overview")


# =============================================================================
# Quantitative validation
# =============================================================================
def validation_row(c):
    re, F = c["re"], c["F"]
    M, N = F.shape
    psi_min, xc, yc = ns.vortex_centre(ns.stream_function(F), F.xf, F.yf, sign=-1)
    row = {"Re": re, "grid": f"{N}x{M}", "psi_min": psi_min, "x_c": xc, "y_c": yc}
    ref = ghia(re)
    if ref is not None:
        (y, u), (x, v) = c["profiles"]
        (gy, gu), (gx, gv) = ref
        iu = (gy > 0) & (gy < 1)
        iv = (gx > 0) & (gx < 1)
        eu = np.interp(gy[iu], y, u) - gu[iu]
        ev = np.interp(gx[iv], x, v) - gv[iv]
        g_psi, g_x, g_y = GHIA_PRIMARY[re]
        row.update({
            "psi_min_ghia": g_psi, "x_c_ghia": g_x, "y_c_ghia": g_y,
            "psi_min_rel_err_pct": 100 * abs(psi_min - g_psi) / abs(g_psi),
            "u_rms_err": np.sqrt(np.mean(eu ** 2)), "u_max_err": np.max(np.abs(eu)),
            "v_rms_err": np.sqrt(np.mean(ev ** 2)), "v_max_err": np.max(np.abs(ev)),
        })
    return row


def main():
    ns.apply_style()
    re_list = ns.list_cases(SCENARIO, "Re")
    if not re_list:
        raise SystemExit(f"No output/{SCENARIO}/Re_<N>/ folders found — run the solver first.")
    print(f"Lid-driven cavity, cases: {re_list}")

    cases = []
    for re in re_list:
        paths = ns.case_paths(SCENARIO, "Re", re)
        F = ns.load_staggered(paths.fields)
        profiles = load_profiles(paths)
        plot_flow_fields(re, F, paths)
        plot_centerlines(re, F, profiles, paths)
        cases.append({"re": re, "F": F, "profiles": profiles})

    out = os.path.join(ns.scenario_dir(SCENARIO), "figures")
    resolved = [c for c in cases if min(c["F"].shape) >= MIN_CELLS_SUMMARY]
    skipped = [c["re"] for c in cases if c not in resolved]
    print(f"Summary figures (grids >= {MIN_CELLS_SUMMARY} cells; omitted Re = {skipped}):")
    plot_centerline_comparison(resolved, out)
    plot_centerline_trend(resolved, out)
    plot_streamfunction_overview(resolved, out)

    df = pd.DataFrame([validation_row(c) for c in cases])
    csv = os.path.join(ns.scenario_dir(SCENARIO), "validation_summary.csv")
    df.to_csv(csv, index=False, float_format="%.6g")
    print(f"  {os.path.relpath(csv, ns.PROJECT_ROOT)}\n")
    cols = [c for c in ["Re", "grid", "psi_min", "psi_min_ghia", "psi_min_rel_err_pct",
                        "u_rms_err", "v_rms_err"] if c in df]
    print(df[cols].to_string(index=False, float_format=lambda v: f"{v:.4f}"))


if __name__ == "__main__":
    main()
