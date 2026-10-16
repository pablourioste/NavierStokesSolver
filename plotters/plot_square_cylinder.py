"""
plot_square_cylinder.py — Square Cylinder post-processing and validation.

Usage:
    python3 plotters/plot_square_cylinder.py [scenario]

`scenario` is the output/<scenario>/ folder (default: square_cylinder); a
sub-folder such as square_cylinder_breuer/grid_study/D20 also works. The
benchmark variant `square_cylinder_breuer` (input/square_cylinder_breuer/)
reproduces Breuer et al. (2000) and is compared against it; the default
scenario is compared, for orientation only, against Sohankar et al. (1998).
See docs/scenarios/square_cylinder_validation.md.

Writes, per case Re_<N>/figures/:
    vorticity              instantaneous vorticity (last time step)
    force_coefficients     C_D(t*) and C_L(t*)
    lift_spectrum          amplitude spectrum of C_L vs St (unsteady cases)
    recirculation          reverse-flow region and centreline velocity
                           (steady cases), with the reference L_r
and per scenario:
    figures/strouhal              St vs Re (>= 2 unsteady cases)
    figures/recirculation_length  L_r/D vs Re with the reference fit (steady cases)
    force_summary.csv             C_D, C_L, St, L_r and reference values per case

Coefficients and Re use the inflow velocity scale U of the scenario: the
uniform U_inf, or the centreline U_max of a parabolic inflow (as in Breuer).
"""
import os
import sys

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.colors import TwoSlopeNorm
from matplotlib.lines import Line2D

import nsplot as ns

# =============================================================================
# Reference data — see docs/scenarios/square_cylinder_validation.md
# =============================================================================
REFERENCES = {
    # Sohankar, Norberg & Davidson (1998), IJNMF 26:39-56, Table III:
    # uniform inflow, slip side walls, blockage 5 %. Different set-up from
    # input/square_cylinder/, so indicative only.
    "square_cylinder": {
        "label": r"Sohankar et al. (1998), $\beta = 5\,\%$",
        "unsteady": {100: (0.146, 1.460, 0.139),     # Re -> (St, mean C_D, rms C_L)
                     150: (0.162, 1.417, 0.160),
                     200: (0.165, 1.427, 0.232)},
        "steady": {},
        "lr_fit": None,
        "re_cr": 51.2,
    },
    # Breuer, Bernsdorf, Zeiser & Durst (2000), IJHFF 21:186-196: blockage 1/8,
    # parabolic inflow, Re and coefficients with U_max. Steady values at Re = 40
    # via Korukcu (2020) Table 2; L_r fit and Re_cr via Lamura et al. (2001);
    # St(100) via Chen, Xu & Li (2015).
    "square_cylinder_breuer": {
        "label": "Breuer et al. (2000)",
        "unsteady": {100: (0.137, np.nan, np.nan)},
        "steady": {40: (1.708, 2.151)},              # Re -> (C_D, L_r/D)
        "lr_fit": (-0.065, 0.0554),                  # L_r/D = a + b Re, steady range
        "re_cr": 60.0,
    },
}

SCENARIO = "square_cylinder"
REF = REFERENCES[SCENARIO]
OBST, U_INF, D = {}, 1.0, 1.0


def read_keyvals(path):
    vals = {}
    with open(path) as f:
        for line in f:
            parts = line.split("=")
            if len(parts) == 2 and not line.lstrip().startswith("#"):
                try:
                    vals[parts[0].strip()] = float(parts[1].split()[0])
                except (ValueError, IndexError):
                    pass
    return vals


def set_scenario(name):
    global SCENARIO, REF, OBST, U_INF, D
    SCENARIO = name                       # output/<name>/, may be a sub-folder
    base = name.split("/")[0]             # e.g. square_cylinder_breuer/grid_study/D20
    REF = REFERENCES.get(base, {"label": None, "unsteady": {}, "steady": {},
                                "lr_fit": None, "re_cr": None})
    inp = os.path.join(ns.PROJECT_ROOT, "input", base)
    OBST = read_keyvals(os.path.join(inp, "obstacle.txt"))
    U_INF = read_keyvals(os.path.join(inp, "sim_params.txt")).get("U_inf", 1.0)
    D = OBST["obstacle_y1"] - OBST["obstacle_y0"]    # cylinder side, reference length


def mask_solid(F):
    """Zero velocity nodes strictly inside the obstacle. Older runs left
    non-physical values there (never used by the fluid stencils); they would
    corrupt the stream function, which is integrated across the body."""
    x0, x1, y0, y1 = (OBST[k] for k in ("obstacle_x0", "obstacle_x1",
                                         "obstacle_y0", "obstacle_y1"))
    e = 1e-9
    F.u[np.ix_((F.yc > y0) & (F.yc < y1), (F.xf > x0 + e) & (F.xf < x1 - e))] = 0.0
    F.v[np.ix_((F.yf > y0 + e) & (F.yf < y1 - e), (F.xc > x0) & (F.xc < x1))] = 0.0
    return F


def title(re, F=None):
    grid = ""
    if F is not None:
        M, N = F.shape
        grid = rf", grid ${N}\times{M}$"
    setup = " (Breuer set-up)" if SCENARIO.startswith("square_cylinder_breuer") else ""
    return rf"Square cylinder{setup}, $Re = {re}${grid}"


def obstacle_patch():
    return patches.Rectangle((OBST["obstacle_x0"] / D, OBST["obstacle_y0"] / D),
                             (OBST["obstacle_x1"] - OBST["obstacle_x0"]) / D, 1.0,
                             facecolor="0.55", edgecolor="k", linewidth=0.6, zorder=4)


# =============================================================================
# Fields
# =============================================================================
def plot_vorticity(re, F, paths):
    omega = ns.vorticity(F, v_west=0.0, v_east=None) * D / U_INF
    X, Y = np.meshgrid(F.xf / D, F.yf / D)
    lim = float(f"{np.percentile(np.abs(omega), 99):.1g}")
    levels = np.linspace(-lim, lim, 21)

    Lx, Ly = F.xf[-1] / D, F.yf[-1] / D
    x_max = min(Lx, OBST["obstacle_x1"] / D + 25)       # near wake only
    width = ns.TEXTWIDTH
    fig, ax = plt.subplots(figsize=(width, width * Ly / x_max + 0.75))
    cf = ax.contourf(X, Y, np.clip(omega, -lim, lim), levels=levels,
                     cmap=ns.CMAP_DIVERGING, extend="both")
    ax.add_patch(obstacle_patch())
    ax.set_xlim(0, x_max)
    ax.set_ylim(0, Ly)
    ax.set_aspect("equal")
    ax.minorticks_off()
    ax.set_xlabel(r"$x/D$")
    ax.set_ylabel(r"$y/D$")
    ns.colorbar(fig, cf, ax, r"$\omega_z D/U$", ticks=[-lim, 0, lim])
    ax.set_title(title(re, F) + ", final time step", loc="left")
    ns.save(fig, paths.figures, "vorticity")


def centreline_u(F):
    """u on y = H/2 (mean of the two u-rows straddling it) at x = xf."""
    yc = 0.5 * (F.yf[0] + F.yf[-1])
    j = np.searchsorted(F.yc, yc)
    w = (yc - F.yc[j - 1]) / (F.yc[j] - F.yc[j - 1])
    return F.xf, (1 - w) * F.u[j - 1] + w * F.u[j]


def recirculation_length(F):
    """L_r/D: rear face of the cylinder to the first downstream zero
    crossing (negative -> positive) of u on the centreline."""
    x, u = centreline_u(F)
    x1 = OBST["obstacle_x1"]
    k = np.where(x > x1 + 1e-12)[0]
    xs, us = x[k], u[k]
    neg = np.where(us < 0)[0]
    if len(neg) == 0:
        return 0.0
    after = np.where((us[:-1] < 0) & (us[1:] >= 0))[0]
    if len(after) == 0:
        return np.nan
    i = after[0]
    xr = xs[i] - us[i] * (xs[i + 1] - xs[i]) / (us[i + 1] - us[i])
    return (xr - x1) / D


def plot_recirculation(re, F, lr, paths):
    x1 = OBST["obstacle_x1"] / D
    ref = REF["steady"].get(re)
    fig, (ax_f, ax_c) = plt.subplots(1, 2, figsize=(ns.TEXTWIDTH, 2.5),
                                     gridspec_kw={"width_ratios": [1.25, 1]})
    fig.subplots_adjust(wspace=0.62)

    # (a) streamwise velocity (reverse flow in blue) + streamlines
    xc, yc = F.xf / D, F.yc / D
    U = F.u / U_INF
    Xu, Yu = np.meshgrid(xc, yc)
    norm = TwoSlopeNorm(vmin=-0.3, vcenter=0.0, vmax=1.5)
    cf = ax_f.contourf(Xu, Yu, np.clip(U, -0.3, 1.5), levels=np.linspace(-0.3, 1.5, 19),
                       cmap=ns.CMAP_DIVERGING, norm=norm)
    ax_f.contour(Xu, Yu, U, levels=[0.0], colors="k", linewidths=0.8)
    psi = ns.stream_function(F) / (U_INF * D)
    Xp, Yp = np.meshgrid(F.xf / D, F.yf / D)
    jm = len(F.yf) // 2
    i_body = np.searchsorted(F.xf, 0.5 * (OBST["obstacle_x0"] + OBST["obstacle_x1"]))
    psi_body = psi[jm, i_body]
    ax_f.contour(Xp, Yp, psi, levels=psi_body + np.linspace(-1.2, 1.2, 13),
                 colors="k", linewidths=0.4, linestyles="solid")
    ax_f.contour(Xp, Yp, psi, levels=psi_body + np.array([-0.01, -0.003, 0.003, 0.01]),
                 colors="k", linewidths=0.4, linestyles="solid")
    ax_f.add_patch(obstacle_patch())
    ax_f.set_xlim(x1 - 3, x1 + max(4, 1.6 * lr))
    yc_mid = 0.5 * (F.yf[0] + F.yf[-1]) / D
    ax_f.set_ylim(yc_mid - 2.0, yc_mid + 2.0)
    ax_f.set_aspect("equal")
    ax_f.minorticks_off()
    ax_f.set_xlabel(r"$x/D$")
    ax_f.set_ylabel(r"$y/D$")
    ns.colorbar(fig, cf, ax_f, r"$u/U$", ticks=[-0.3, 0, 0.5, 1.0, 1.5])
    ns.panel_title(ax_f, "a", r"Streamlines, $u = 0$ contour")

    # (b) centreline velocity behind the cylinder
    x, u = centreline_u(F)
    m = x / D > x1
    ax_c.plot(x[m] / D - x1, u[m] / U_INF, color="k", lw=0.9, label="Present")
    ax_c.axhline(0, color="0.6", lw=0.4, ls=":", zorder=0)
    ax_c.plot([lr], [0], ls="none", marker="o", mfc="k", mec="k", ms=4,
              label=rf"$L_r/D = {lr:.3f}$")
    if ref is not None:
        ax_c.plot([ref[1]], [0], ls="none", marker="s", mfc="white", mec="k", ms=6,
                  mew=0.8, label=rf"{REF['label'].split(' (')[0]}, ${ref[1]:.3f}$", zorder=1)
    ax_c.set_xlim(0, max(6, 2.5 * lr))
    ax_c.set_xlabel(r"$(x - x_{\mathrm{rear}})/D$")
    ax_c.set_ylabel(r"$u/U$")
    xr = x[m] / D - x1
    shown = xr <= max(6, 2.5 * lr)
    ax_c.set_ylim(min(-0.2, 1.3 * np.min(u[m]) / U_INF),
                  max(0.6, 1.25 * np.max(u[m][shown]) / U_INF))
    ax_c.legend(loc="upper left", fontsize=7.5)
    ns.panel_title(ax_c, "b", r"Wake centreline, $y = H/2$")
    ns.heading(fig, [ax_f, ax_c], title(re, F) + ": steady recirculation")
    ns.save(fig, paths.figures, "recirculation")


# =============================================================================
# Forces
# =============================================================================
def saturated_window(df, min_periods=5, st_guess=0.15):
    """Part of the record where the lift oscillation has saturated: from the
    first time |C_L - mean| reaches 90 % of its final amplitude. Returns None
    if that leaves fewer than `min_periods` shedding periods."""
    t, cl = df["t_star"].values, df["Cl"].values
    T = t[-1]
    last = t >= 0.85 * T
    mean = cl[last].mean()
    amp = np.max(np.abs(cl[last] - mean))
    if amp <= 0:
        return None
    t_sat = t[np.argmax(np.abs(cl - mean) >= 0.9 * amp)]
    if (T - t_sat) * st_guess < min_periods:
        return None
    return df[df["t_star"] >= t_sat]


def strouhal_zero_crossing(tail):
    """St from the mean period between upward mean-crossings of C_L."""
    t, cl = tail["t_star"].values, tail["Cl"].values - tail["Cl"].mean()
    k = np.where((cl[:-1] < 0) & (cl[1:] >= 0))[0]
    if len(k) < 3:
        return np.nan
    tc = t[k] - cl[k] * (t[k + 1] - t[k]) / (cl[k + 1] - cl[k])
    return 1.0 / np.mean(np.diff(tc))


def flow_regime(df):
    """'steady' (lift at round-off level), 'developing' (lift oscillation
    still growing: wake instability not saturated) or 'periodic'."""
    t = df["t_star"].values
    cl = df["Cl"].values
    T = t[-1]
    q3 = cl[(t >= 0.50 * T) & (t < 0.75 * T)]
    q4 = cl[t >= 0.75 * T]
    a3 = np.max(np.abs(q3 - q3.mean())) if len(q3) else 0.0
    a4 = np.max(np.abs(q4 - q4.mean())) if len(q4) else 0.0
    if a4 < 1e-4:
        return "steady"
    if a3 == 0.0 or a4 / a3 > 1.2:
        return "developing"
    return "periodic"


def load_forces(paths):
    path = os.path.join(paths.data, "forces.csv")
    if not os.path.exists(path):
        return None
    df = pd.read_csv(path)
    df["t_star"] = df["time"] * U_INF / D
    # The impulsive start (fluid at rest, inflow switched on at t = 0) gives
    # a pressure spike over the first steps that dwarfs the physical range.
    return df[df["t_star"] >= 0.05 * df["t_star"].iloc[-1]].reset_index(drop=True)


def plot_forces(re, df, paths, steady):
    fig, (ax_d, ax_l) = plt.subplots(2, 1, sharex=True, figsize=(ns.TEXTWIDTH, 3.0))
    fig.subplots_adjust(hspace=0.12)
    ax_d.plot(df["t_star"], df["Cd"], color="k", lw=0.8, label="Present")
    ax_l.plot(df["t_star"], df["Cl"], color=ns.CMAP_SEQUENTIAL(1.0), lw=0.8)
    ax_l.axhline(0, color="0.6", lw=0.4, ls=":", zorder=0)
    ref = REF["steady"].get(re) if steady else None
    if ref is not None:
        ax_d.axhline(ref[0], color="0.4", lw=0.8, ls="--",
                     label=rf"{REF['label']}, $C_D = {ref[0]:.3f}$")
        ax_d.legend(loc="upper right")
    if steady:   # keep the converged value readable, not the start-up range
        tail = df["Cd"].iloc[len(df) // 4:]
        span = max(tail.max() - tail.min(), 0.05 * abs(tail.iloc[-1]))
        lo = min(tail.min(), ref[0] if ref else tail.min()) - 0.3 * span
        hi = max(tail.max(), ref[0] if ref else tail.max()) + 0.3 * span
        ax_d.set_ylim(lo, hi)
    ax_d.set_ylabel(r"$C_D$")
    ax_l.set_ylabel(r"$C_L$")
    ax_l.set_xlabel(r"$t\,U/D$")
    ax_l.set_xlim(df["t_star"].iloc[0], df["t_star"].iloc[-1])
    ax_d.set_title(title(re), loc="left")
    ns.save(fig, paths.figures, "force_coefficients")


def lift_spectrum(tail):
    """One-sided amplitude spectrum of C_L over the window `tail` (mean
    removed, Hann window). Returns (St, amplitude, St_peak)."""
    t = tail["t_star"].values
    cl = tail["Cl"].values - tail["Cl"].mean()
    if len(cl) < 16:
        return None
    dt = np.median(np.diff(t))                 # in D/U units -> frequency = St
    w = np.hanning(len(cl))
    amp = 2 * np.abs(np.fft.rfft(cl * w)) / w.sum()
    st = np.fft.rfftfreq(len(cl), dt)
    k = 1 + np.argmax(amp[1:])
    if 0 < k < len(amp) - 1:                   # parabolic peak refinement
        a, b, c = np.log(amp[k - 1:k + 2] + 1e-300)
        k = k + 0.5 * (a - c) / (a - 2 * b + c)
    return st, amp, float(np.interp(k, np.arange(len(st)), st))


def plot_lift_spectrum(re, spec, t_window, paths):
    st, amp, st_peak = spec
    fig, ax = plt.subplots(figsize=(0.55 * ns.TEXTWIDTH, 2.3))
    ax.plot(st, amp, color="k", lw=0.8, marker="o", ms=2.2, mfc="k")
    ax.axvline(st_peak, color="k", lw=0.5, ls="--", label=rf"Present, ${st_peak:.3f}$")
    if re in REF["unsteady"]:
        st_ref = REF["unsteady"][re][0]
        ax.axvline(st_ref, color="0.45", lw=0.8, ls=":",
                   label=rf"{REF['label'].split(',')[0]}, ${st_ref:.3f}$")
    ax.legend(loc="center right", title=r"Peak $St$", title_fontsize=8, fontsize=7.5,
              handlelength=1.3, handletextpad=0.4)
    ax.set_xlim(0, 0.5)
    ax.set_ylim(bottom=0)
    ax.set_xlabel(r"$St = f D/U$")
    ax.set_ylabel(r"$|\hat{C}_L|$")
    ax.set_title(rf"Lift spectrum, $Re = {re}$, ${t_window[0]:.0f} \le tU/D \le "
                 rf"{t_window[1]:.0f}$" "\n" rf"($\Delta St = {1 / (t_window[1] - t_window[0]):.3f}$)",
                 loc="left")
    ns.save(fig, paths.figures, "lift_spectrum")


# =============================================================================
# Cross-case summaries
# =============================================================================
def plot_strouhal(rows, out_dir):
    """St vs Re across the sweep; needs at least two Reynolds numbers."""
    rows = [r for r in rows if np.isfinite(r.get("St_fft", np.nan))]
    if len(rows) < 2:
        print("  figures/strouhal skipped (needs >= 2 Reynolds numbers with shedding)")
        return
    re = np.array([r["Re"] for r in rows])
    st = np.array([r["St_fft"] for r in rows])
    fig, ax = plt.subplots(figsize=(0.5 * ns.TEXTWIDTH, 2.3))
    ref_re = sorted(REF["unsteady"])
    if ref_re:
        ax.plot(ref_re, [REF["unsteady"][r][0] for r in ref_re], ls="none", marker="s",
                mfc="white", mec="k", ms=5, mew=0.7, label=REF["label"])
    ax.plot(re, st, color="k", marker="o", ms=3.5, lw=0.9, label="Present")
    lo = min([re.min()] + ref_re) - 25
    hi = max([re.max()] + ref_re) + 25
    ax.set_xlim(lo, hi)
    ax.set_ylim(min(0.10, st.min() - 0.01), max(0.20, st.max() + 0.01))
    ax.set_xlabel(r"$Re$")
    ax.set_ylabel(r"$St = f D/U$")
    ax.legend(loc="upper right")
    ns.save(fig, out_dir, "strouhal")


def plot_recirculation_length(rows, out_dir):
    rows = [r for r in rows if np.isfinite(r.get("Lr_D", np.nan))]
    if not rows or REF["lr_fit"] is None:
        return
    fig, ax = plt.subplots(figsize=(0.5 * ns.TEXTWIDTH, 2.3))
    a, b = REF["lr_fit"]
    re_fit = np.linspace(5, REF["re_cr"], 50)
    ax.plot(re_fit, a + b * re_fit, color="0.45", lw=0.9, ls="--",
            label=REF["label"] + ", fit")
    ref_re = sorted(REF["steady"])
    ax.plot(ref_re, [REF["steady"][r][1] for r in ref_re], ls="none", marker="s",
            mfc="white", mec="k", ms=5, mew=0.7, label=REF["label"])
    ax.plot([r["Re"] for r in rows], [r["Lr_D"] for r in rows], ls="none", marker="o",
            mfc="k", mec="k", ms=3.5, label="Present")
    ax.set_xlim(0, REF["re_cr"] + 5)
    ax.set_ylim(bottom=0)
    ax.set_xlabel(r"$Re$")
    ax.set_ylabel(r"$L_r/D$")
    ax.legend(loc="upper left")
    ns.save(fig, out_dir, "recirculation_length")


# =============================================================================
def main():
    set_scenario(sys.argv[1] if len(sys.argv) > 1 else "square_cylinder")
    ns.apply_style()
    re_list = ns.list_cases(SCENARIO, "Re")
    if not re_list:
        raise SystemExit(f"No output/{SCENARIO}/Re_<N>/ folders found — run the solver "
                         f"with input/scenario.txt = {SCENARIO} first.")
    print(f"{SCENARIO}, cases: {re_list}")

    rows = []
    for re in re_list:
        paths = ns.case_paths(SCENARIO, "Re", re)
        row = {"Re": re}
        F = None
        if os.path.exists(os.path.join(paths.fields, "u_velocity.vtk")):
            F = mask_solid(ns.load_staggered(paths.fields))
            plot_vorticity(re, F, paths)
        else:
            print(f"  Re_{re}: no fields/ (run did not finish) — field plots skipped")

        df = load_forces(paths)
        if df is not None and len(df) > 1:
            regime = flow_regime(df)
            steady = regime == "steady"
            tail = saturated_window(df) if regime == "periodic" else None
            if regime == "periodic" and tail is None:
                regime = "developing"             # saturated, but too few periods
            if tail is None:
                tail = df[df["t_star"] >= 0.5 * df["t_star"].iloc[-1]]
            cl_rms = np.sqrt(np.mean((tail["Cl"] - tail["Cl"].mean()) ** 2))
            row.update({"t_star_end": df["t_star"].iloc[-1], "regime": regime,
                        "t_star_stats_from": tail["t_star"].iloc[0],
                        "Cd_mean": tail["Cd"].mean(), "Cl_mean": tail["Cl"].mean(),
                        "Cl_rms": cl_rms})
            if steady:
                row["Cd_mean"] = df["Cd"].iloc[-1]       # converged value
            plot_forces(re, df, paths, steady)
            spec = lift_spectrum(tail) if regime == "periodic" else None
            if spec is not None:
                plot_lift_spectrum(re, spec, (tail["t_star"].iloc[0], tail["t_star"].iloc[-1]),
                                   paths)
                row["St_fft"] = spec[2]
                row["St_zc"] = strouhal_zero_crossing(tail)
            if steady and F is not None:
                row["Lr_D"] = recirculation_length(F)
                plot_recirculation(re, F, row["Lr_D"], paths)

        st_path = os.path.join(paths.data, "strouhal.csv")
        if os.path.exists(st_path):
            st = float(pd.read_csv(st_path).iloc[0]["St"])
            # 0 = no shedding detected; steady runs only give noise crossings
            row["St_solver"] = st if st > 0 and row.get("regime") != "steady" else np.nan
        if re in REF["unsteady"]:
            row.update(dict(zip(["St_ref", "Cd_mean_ref", "Cl_rms_ref"], REF["unsteady"][re])))
        if re in REF["steady"]:
            row.update(dict(zip(["Cd_ref", "Lr_D_ref"], REF["steady"][re])))
        rows.append(row)

    out = os.path.join(ns.scenario_dir(SCENARIO), "figures")
    plot_strouhal(rows, out)
    plot_recirculation_length(rows, out)
    df = pd.DataFrame(rows)
    csv = os.path.join(ns.scenario_dir(SCENARIO), "force_summary.csv")
    df.to_csv(csv, index=False, float_format="%.6g")
    print(f"  {os.path.relpath(csv, ns.PROJECT_ROOT)}  (periodic: statistics over the "
          "saturated part of the record, from t_star_stats_from; steady: final values)\n")
    print(df.to_string(index=False, float_format=lambda v: f"{v:.4g}"))


if __name__ == "__main__":
    main()
