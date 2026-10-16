"""
nsplot.py — shared figure style and I/O for every plotter in this folder.

Style
    apply_style()   LaTeX (Computer Modern) typesetting, journal-style axes.
                    Falls back to matplotlib's Computer Modern mathtext if no
                    LaTeX installation is found.
    save(fig, ...)  writes <name>.pdf (vector, for LaTeX documents) and
                    <name>.png (300 dpi, for the README / slides).

Output layout (written by the solver, see src/main.cpp)
    output/<scenario>/mesh/                    mesh_P/U/V.vtk, obstacle_mask.vtk
    output/<scenario>/<Re_N|Ra_N>/fields/      VTK fields
    output/<scenario>/<Re_N|Ra_N>/data/        CSV profiles / metrics
    output/<scenario>/<Re_N|Ra_N>/figures/     per-case figures (plotters)
    output/<scenario>/figures/                 cross-case figures (plotters)

Staggered fields
    load_staggered() reads pressure/u/v(/temperature) VTK files and places each
    unknown at its true MAC-grid location (u on vertical faces, v on horizontal
    faces, p and T at cell centres). stream_function() and vorticity() are the
    exact discrete operators on that grid, so no interpolation is involved.
"""
import os
import re
import shutil
from dataclasses import dataclass
from typing import Optional

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap

PLOTTERS_DIR = os.path.dirname(os.path.abspath(__file__))
PROJECT_ROOT = os.path.dirname(PLOTTERS_DIR)
OUTPUT_ROOT = os.path.join(PROJECT_ROOT, "output")

# Width of the text block of an A4 thesis with ~2.5 cm margins. Figures are
# drawn at this physical size so that, included with width=\textwidth, the
# fonts come out at the document's own 9-10 pt size (no rescaling).
TEXTWIDTH = 6.3  # in

# Diverging map for signed fields (vorticity, pressure, lift): two hues and
# a neutral midpoint.
CMAP_DIVERGING = "RdBu_r"


def _truncate(cmap_name, lo, hi, n=256):
    base = plt.get_cmap(cmap_name)
    return LinearSegmentedColormap.from_list(
        f"{cmap_name}_{lo:.2f}_{hi:.2f}", base(np.linspace(lo, hi, n)))


# Single-hue sequential map for ordered parameters (Re, Ra). The palest end is
# cut off so that every curve keeps enough contrast against a white page.
CMAP_SEQUENTIAL = _truncate("Blues", 0.38, 1.0)
CMAP_THERMAL = "RdBu_r"   # hot (red) wall on the left, cold (blue) on the right

# Marker shapes cycled alongside the sequential colour, so series stay
# distinguishable in grayscale print and for colour-blind readers.
MARKERS = ["o", "s", "^", "D", "v", "p", "h", "<", ">"]


def apply_style():
    has_latex = shutil.which("latex") is not None and shutil.which("dvipng") is not None
    rc = {
        "font.family": "serif",
        "font.size": 9,
        "axes.labelsize": 9,
        "axes.titlesize": 9,
        "legend.fontsize": 8,
        "xtick.labelsize": 8,
        "ytick.labelsize": 8,
        "axes.linewidth": 0.6,
        "axes.titlepad": 4,
        "axes.labelpad": 2,
        "lines.linewidth": 1.0,
        "lines.markersize": 3.5,
        "lines.markeredgewidth": 0.6,
        "xtick.direction": "in",
        "ytick.direction": "in",
        "xtick.top": True,
        "ytick.right": True,
        "xtick.minor.visible": True,
        "ytick.minor.visible": True,
        "xtick.major.size": 3.0,
        "ytick.major.size": 3.0,
        "xtick.minor.size": 1.6,
        "ytick.minor.size": 1.6,
        "xtick.major.width": 0.6,
        "ytick.major.width": 0.6,
        "xtick.minor.width": 0.4,
        "ytick.minor.width": 0.4,
        "xtick.major.pad": 2.5,
        "ytick.major.pad": 2.5,
        "axes.grid": False,
        "legend.frameon": False,
        "legend.handlelength": 1.8,
        "legend.borderaxespad": 0.3,
        "figure.dpi": 150,
        "savefig.dpi": 300,
        "savefig.bbox": "tight",
        "savefig.pad_inches": 0.02,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        "axes.formatter.use_mathtext": True,
        "axes.unicode_minus": True,
    }
    if has_latex:
        rc.update({
            "text.usetex": True,
            "font.serif": ["Computer Modern Roman"],
            "text.latex.preamble": r"\usepackage{amsmath}\usepackage{bm}",
        })
    else:
        rc.update({
            "text.usetex": False,
            "mathtext.fontset": "cm",
            "font.serif": ["CMU Serif", "Latin Modern Roman", "DejaVu Serif"],
        })
    mpl.rcParams.update(rc)


def save(fig, out_dir, name):
    os.makedirs(out_dir, exist_ok=True)
    stem = os.path.join(out_dir, name)
    fig.savefig(stem + ".pdf")
    fig.savefig(stem + ".png", dpi=300)
    plt.close(fig)
    print(f"  {os.path.relpath(stem, PROJECT_ROOT)}.{{pdf,png}}")


def plain_axes(ax):
    """Field plots: no minor ticks, ticks only at domain bounds and centre."""
    ax.minorticks_off()
    ax.set_aspect("equal")


def colorbar(fig, mappable, ax, label, ticks=None):
    """Colour bar exactly as tall as `ax`, placed just to its right."""
    ax.apply_aspect()
    w_in = ax.get_position().width * fig.get_figwidth()
    cax = ax.inset_axes([1.0 + 0.06 / w_in, 0.0, 0.075 / w_in, 1.0])
    cb = fig.colorbar(mappable, cax=cax, ticks=ticks)
    cb.set_label(label)
    cb.ax.tick_params(which="both", direction="in", width=0.5, length=2.5, right=True, left=False)
    cb.ax.minorticks_off()
    cb.outline.set_linewidth(0.5)
    return cb


def unit_ticks(ax):
    """Ticks at 0, 1/2, 1 labelled without trailing zeros (no clash at the
    origin corner)."""
    ax.set_xticks([0, 0.5, 1], ["0", "0.5", "1"])
    ax.set_yticks([0, 0.5, 1], ["0", "0.5", "1"])


def top_of(fig, axes):
    """Top edge (figure fraction) of a set of axes after aspect is applied."""
    tops = []
    for ax in np.ravel(axes):
        if ax.get_visible():
            ax.apply_aspect()
            tops.append(ax.get_position().y1)
    return max(tops)


def heading(fig, axes, text, gap=0.30):
    """Figure heading placed `gap` inches above the top of `axes` (clears
    panel titles), instead of matplotlib's fixed suptitle position."""
    fig.text(0.5, top_of(fig, axes) + gap / fig.get_figheight(), text,
             ha="center", va="bottom", fontsize=9)


def legend_above(fig, axes, handles, ncol, gap=0.22):
    """Figure legend centred `gap` inches above the top of `axes`."""
    return fig.legend(handles=handles, loc="lower center", ncol=ncol,
                      bbox_to_anchor=(0.5, top_of(fig, axes) + gap / fig.get_figheight()),
                      columnspacing=1.4)


def panel_title(ax, letter, text):
    ax.set_title(rf"({letter}) {text}", loc="left")


def sci_tex(x):
    """1e4 -> 10^{4}, 2.5e5 -> 2.5\\times10^{5} (LaTeX, no $)."""
    exp = int(np.floor(np.log10(abs(x))))
    mant = x / 10 ** exp
    if abs(mant - 1.0) < 1e-9:
        return rf"10^{{{exp}}}"
    return rf"{mant:g}\times 10^{{{exp}}}"


# =============================================================================
# Case discovery
# =============================================================================
def scenario_dir(scenario):
    return os.path.join(OUTPUT_ROOT, scenario)


def list_cases(scenario, prefix):
    """Sorted numeric values of output/<scenario>/<prefix>_<N>/ folders."""
    root = scenario_dir(scenario)
    if not os.path.isdir(root):
        return []
    vals = []
    for name in os.listdir(root):
        m = re.fullmatch(rf"{prefix}_(\d+)", name)
        if m and os.path.isdir(os.path.join(root, name)):
            vals.append(int(m.group(1)))
    return sorted(vals)


@dataclass
class CasePaths:
    root: str

    @property
    def fields(self):
        return os.path.join(self.root, "fields")

    @property
    def data(self):
        return os.path.join(self.root, "data")

    @property
    def figures(self):
        return os.path.join(self.root, "figures")


def case_paths(scenario, prefix, value):
    return CasePaths(os.path.join(scenario_dir(scenario), f"{prefix}_{value}"))


# =============================================================================
# VTK I/O
# =============================================================================
def read_vtk(path):
    """ASCII RECTILINEAR_GRID with one CELL_DATA scalar, as written by the
    solver. Returns (x_edges, y_edges, field[ny_cells, nx_cells])."""
    with open(path) as f:
        tokens = f.read().split()
    coords = {}
    values = None
    i = 0
    while i < len(tokens):
        t = tokens[i]
        if t in ("X_COORDINATES", "Y_COORDINATES", "Z_COORDINATES"):
            n = int(tokens[i + 1])
            coords[t[0]] = np.array(tokens[i + 3:i + 3 + n], dtype=float)
            i += 3 + n
            continue
        if t == "LOOKUP_TABLE":
            values = np.array(tokens[i + 2:], dtype=float)
            break
        i += 1
    xe, ye = coords["X"], coords["Y"]
    nx, ny = len(xe) - 1, len(ye) - 1
    if values is None or values.size < nx * ny:
        raise ValueError(f"{path}: expected {nx * ny} cell values")
    return xe, ye, values[:nx * ny].reshape(ny, nx)


@dataclass
class StaggeredField:
    xf: np.ndarray            # P-cell faces in x  (N+1)
    yf: np.ndarray            # P-cell faces in y  (M+1)
    xc: np.ndarray            # P-cell centres     (N)
    yc: np.ndarray            # P-cell centres     (M)
    p: np.ndarray             # (M, N)     at (xc, yc)
    u: np.ndarray             # (M, N+1)   at (xf, yc)
    v: np.ndarray             # (M+1, N)   at (xc, yf)
    T: Optional[np.ndarray] = None  # (M, N) at (xc, yc)

    @property
    def shape(self):
        return len(self.yc), len(self.xc)

    def u_centre(self):
        return 0.5 * (self.u[:, :-1] + self.u[:, 1:])

    def v_centre(self):
        return 0.5 * (self.v[:-1, :] + self.v[1:, :])


def load_staggered(fields_dir):
    xf, yf, p = read_vtk(os.path.join(fields_dir, "pressure.vtk"))
    _, _, u = read_vtk(os.path.join(fields_dir, "u_velocity.vtk"))
    _, _, v = read_vtk(os.path.join(fields_dir, "v_velocity.vtk"))
    M, N = p.shape
    if u.shape != (M, N + 1) or v.shape != (M + 1, N):
        raise ValueError(f"{fields_dir}: unexpected staggered shapes "
                         f"p{p.shape} u{u.shape} v{v.shape}")
    T = None
    t_path = os.path.join(fields_dir, "temperature.vtk")
    if os.path.exists(t_path):
        T = read_vtk(t_path)[2]
    return StaggeredField(xf, yf, 0.5 * (xf[:-1] + xf[1:]), 0.5 * (yf[:-1] + yf[1:]),
                          p, u, v, T)


# =============================================================================
# Derived quantities on the MAC grid
# =============================================================================
def stream_function(F):
    """psi at P-cell corners (M+1, N+1), with u = dpsi/dy, v = -dpsi/dx and
    psi(x0, y0) = 0. Integrates -v along the bottom row, then u up every
    vertical face line — exact for a discretely divergence-free field."""
    dx = np.diff(F.xf)
    dy = np.diff(F.yf)
    M, N = F.shape
    psi = np.zeros((M + 1, N + 1))
    psi[0, 1:] = -np.cumsum(F.v[0, :] * dx)
    psi[1:, :] = psi[0, :] + np.cumsum(F.u * dy[:, None], axis=0)
    return psi


def vorticity(F, u_south=0.0, u_north=0.0, v_west=0.0, v_east=0.0):
    """omega_z = dv/dx - du/dy at P-cell corners (M+1, N+1). The wall
    tangential velocities close the one-sided differences on the boundary
    (same image-ghost treatment as the solver); None means an open boundary
    (zero normal gradient, e.g. an outlet)."""
    M, N = F.shape
    v_west = F.v[:, :1] if v_west is None else np.full((M + 1, 1), v_west)
    v_east = F.v[:, -1:] if v_east is None else np.full((M + 1, 1), v_east)
    # dv/dx on vertical face lines, v from cell centres; walls use v_west/east
    xs = np.concatenate([[F.xf[0]], F.xc, [F.xf[-1]]])
    vs = np.concatenate([v_west, F.v, v_east], axis=1)
    dvdx = np.diff(vs, axis=1) / np.diff(xs)[None, :]
    ys = np.concatenate([[F.yf[0]], F.yc, [F.yf[-1]]])
    us = np.concatenate([np.full((1, N + 1), u_south), F.u, np.full((1, N + 1), u_north)], axis=0)
    dudy = np.diff(us, axis=0) / np.diff(ys)[:, None]
    return dvdx - dudy


def vortex_centre(psi, xf, yf, sign=-1):
    """Location and value of the extremum of psi (min for sign=-1), refined
    with a 2D quadratic fit over the 3x3 neighbourhood of the extremal node."""
    j, i = np.unravel_index(np.argmax(sign * psi), psi.shape)
    if 0 < i < len(xf) - 1 and 0 < j < len(yf) - 1:
        X, Y = np.meshgrid(xf[i - 1:i + 2], yf[j - 1:j + 2])
        Z = psi[j - 1:j + 2, i - 1:i + 2]
        A = np.column_stack([np.ones(9), X.ravel(), Y.ravel(), X.ravel() ** 2,
                             X.ravel() * Y.ravel(), Y.ravel() ** 2])
        c = np.linalg.lstsq(A, Z.ravel(), rcond=None)[0]
        H = np.array([[2 * c[3], c[4]], [c[4], 2 * c[5]]])
        try:
            x0, y0 = np.linalg.solve(H, -c[1:3])
            if xf[i - 1] <= x0 <= xf[i + 1] and yf[j - 1] <= y0 <= yf[j + 1]:
                val = c @ [1, x0, y0, x0 ** 2, x0 * y0, y0 ** 2]
                return float(val), float(x0), float(y0)
        except np.linalg.LinAlgError:
            pass
    return float(psi[j, i]), float(xf[i]), float(yf[j])
