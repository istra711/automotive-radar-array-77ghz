"""Write the VisualopenEMS project files and the construction previews.

Usage::

    python make_project.py                # all four .fcsim.json files
    python make_project.py --preview      # additionally render the PNGs
    python make_project.py --only broadside

The four files differ only in the excitation (see ``array_model.VARIANTS``):

=========================  ==========================================
patch_element_77ghz        one patch at the board centre
patch_array_4x4_77ghz      16 patches, P6 excited (embedded element)
patch_array_4x4_broadside  all 16 ports in phase
patch_array_4x4_scan30     all 16 ports, progressive phase, +30 deg
=========================  ==========================================

The JSON is human readable and loads in VisualopenEMS via ``File > Open``.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import array_model as am

HERE = Path(__file__).resolve().parent
FIG = HERE / "docs" / "figures"

#: variant -> output file stem
PROJECTS = {
    "element": "patch_element_77ghz",
    "embedded": "patch_array_4x4_77ghz",
    "broadside": "patch_array_4x4_broadside",
    "scan30": "patch_array_4x4_scan30",
}


def write_project(variant: str) -> Path:
    stem = PROJECTS[variant]
    data = am.project_dict(name=stem, variant=variant)
    path = HERE / f"{stem}.fcsim.json"
    path.write_text(json.dumps(data, indent=2, ensure_ascii=False) + "\n",
                    encoding="utf-8")
    return path


# ---------------------------------------------------------------------------
# Previews
# ---------------------------------------------------------------------------

_CU = "#c98a2e"
_CU_EDGE = "#6b4410"
_SUB = "#cfe0f2"
_SUB_EDGE = "#315a86"


def _dims(ax, x0: float, y0: float, x1: float, y1: float, text: str,
          color: str = "#222222") -> None:
    """Horizontal dimension arrow between x0 and x1 at height y1."""
    ax.annotate("", xy=(x0, y1), xytext=(x1, y1),
                arrowprops=dict(arrowstyle="<->", color=color, lw=1.0))
    ax.text((x0 + x1) / 2.0, y1 + 0.18, text, ha="center", fontsize=8.5,
            color=color)


def render_layout(path: Path) -> Path:
    """Top view of the 4x4 array with pitch / patch / board dimensions."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    half = am.BOARD_MM / 2.0
    fig, ax = plt.subplots(figsize=(8.6, 8.0), dpi=160)

    ax.add_patch(plt.Rectangle((-half, -half), am.BOARD_MM, am.BOARD_MM,
                               facecolor=_SUB, edgecolor=_SUB_EDGE,
                               linewidth=1.4, zorder=0))
    ax.text(0, half - 0.55, f"RT/duroid 5880 · h = {am.SUBSTRATE_H_MM} mm",
            ha="center", fontsize=8.5, color="#1b3a5c", zorder=6)

    for j in range(am.N_Y):
        for i in range(am.N_X):
            x0, y0, x1, y1 = am.patch_rect(i, j)
            ax.add_patch(plt.Rectangle((x0, y0), x1 - x0, y1 - y0,
                                       facecolor=_CU, edgecolor=_CU_EDGE,
                                       linewidth=1.0, zorder=3))
            n = am.port_number(i, j)
            px, py = am.port_xy(n)
            ax.plot([px], [py], marker="o", ms=4.6, mfc="#d62728",
                    mec="#7f1113", zorder=5)
            ax.text(px + 0.16, py - 0.30, f"P{n}", fontsize=6.6,
                    color="#7f1113", zorder=5)

    # dimensions
    _dims(ax, -am.PITCH_MM / 2, 0, am.PITCH_MM / 2, -half - 1.15,
          f"d = {am.PITCH_MM:.3f} mm = $\\lambda_0$/2")
    _dims(ax, am.element_xy(1, 1)[0] - am.PATCH_L_MM / 2, 0,
          am.element_xy(1, 1)[0] + am.PATCH_L_MM / 2, half + 0.75,
          f"L = {am.PATCH_L_MM:.3f} mm")
    ax.annotate("", xy=(half + 0.75, -am.PATCH_W_MM / 2),
                xytext=(half + 0.75, am.PATCH_W_MM / 2),
                arrowprops=dict(arrowstyle="<->", color="#222222", lw=1.0))
    ax.text(half + 1.0, 0, f"W = {am.PATCH_W_MM:.3f} mm", rotation=90,
            va="center", fontsize=8.5, color="#222222")
    _dims(ax, -half, -half - 1.9, half, -half - 1.9,
          f"board {am.BOARD_MM:g} x {am.BOARD_MM:g} mm")

    # probe callout
    px, py = am.port_xy(6)
    ax.text(-half - 0.1, -half - 2.6,
            "red dots: 50 $\\Omega$ probe feeds (lumped ports, ground to patch)",
            fontsize=8.4, color="#7f1113", ha="left")

    ax.set_aspect("equal")
    ax.set_xlim(-half - 2.6, half + 2.6)
    ax.set_ylim(-half - 3.0, half + 1.6)
    ax.set_xlabel("x (mm)")
    ax.set_ylabel("y (mm)")
    ax.set_title("4x4 patch array, 77 GHz - top view (E-plane = xz)",
                 fontsize=11.5)
    ax.grid(True, which="both", color="#dddddd", linewidth=0.5, zorder=1)
    fig.tight_layout()
    fig.savefig(path)
    plt.close(fig)
    return path


def render_element(path: Path) -> Path:
    """Zoom on one element: patch outline, probe inset, polarisation."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    x0, y0, x1, y1 = am.patch_rect(1, 1)
    yc = 0.5 * (y0 + y1)
    px, py = am.port_xy(am.port_number(1, 1))
    fig, (ax, axz) = plt.subplots(1, 2, figsize=(10.6, 5.2), dpi=160,
                                  gridspec_kw={"width_ratios": [3, 2]})

    ax.add_patch(plt.Rectangle((x0 - 0.7, y0 - 0.7),
                               (x1 - x0) + 1.4, (y1 - y0) + 1.4,
                               facecolor=_SUB, edgecolor=_SUB_EDGE,
                               linewidth=1.1, zorder=0))
    ax.add_patch(plt.Rectangle((x0, y0), x1 - x0, y1 - y0, facecolor=_CU,
                               edgecolor=_CU_EDGE, linewidth=1.4, zorder=3))
    ax.plot([px], [py], marker="o", ms=7.0, mfc="#d62728", mec="#7f1113",
            zorder=5)

    _dims(ax, x0, 0, x1, y1 + 0.95,
          f"L = {am.PATCH_L_MM:.3f} mm (resonant, x)")
    ax.annotate("", xy=(x1 + 0.55, y0), xytext=(x1 + 0.55, y1),
                arrowprops=dict(arrowstyle="<->", color="#222222", lw=1.0))
    ax.text(x1 + 0.76, yc, f"W = {am.PATCH_W_MM:.3f} mm", rotation=90,
            va="center", fontsize=8.5)
    ax.annotate("", xy=(x0 - 0.1, y1 + 1.62), xytext=(x0 - 1.1, y1 + 1.62),
                arrowprops=dict(arrowstyle="-|>", color="#1b3a5c", lw=1.6))
    ax.text(x0 + 0.02, y1 + 1.60, "E $\\parallel$ x  ($\\varphi$ = 0 cut)",
            fontsize=8.4, color="#1b3a5c", va="center")

    # second panel: feed zoom (the 5.8 % inset is invisible on the whole patch)
    zi = axz
    zi.add_patch(plt.Rectangle((x0 - 0.6, y0 - 0.6), (x1 - x0) + 1.2,
                               (y1 - y0) + 1.2, facecolor=_SUB,
                               edgecolor=_SUB_EDGE, linewidth=0.8, zorder=0))
    zi.add_patch(plt.Rectangle((x0, y0), x1 - x0, y1 - y0, facecolor=_CU,
                               edgecolor=_CU_EDGE, linewidth=1.0, zorder=3))
    zi.plot([px], [py], marker="o", ms=8.0, mfc="#d62728", mec="#7f1113",
            zorder=5)
    zi.annotate("", xy=(x0, py + 0.13), xytext=(px, py + 0.13),
                arrowprops=dict(arrowstyle="<|-|>", color="#7f1113", lw=1.2,
                                mutation_scale=9))
    zi.text(0.5 * (x0 + px), py + 0.155,
            f"inset {am.FEED_INSET_MM:.3f} mm", color="#7f1113",
            fontsize=8.2, ha="center", va="bottom")
    zi.annotate("patch edge (x = -L/2)", xy=(x0, py - 0.10),
                xytext=(x0 + 0.02, py - 0.19), color="#315a86",
                fontsize=7.8, ha="left",
                arrowprops=dict(arrowstyle="->", color="#315a86", lw=0.8))
    zi.text(px + 0.02, py + 0.02, "50 $\\Omega$ probe", color="#7f1113",
            fontsize=8.0, ha="left", va="bottom")
    zi.set_xlim(x0 - 0.10, x0 + 0.42)
    zi.set_ylim(py - 0.26, py + 0.26)
    zi.set_aspect("equal")
    zi.set_xticks([])
    zi.set_yticks([])
    zi.set_title("feed zoom", fontsize=10)
    zi.grid(True, color="#dddddd", linewidth=0.5)

    ax.set_aspect("equal")
    ax.set_xlim(x0 - 1.3, x1 + 2.5)
    ax.set_ylim(y0 - 1.9, y1 + 2.1)
    ax.set_xlabel("x (mm)")
    ax.set_ylabel("y (mm)")
    ax.set_title("Single element: probe-fed rectangular patch", fontsize=11.5)
    ax.grid(True, which="both", color="#dddddd", linewidth=0.5, zorder=1)
    fig.tight_layout()
    fig.savefig(path)
    plt.close(fig)
    return path


def render_stackup(path: Path) -> Path:
    """Cross-section through a patch row (y = -0.973 mm), z exaggerated."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.patches import Rectangle

    h = am.SUBSTRATE_H_MM
    # cut through the patch row j = 1 (y = -0.973 mm)
    row_y = am.element_xy(0, 1)[1]
    half = am.BOARD_MM / 2.0
    z_sub, z_top = 0.0, h
    # z is drawn with an exaggeration factor so a 0.254 mm slab is visible
    zex = 6.0

    fig, ax = plt.subplots(figsize=(9.4, 4.6), dpi=160)

    def Z(z: float) -> float:
        return z * zex

    # substrate + ground as one band (x from -half..half)
    ax.add_patch(Rectangle((-half, Z(0)), am.BOARD_MM, Z(h),
                           facecolor=_SUB, edgecolor=_SUB_EDGE, lw=1.3,
                           zorder=2))
    ax.add_patch(Rectangle((-half, Z(0) - 0.16), am.BOARD_MM, 0.16,
                           facecolor=_CU, edgecolor=_CU_EDGE, lw=1.0,
                           zorder=3))

    # patches of the row
    for i in range(am.N_X):
        x0, _, x1, _ = am.patch_rect(i, 1)
        ax.add_patch(Rectangle((x0, Z(h)), x1 - x0, 0.14, facecolor=_CU,
                               edgecolor=_CU_EDGE, lw=1.0, zorder=4))
        px, _ = am.port_xy(am.port_number(i, 1))
        ax.plot([px, px], [Z(0), Z(h)], color="#d62728", lw=2.4, zorder=5)
        ax.plot([px], [Z(h) / 2], marker="o", ms=4.0, mfc="#d62728",
                mec="#7f1113", zorder=6)

    # annotations
    ax.annotate("ground plane (Cu 35 $\\mu$m, PEC)",
                xy=(-4.0, Z(0) - 0.08), xytext=(-5.9, -1.25),
                fontsize=8.8,
                arrowprops=dict(arrowstyle="->", lw=0.9, color="#333333"))
    ax.annotate(f"RT/duroid 5880\n$\\epsilon_r$ = {am.SUBSTRATE_EPS_R}, "
                f"tan $\\delta$ = {am.SUBSTRATE_TAND}\nh = {h} mm",
                xy=(-4.35, Z(h / 2)), xytext=(-8.55, 0.35), fontsize=8.8,
                ha="left", va="center",
                arrowprops=dict(arrowstyle="->", lw=0.9, color="#333333"))
    ax.annotate("patch (PEC sheet)",
                xy=(am.element_xy(2, 1)[0], Z(h) + 0.12),
                xytext=(0.1, 2.55), fontsize=8.8, ha="center",
                arrowprops=dict(arrowstyle="->", lw=0.9, color="#333333"))
    ax.annotate("50 $\\Omega$ probe\n(lumped port)",
                xy=(am.port_xy(1)[0], Z(h / 2)),
                xytext=(-3.2, 2.95), fontsize=8.8, color="#7f1113",
                ha="center",
                arrowprops=dict(arrowstyle="->", lw=0.9, color="#7f1113"))

    # z axis ticks
    for z, lab in ((0.0, "z = 0"), (h, f"z = h = {h} mm")):
        ax.axhline(Z(z), color="#999999", lw=0.5, ls=":")
    ax.text(half + 0.35, Z(0), "z = 0", va="center", fontsize=8.2,
            color="#555555")
    ax.text(half + 0.35, Z(h), f"z = h\n({h} mm)", va="center", fontsize=8.2,
            color="#555555")
    ax.text(half + 0.35, Z(h) + 0.55, "patches", va="center", fontsize=8.2,
            color="#555555")

    ax.set_aspect("equal")
    ax.set_xlim(-half - 2.9, half + 2.6)
    ax.set_ylim(-1.8, 3.4)
    ax.set_xlabel("x (mm)")
    ax.set_yticks([])
    ax.set_title("Stack-up along a patch row (y = "
                 f"{row_y:.3f} mm) - z exaggerated 6x", fontsize=11.5)
    ax.grid(True, axis="x", which="both", color="#dddddd", linewidth=0.5)
    fig.tight_layout()
    fig.savefig(path)
    plt.close(fig)
    return path


def render_array3d(path: Path) -> Path:
    """Isometric 3D view of the complete 4x4 array."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import numpy as np
    from mpl_toolkits.mplot3d.art3d import Poly3DCollection

    half = am.BOARD_MM / 2.0
    h = am.SUBSTRATE_H_MM
    # the substrate is 0.254 mm thick - exaggerate z for visibility
    zex = 8.0

    fig = plt.figure(figsize=(9.6, 7.4), dpi=160)
    ax = fig.add_subplot(111, projection="3d")
    # matplotlib's painter algorithm mis-sorts large coplanar boxes; draw in
    # explicit order (ground -> substrate -> patches) instead
    try:
        ax.computed_zorder = False
    except AttributeError:
        pass

    def box(x0, y0, x1, y1, z0, z1, color, alpha=1.0, lw=0.0, zorder=1,
            edge="none"):
        z0, z1 = z0 * zex, z1 * zex
        verts = [
            [(x0, y0, z0), (x1, y0, z0), (x1, y1, z0), (x0, y1, z0)],
            [(x0, y0, z1), (x1, y0, z1), (x1, y1, z1), (x0, y1, z1)],
            [(x0, y0, z0), (x1, y0, z0), (x1, y0, z1), (x0, y0, z1)],
            [(x0, y1, z0), (x1, y1, z0), (x1, y1, z1), (x0, y1, z1)],
            [(x0, y0, z0), (x0, y1, z0), (x0, y1, z1), (x0, y0, z1)],
            [(x1, y0, z0), (x1, y1, z0), (x1, y1, z1), (x1, y0, z1)],
        ]
        pc = Poly3DCollection(verts, facecolor=color, edgecolor=edge,
                              linewidth=lw, alpha=alpha, zorder=zorder)
        ax.add_collection3d(pc)

    # distinct colours: dark bronze ground, pale blue dielectric,
    # bright gold patches so the 4x4 pattern is unmistakable
    gnd_col = "#7a4a12"
    cu_col = "#f0b323"
    cu_edge = "#5a3a08"

    # ground plane first (bottom), then the translucent substrate on top of it
    box(-half, -half, half, half, -0.06, 0.0, gnd_col, alpha=1.0)
    box(-half, -half, half, half, 0.0, h, "#b9d2ea", alpha=0.80, lw=0.9,
        edge="#315a86")

    # patches, deliberately lifted by 0.03 mm (drawing only) and outlined
    for j in range(am.N_Y):
        for i in range(am.N_X):
            x0, y0, x1, y1 = am.patch_rect(i, j)
            box(x0, y0, x1, y1, h, h + 0.03, cu_col, alpha=1.0,
                lw=0.9, edge=cu_edge, zorder=6)
            px, py = am.port_xy(am.port_number(i, j))
            ax.plot([px, px], [py, py], [0, h * zex], color="#b3121b",
                    lw=2.4, zorder=8)
            ax.plot([px], [py], [h * zex], marker="o", ms=2.6,
                    mfc="#ff2d2d", mec="#7f1113", zorder=9)

    ax.set_xlim(-half - 1, half + 1)
    ax.set_ylim(-half - 1, half + 1)
    ax.set_zlim(-1.2, (h + 0.03) * zex + 2.6)
    try:
        ax.set_box_aspect((1, 1, 0.42))
    except Exception:
        pass
    ax.view_init(elev=34, azim=-58)
    ax.set_xlabel("x (mm)", labelpad=8)
    ax.set_ylabel("y (mm)", labelpad=8)
    ax.set_zlabel("z (mm)", labelpad=6)
    ax.set_title("4x4 patch array on RT/duroid 5880 - 3D view\n"
                 "(z exaggerated 8x; gold = patches, red = probe ports, "
                 "bronze = ground)", fontsize=11.0)
    fig.tight_layout()
    fig.savefig(path)
    plt.close(fig)
    return path


PREVIEWS = {
    "layout": (render_layout, "fig01_layout.png"),
    "element": (render_element, "fig02_element.png"),
    "stackup": (render_stackup, "fig03_stackup.png"),
    "array3d": (render_array3d, "fig04_array3d.png"),
}


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--only", choices=tuple(PROJECTS), default=None,
                    help="write a single project file")
    ap.add_argument("--preview", action="store_true",
                    help="render docs/figures/fig01..04 construction PNGs")
    args = ap.parse_args()

    written: list[Path] = []
    variants = [args.only] if args.only else list(PROJECTS)
    for v in variants:
        written.append(write_project(v))

    if args.preview:
        FIG.mkdir(parents=True, exist_ok=True)
        for render, name in PREVIEWS.values():
            written.append(render(FIG / name))

    print(am.summary())
    print()
    for p in written:
        print(f"wrote {p}  ({p.stat().st_size} bytes)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
