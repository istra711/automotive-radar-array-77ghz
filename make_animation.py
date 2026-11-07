"""Build the 3D beam-steering animation from the full-wave sweep runs.

Reads the archived 77 GHz gain patterns in ``docs/sweep/scan_*.npz`` (written
by ``run_sweep.py``) and renders one 3D lobe per steering angle with a fixed
camera, a fixed colour scale and a fixed radial scale, then assembles the
frames into ``docs/figures/anim_beam_steer.gif``.

Usage::

    python make_animation.py
    python make_animation.py --duration 300 --size 8.2
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from PIL import Image

HERE = Path(__file__).resolve().parent
SWEEP = HERE / "docs" / "sweep"
FIG = HERE / "docs" / "figures"
DB_FLOOR = 30.0          # colour/radius dynamic range in dB


def load_frames() -> list[dict]:
    frames = []
    for p in sorted(SWEEP.glob("scan_*.npz")):
        d = np.load(p)
        frames.append({
            "file": p,
            "steer": float(d["steer_deg"]),
            "theta": d["theta"].astype(float),
            "phi": d["phi"].astype(float),
            "gain": d["gain77"].astype(float),
            "peak_theta": float(d["peak_theta_deg"]),
            "peak_phi": float(d["peak_phi_deg"]),
            "peak_gain": float(d["peak_gain_dbi"]),
            "dmax": float(d["dmax77_dbi"]),
        })
    frames.sort(key=lambda f: f["steer"])
    return frames


def render(frame: dict, gmax: float, vmin: float, path: Path,
           size: float) -> None:
    th = np.radians(frame["theta"])
    ph = np.radians(frame["phi"])
    TH, PH = np.meshgrid(th, ph, indexing="ij")
    G = np.clip(frame["gain"], vmin, gmax)

    amp = 10.0 ** ((G - gmax) / 20.0)
    X = amp * np.sin(TH) * np.cos(PH)
    Y = amp * np.sin(TH) * np.sin(PH)
    Z = amp * np.cos(TH)

    norm = Normalize(vmin, gmax)
    cmap = plt.get_cmap("jet")

    fig = plt.figure(figsize=(size, size * 0.92), dpi=110)
    ax = fig.add_subplot(111, projection="3d")
    ax.plot_surface(X, Y, Z, facecolors=cmap(norm(G)), rstride=2, cstride=2,
                    linewidth=0, antialiased=False, shade=False)

    ax.set_box_aspect((1, 1, 0.9))
    lim = 1.06
    ax.set_xlim(-lim, lim)
    ax.set_ylim(-lim, lim)
    ax.set_zlim(-lim, lim)
    ax.set_xlabel("x", labelpad=-4)
    ax.set_ylabel("y", labelpad=-4)
    ax.set_zlabel("z", labelpad=-6)
    ax.tick_params(labelsize=7, pad=-2)
    ax.view_init(elev=22, azim=-132)
    ax.set_title(
        f"Beam steering, full wave, 77 GHz\n"
        f"$\\theta_0$ = {frame['steer']:+.0f}$^\\circ$   "
        f"peak $\\theta$ = {frame['peak_theta']:.0f}$^\\circ$, "
        f"$\\varphi$ = {frame['peak_phi']:.0f}$^\\circ$   "
        f"({frame['peak_gain']:.2f} dBi)",
        fontsize=10)

    m = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    m.set_array([])
    fig.colorbar(m, ax=ax, shrink=0.62, pad=0.02, label="gain (dBi)")
    fig.tight_layout(pad=0.4)
    fig.savefig(path)
    plt.close(fig)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--duration", type=int, default=260,
                    help="frame duration in ms")
    ap.add_argument("--size", type=float, default=7.6,
                    help="figure width in inches")
    ap.add_argument("--colors", type=int, default=256)
    args = ap.parse_args()

    frames = load_frames()
    if len(frames) < 2:
        raise SystemExit(f"need at least two sweep runs, found {len(frames)} "
                         f"in {SWEEP}")

    gmax = max(f["gain"].max() for f in frames)
    vmin = gmax - DB_FLOOR
    print(f"{len(frames)} frames, colour scale {vmin:.1f} .. {gmax:.1f} dBi, "
          f"steering {frames[0]['steer']:+.0f} .. {frames[-1]['steer']:+.0f} deg")

    FIG.mkdir(parents=True, exist_ok=True)
    tmp = HERE / "sweep_frames"
    tmp.mkdir(exist_ok=True)
    paths = []
    for i, fr in enumerate(frames):
        p = tmp / f"frame_{i:02d}.png"
        render(fr, gmax, vmin, p, args.size)
        paths.append(p)
        print(f"  frame {i+1:2d}/{len(frames)}  theta0={fr['steer']:+.0f}  "
              f"peak {fr['peak_theta']:.0f} deg @ {fr['peak_gain']:.2f} dBi")

    # static fallback frame: the +30 deg lobe is the most illustrative one
    pick = min(range(len(frames)), key=lambda i: abs(frames[i]["steer"] - 30.0))
    fallback = FIG / "anim_beam_steer_frame.png"
    Image.open(paths[pick]).save(fallback)
    print(f"wrote {fallback}")

    imgs = [Image.open(p).convert("RGB") for p in paths]
    quants = [im.quantize(colors=args.colors, method=Image.MEDIANCUT)
              for im in imgs]
    gif = FIG / "anim_beam_steer.gif"
    quants[0].save(gif, save_all=True, append_images=quants[1:],
                   duration=args.duration, loop=0, optimize=True,
                   disposal=2)
    size_mb = gif.stat().st_size / 1e6
    print(f"wrote {gif}  ({size_mb:.2f} MB, {len(quants)} frames, "
          f"{args.duration} ms)")

    # static companion: signed E-plane cuts (phi=0 for +x, phi=180 for -x)
    fig, ax = plt.subplots(figsize=(8.6, 4.6), dpi=150)
    cmap = plt.get_cmap("coolwarm")
    snorm = Normalize(min(f["steer"] for f in frames),
                      max(f["steer"] for f in frames))
    for fr in frames:
        if fr["steer"] < 0:
            j = int(np.argmin(np.abs(fr["phi"] - 180.0)))
            x = -fr["theta"]
            xpk = -fr["peak_theta"]
        else:
            j = int(np.argmin(np.abs(fr["phi"] - 0.0)))
            x = fr["theta"]
            xpk = fr["peak_theta"]
        col = cmap(snorm(fr["steer"]))
        ax.plot(x, fr["gain"][:, j], color=col, lw=1.6)
        ax.plot([xpk], [fr["peak_gain"]], "o", ms=3.5, color=col)
    ax.axhline(gmax - 3.0, color="#7f8c8d", lw=0.8, ls=":")
    ax.set_xlabel("signed theta (deg)  -  E plane "
                  "($\\varphi$ = 0 for $\\theta_0>0$, $\\varphi$ = 180 for "
                  "$\\theta_0<0$)")
    ax.set_ylabel("gain (dBi), 77 GHz")
    ax.set_title("Full-wave beam steering sweep - E-plane cuts")
    ax.set_xlim(-90, 90)
    ax.set_ylim(gmax - 30, gmax + 3)
    ax.grid(True, alpha=0.4)
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=snorm)
    sm.set_array([])
    fig.colorbar(sm, ax=ax, label="steering angle $\\theta_0$ (deg)")
    fig.tight_layout()
    sweep_fig = FIG / "fig15_beam_sweep.png"
    fig.savefig(sweep_fig)
    plt.close(fig)
    print(f"wrote {sweep_fig}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
