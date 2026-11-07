"""Diagnostic: spectral content of the port signals (usable DFT band).

The Gaussian excitation is band limited (fc = 0.5*(fstop-fstart)); outside
that band Zin(f) from the raw port signals is noise.  This prints |U(f)| and
|I(f)| so we know which part of an impedance curve is trustworthy.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


def load_signal(path: Path) -> tuple[np.ndarray, np.ndarray]:
    """Read an openEMS port signal file (columns: t u / t i, '%' header)."""
    t, v = [], []
    with open(path, "r", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if line.startswith("%"):
                continue
            parts = line.split()
            if len(parts) >= 2:
                t.append(float(parts[0]))
                v.append(float(parts[1]))
    return np.asarray(t), np.asarray(v)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("work", help="work directory (contains simulation/)")
    ap.add_argument("--port", type=int, default=1)
    ap.add_argument("--fmin", type=float, default=40e9)
    ap.add_argument("--fmax", type=float, default=120e9)
    args = ap.parse_args()

    sim = Path(args.work) / "simulation"
    tu, ut = load_signal(sim / f"port_ut_{args.port}")
    ti, it = load_signal(sim / f"port_it_{args.port}")

    print(f"ut: {len(ut)} samples, dt = {np.diff(tu).mean()*1e12:.3f} ps, "
          f"span {tu[-1]*1e12:.1f} ps")
    print(f"it: {len(it)} samples, dt = {np.diff(ti).mean()*1e12:.3f} ps, "
          f"span {ti[-1]*1e12:.1f} ps")
    print(f"ut envelope: max |u| = {np.abs(ut).max():.4g}, "
          f"last 10% max = {np.abs(ut[-len(ut)//10:]).max():.4g}")
    print(f"it envelope: max |i| = {np.abs(it).max():.4g}, "
          f"last 10% max = {np.abs(it[-len(it)//10:]).max():.4g}")

    freqs = np.linspace(args.fmin, args.fmax, 1601)
    u = np.array([np.trapz(ut * np.exp(-1j * 2 * np.pi * f * tu), tu)
                  for f in freqs])
    i = np.array([np.trapz(it * np.exp(-1j * 2 * np.pi * f * ti), ti)
                  for f in freqs])
    um, im = np.abs(u), np.abs(i)
    # relative to the in-band peak (74..80 GHz)
    band = (freqs >= 74e9) & (freqs <= 80e9)
    uref, iref = um[band].max(), im[band].max()
    print("\n   f [GHz]   |U| [dB rel band peak]   |I| [dB rel band peak]")
    for f in np.arange(args.fmin / 1e9, args.fmax / 1e9 + 1e-9, 2.0):
        k = int(np.argmin(np.abs(freqs - f * 1e9)))
        print(f"   {freqs[k]/1e9:7.1f}   {20*np.log10(um[k]/uref):8.2f}"
              f"            {20*np.log10(im[k]/iref):8.2f}")

    fig, ax = plt.subplots(1, 2, figsize=(11.4, 4.4), dpi=140)
    ax[0].plot(freqs / 1e9, 20 * np.log10(um / uref), label="|U|")
    ax[0].plot(freqs / 1e9, 20 * np.log10(im / iref), label="|I|")
    ax[0].axhline(-20, color="k", lw=0.7, ls="--")
    ax[0].axhline(-40, color="k", lw=0.7, ls=":")
    ax[0].set_xlabel("f [GHz]")
    ax[0].set_ylabel("spectrum [dB rel. in-band peak]")
    ax[0].set_title("port signal spectra (usable DFT band)")
    ax[0].legend()
    ax[0].grid(alpha=0.3)

    z = u / i
    ax[1].plot(freqs / 1e9, z.real, label="Re Z")
    ax[1].plot(freqs / 1e9, z.imag, label="Im Z")
    ax[1].axhline(0, color="k", lw=0.7)
    ax[1].axhline(50, color="gray", lw=0.7, ls="--")
    ax[1].set_xlabel("f [GHz]")
    ax[1].set_ylabel("Zin [ohm]")
    ax[1].set_title("Zin (greyed outside -20 dB excitation band)")
    ax[1].legend()
    ax[1].grid(alpha=0.3)

    # mask the region where either spectrum is below -20 dB
    mask = (20 * np.log10(um / uref) < -20) | (20 * np.log10(im / iref) < -20)
    for a in (ax[1],):
        a.fill_between(freqs / 1e9, -400, 400, where=mask, color="0.9",
                       zorder=0)

    out = Path(args.work) / "spectrum_check.png"
    fig.tight_layout()
    fig.savefig(out)
    print(f"\nwrote {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
