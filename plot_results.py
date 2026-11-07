"""Result figures + docs/results.json for the 4x4 patch array study.

Reads the finished runs and writes ``docs/figures/fig05..fig12.png`` (the
construction previews are fig01..fig04 from ``make_project.py``) and
``docs/results.json``, the number source for ``build_docs.py``.

Runs consumed
-------------
============================  =============================================
``work_element``              final element, production band 74..80 GHz
``work_patch_element_77ghz_wideband``  the *before* design (Hammerstad
                              L = 1.1383 mm, inset 0.30 L), 47..107 GHz
``work_probe_only_77ghz``     open probe control -> feed line model
``work_inset005_77ghz``       inset control (R_edge)
``work_cand_L1p1610_s0p050``  retune verification run (identical geometry)
``work_array_embedded``       P6 excited -> coupling column S*n6 +
                              embedded element pattern
``work_array_broadside``      16 ports in phase
``work_array_scan30``         progressive phase, beam at +30 deg
============================  =============================================

Post-processing notes
---------------------
*  **Efficiency / gain of the driven arrays.**  The bundled
  ``postproc.nf2ff._load_pin_spectral`` takes the accepted power of the
  *largest* port only - correct for the single-driven runs (element,
  embedded), but 16x too small for broadside/scan30, which inflates
  eta/gain by ~12 dB.  This script sums the accepted power of **all driven
  ports** (same pulse-DFT convention ``2*sum(x*exp(-j2*pi*f*t))*dt`` as the
  NF2FF Prad, so the ratio stays unitless) and reports corrected values.
  Dmax and the patterns themselves are unaffected (pure field integrals).
*  **Active reflection of the phased runs.**  One FDTD run measures the
  waves referenced to the *first* excited port: the stored column is
  ``S_k1 = uf_ref_k/uf_inc_1``.  With excitation phases phi_k the active
  reflection of port k is ``Gamma_k = S_k1 * exp(-j*(phi_k - phi_1))``.

Usage::

    python plot_results.py
"""

from __future__ import annotations

import json
import math
from pathlib import Path

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize

import array_model as am
from analyze_zin import load_signal

HERE = Path(__file__).resolve().parent
FIG = HERE / "docs" / "figures"
F77 = 77.0e9
DB_FLOOR = -30.0                 # pattern dynamic range in the figures

RUNS = {
    "element": "work_element",
    "before": "work_patch_element_77ghz_wideband",
    "probe": "work_probe_only_77ghz",
    "inset005": "work_inset005_77ghz",
    "candA": "work_cand_L1p1610_s0p050",
    "embedded": "work_array_embedded",
    "broadside": "work_array_broadside",
    "scan30": "work_array_scan30",
}


# ---------------------------------------------------------------------------
# loaders
# ---------------------------------------------------------------------------

def load_s(work: str, key: str = "S11") -> tuple[np.ndarray, np.ndarray]:
    d = json.loads((HERE / work / "sparams" / "sparams.json")
                   .read_text(encoding="utf-8"))
    f = np.array(d["freq"])
    s = np.array(d[key]["real"]) + 1j * np.array(d[key]["imag"])
    return f, s


def load_nf(work: str) -> dict:
    d = json.loads((HERE / work / "nf2ff" / "nf2ff.json")
                   .read_text(encoding="utf-8"))
    out = {
        "freq": np.array(d["frequencies"]),
        "theta": np.array(d["theta_deg"]),
        "phi": np.array(d["phi_deg"]),
        "D": np.array(d["D_dBi"]),
        "gain": np.array(d["gain_dBi"]),
        "Dmax": np.array(d["Dmax_dBi"]),
        "eff": np.array(d["efficiency"]),
        "Prad": np.array(d["Prad"]),
    }
    return out


def s11_db(f: np.ndarray, s: np.ndarray) -> np.ndarray:
    return 20.0 * np.log10(np.clip(np.abs(s), 1e-12, None))


def z_of(f: np.ndarray, s: np.ndarray) -> np.ndarray:
    return 50.0 * (1.0 + s) / (1.0 - s)


def f77_index(nf: dict) -> int:
    return int(np.argmin(np.abs(nf["freq"] - F77)))


def s77_index(f: np.ndarray) -> int:
    return int(np.argmin(np.abs(f - F77)))


def bw_10db(f: np.ndarray, s: np.ndarray) -> tuple[float, float, float]:
    """(-10 dB band around 77 GHz): (f_low, f_high, bandwidth) in Hz."""
    db = s11_db(f, s)
    k0 = s77_index(f)
    if db[k0] > -10.0:
        return float("nan"), float("nan"), float("nan")
    lo = k0
    while lo > 0 and db[lo - 1] <= -10.0:
        lo -= 1
    hi = k0
    while hi < len(f) - 1 and db[hi + 1] <= -10.0:
        hi += 1
    return float(f[lo]), float(f[hi]), float(f[hi] - f[lo])


def dft_pulse(x: np.ndarray, t: np.ndarray, f: np.ndarray) -> np.ndarray:
    """openEMS pulse-DFT convention (matches the NF2FF Prad units)."""
    dt = float(np.mean(np.diff(t)))
    ph = np.exp(-1j * 2.0 * np.pi * f[None, :] * t[:, None])
    return 2.0 * (x[:, None] * ph).sum(axis=0) * dt


def accepted_power(work: str, freqs: np.ndarray) -> np.ndarray:
    """P_acc(f) = sum over ALL driven ports of 0.5*Re(U*conj(I)).

    The DFT of **both** signals uses the time axis of ``port_ut``.  openEMS
    samples the voltage and the current on slightly different FDTD sub-times
    (the ``port_it`` time column is ~0.03..0.04 dt ahead), so transforming
    each signal on its own axis adds a spurious relative phase.  The reactive
    patch current is large, which biased Pin by several percent (11 % on the
    wideband *before* run).  ``postproc.nf2ff._load_pin_spectral`` uses the
    ``port_ut`` axis for both, hence the same convention is applied here.
    """
    sim = HERE / work / "simulation"
    total = np.zeros(len(freqs), dtype=float)
    n_ports = 0
    for n in range(1, 17):
        ut_p = sim / f"port_ut_{n}"
        it_p = sim / f"port_it_{n}"
        if not (ut_p.is_file() and it_p.is_file()):
            continue
        tu, ut = load_signal(ut_p)
        _, it = load_signal(it_p)
        u = dft_pulse(ut, tu, freqs)
        i = dft_pulse(it, tu, freqs)
        total += 0.5 * np.real(u * np.conj(i))
        n_ports += 1
    if n_ports == 0:
        raise RuntimeError(f"no port signals in {work}")
    return total


def array_factor(theta_deg: np.ndarray, phi_deg: np.ndarray,
                 pos: np.ndarray, phases_deg: np.ndarray) -> np.ndarray:
    """Normalised array-factor power |AF|^2 (peak = 1 at the design beam)."""
    th = np.radians(theta_deg)
    ph = np.radians(phi_deg)
    proj = (pos[:, 0][:, None] * (np.sin(th) * np.cos(ph))[None, :]
            + pos[:, 1][:, None] * (np.sin(th) * np.sin(ph))[None, :])
    total = 2.0 * np.pi / am.LAMBDA0_MM * proj + np.radians(phases_deg)[:, None]
    af = np.exp(1j * total).sum(axis=0) / len(pos)
    return np.abs(af) ** 2


def element_positions() -> tuple[np.ndarray, np.ndarray]:
    pos, ph = [], []
    for j in range(am.N_Y):
        for i in range(am.N_X):
            x, y = am.element_xy(i, j)
            pos.append((x, y))
            ph.append(am.scan_phases_deg(am.SCAN30_THETA_DEG)[am.port_number(i, j)])
    return np.array(pos), np.array(ph)


# ---------------------------------------------------------------------------
# figures
# ---------------------------------------------------------------------------

def fig_element_zin(before: bool = True) -> dict:
    """fig05: Zin(f) + |S11|(f) of the final element vs the before design."""
    f_a, s_a = load_s(RUNS["element"])
    z_a = z_of(f_a, s_a)
    k = s77_index(f_a)
    lo, hi, bw = bw_10db(f_a, s_a)

    fig, ax = plt.subplots(1, 2, figsize=(11.4, 4.3), dpi=140)
    ax[0].plot(f_a / 1e9, z_a.real, color="#1f5fa8", lw=1.8,
               label=r"Re$\,Z_\mathrm{in}$")
    ax[0].plot(f_a / 1e9, z_a.imag, color="#c0392b", lw=1.8, ls="--",
               label=r"Im$\,Z_\mathrm{in}$")
    ax[0].axhline(50, color="gray", lw=0.9, ls=":")
    ax[0].axhline(0, color="k", lw=0.7)
    ax[0].axvline(77, color="#7f8c8d", lw=0.9, ls=":")
    ax[0].plot([77], [z_a[k].real], "o", color="#1f5fa8", ms=5)
    ax[0].plot([77], [z_a[k].imag], "o", color="#c0392b", ms=5)
    ax[0].set_xlabel("Frequency (GHz)")
    ax[0].set_ylabel("Zin (Ohm)")
    ax[0].set_title(f"Input impedance  -  Zin(77) = {z_a[k].real:.1f} "
                    f"{z_a[k].imag:+.1f}j Ohm")
    ax[0].grid(True, alpha=0.4)
    ax[0].legend(loc="best", fontsize=9)

    if before:
        f_b, s_b = load_s(RUNS["before"])
        m = (f_b >= 73.5e9) & (f_b <= 80.5e9)
        ax[1].plot(f_b[m] / 1e9, s11_db(f_b[m], s_b[m]), color="#95a5a6",
                   lw=1.6, ls="--",
                   label="before: Hammerstad L = 1.1383 mm, inset 0.30 L")
    ax[1].plot(f_a / 1e9, s11_db(f_a, s_a), color="#1f5fa8", lw=2.0,
               label="final: L = 1.1610 mm, inset 0.05 L")
    ax[1].axhline(-10, color="#c0392b", lw=1.0, ls=":")
    ax[1].text(74.1, -9.2, "-10 dB", color="#c0392b", fontsize=8)
    ax[1].axvline(77, color="#7f8c8d", lw=0.9, ls=":")
    ax[1].set_xlabel("Frequency (GHz)")
    ax[1].set_ylabel("|S11| (dB)")
    ax[1].set_ylim(-28, 0)
    ax[1].set_title(f"Reflection  -  |S11|(77) = {s11_db(f_a, s_a)[k]:.1f} dB,"
                    f" -10 dB span {bw/1e9:.2f} GHz")
    ax[1].grid(True, alpha=0.4)
    ax[1].legend(loc="best", fontsize=8)

    fig.tight_layout()
    out = FIG / "fig05_element_zin_s11.png"
    fig.savefig(out)
    plt.close(fig)
    print(f"wrote {out}")
    return {
        "zin77_real": float(z_a[k].real), "zin77_imag": float(z_a[k].imag),
        "s11_77_db": float(s11_db(f_a, s_a)[k]),
        "best_s11_db": float(np.min(s11_db(f_a, s_a))),
        "best_s11_ghz": float(f_a[np.argmin(s11_db(f_a, s_a))] / 1e9),
        "bw10_low_ghz": lo / 1e9, "bw10_high_ghz": hi / 1e9,
        "bw10_ghz": bw / 1e9,
        "vswr77": float((1 + abs(s_a[k])) / (1 - abs(s_a[k]))),
    }


def fig_retune_journey() -> dict:
    """fig06: probe control + Zin loci of every design step (target 50 Ohm)."""
    from fit_model import (fit_probe, probe_theta, z_from_signals)

    # probe measurement + TL fit
    f_o, z_o = z_from_signals(RUNS["probe"])
    zp, eps_p, phi0 = fit_probe(f_o, z_o)
    th77 = float(probe_theta(np.array([F77]), eps_p, phi0)[0])

    designs = [
        ("before (inset 0.30 L, L=1.138)", RUNS["before"], "#95a5a6", "--"),
        ("inset 0.05 L, L = 1.1383", RUNS["inset005"], "#e67e22", "-."),
        ("final: L = 1.1610, inset 0.05 L", RUNS["element"], "#1f5fa8", "-"),
    ]

    fig, ax = plt.subplots(1, 2, figsize=(11.4, 4.3), dpi=140)
    band = (f_o >= 50e9) & (f_o <= 105e9)
    pred = -zp / np.tan(probe_theta(f_o[band], eps_p, phi0))
    ax[0].plot(f_o[band] / 1e9, z_o[band].imag, color="#c0392b", lw=1.8,
               label="Im Z measured (probe, open top)")
    ax[0].plot(f_o[band] / 1e9, pred, color="#1f5fa8", lw=1.4, ls="--",
               label=f"open line fit: Zp = {zp:.0f} Ohm, "
                     rf"$\epsilon_p$ = {eps_p:.2f}")
    ax[0].plot(f_o[band] / 1e9, z_o[band].real, color="#7f8c8d", lw=1.2,
               label="Re Z measured")
    ax[0].axhline(0, color="k", lw=0.7)
    ax[0].set_xlabel("Frequency (GHz)")
    ax[0].set_ylabel("Z of the bare feed probe (Ohm)")
    ax[0].set_title(f"Control run probe_only: capacitive "
                    f"({z_o[s77_index(f_o)].imag:.0f}j at 77 GHz - "
                    "a short line, not an inductor)")
    ax[0].grid(True, alpha=0.4)
    ax[0].legend(loc="best", fontsize=8)

    # Zin loci with constant-|S11| circles around 50 Ohm
    for label, work, color, ls in designs:
        f, s = load_s(work)
        z = z_of(f, s)
        m = (f >= 74e9) & (f <= 80e9) if work != RUNS["before"] else \
            (f >= 74e9) & (f <= 80e9)
        ax[1].plot(z[m].real, z[m].imag, color=color, lw=1.6, ls=ls,
                   label=label)
        k = s77_index(f)
        ax[1].plot(z[k].real, z[k].imag, "o", color=color, ms=5)
    # |S11| circles around 50 Ohm: |Z-50|/|Z+50| = 0.316 (-10 dB), 0.1 (-20 dB)
    for db, col in ((-10, "#c0392b"), (-20, "#27ae60")):
        gam = 10 ** (db / 20.0)
        # circle centre/radius in the Z plane for the reflection of a circle
        r_c = 50 * (1 + gam ** 2) / (1 - gam ** 2)
        r_r = 100 * gam / (1 - gam ** 2)
        th = np.linspace(0, 2 * np.pi, 200)
        ax[1].plot(r_c + r_r * np.cos(th), r_r * np.sin(th),
                   color=col, lw=0.9, ls=":")
        ax[1].text(r_c + r_r + 1, 0, f"{db} dB", color=col, fontsize=7,
                   va="center")
    ax[1].plot(50, 0, "*", color="#27ae60", ms=14, zorder=5)
    ax[1].text(52, -6, "50 Ohm", color="#27ae60", fontsize=8)
    ax[1].set_xlabel("Re Zin (Ohm)")
    ax[1].set_ylabel("Im Zin (Ohm)")
    ax[1].set_title("Zin loci, 74-80 GHz (dots = 77 GHz)")
    ax[1].grid(True, alpha=0.4)
    ax[1].legend(loc="upper right", fontsize=8)

    fig.tight_layout()
    out = FIG / "fig06_retune_journey.png"
    fig.savefig(out)
    plt.close(fig)
    print(f"wrote {out}")
    return {"probe_zp_ohm": zp, "probe_eps_p": eps_p,
            "probe_theta77_deg": math.degrees(th77),
            "probe_z77_imag": float(z_o[s77_index(f_o)].imag)}


def _plot_pattern_3d(ax, nf: dict, title: str, fi: int) -> float:
    th = np.radians(nf["theta"])
    ph = np.radians(nf["phi"])
    TH, PH = np.meshgrid(th, ph, indexing="ij")
    G = nf["gain"][fi]
    peak = float(np.max(G))
    Gc = np.clip(G, peak + DB_FLOOR, None)
    amp = 10 ** ((Gc - peak) / 20.0)
    X = amp * np.sin(TH) * np.cos(PH)
    Y = amp * np.sin(TH) * np.sin(PH)
    Z = amp * np.cos(TH)
    norm = Normalize(peak + DB_FLOOR, peak)
    cmap = plt.get_cmap("jet")
    ax.plot_surface(X, Y, Z, facecolors=cmap(norm(Gc)), rstride=1, cstride=1,
                    linewidth=0, antialiased=False, shade=False)
    ax.set_title(title, fontsize=10)
    ax.set_box_aspect((1, 1, 0.85))
    lim = 1.05
    ax.set_xlim(-lim, lim)
    ax.set_ylim(-lim, lim)
    ax.set_zlim(-lim, lim)
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_zlabel("z")
    ax.view_init(elev=24, azim=-132)
    return peak


def fig_element_pattern() -> dict:
    """fig07: 3D gain pattern of the isolated (final) element at 77 GHz."""
    nf = load_nf(RUNS["element"])
    fi = f77_index(nf)
    fig = plt.figure(figsize=(6.4, 5.4), dpi=140)
    ax = fig.add_subplot(111, projection="3d")
    peak = _plot_pattern_3d(
        ax, nf,
        f"Isolated element, gain at 77 GHz (peak {nf['gain'][fi].max():.2f} dBi)",
        fi)
    m = plt.cm.ScalarMappable(
        cmap="jet", norm=Normalize(nf["gain"][fi].max() + DB_FLOOR,
                                   nf["gain"][fi].max()))
    m.set_array([])
    fig.colorbar(m, ax=ax, shrink=0.6, pad=0.02, label="gain (dBi)")
    fig.tight_layout()
    out = FIG / "fig07_element_pattern3d.png"
    fig.savefig(out)
    plt.close(fig)
    print(f"wrote {out}")
    return {"dmax77": float(nf["Dmax"][fi]), "eff77": float(nf["eff"][fi]),
            "gain77": float(np.max(nf["gain"][fi])), "peak77": peak}


def _cut(nf: dict, fi: int, phi_deg: float) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    j = int(round(phi_deg / 180.0 * 90))     # phi grid: 2 deg steps
    return (nf["theta"], nf["D"][fi][:, j], nf["gain"][fi][:, j])


def fig_element_vs_embedded() -> dict:
    """fig08: isolated vs embedded element pattern (E/H plane cuts, 77 GHz)."""
    iso = load_nf(RUNS["element"])
    emb = load_nf(RUNS["embedded"])
    fi = f77_index(iso)

    fig, ax = plt.subplots(1, 2, figsize=(11.0, 4.3), dpi=140,
                           subplot_kw={"projection": "polar"})
    for a, (phi, name) in zip(ax, ((0.0, "E-plane (xz, phi = 0 deg)"),
                                   (90.0, "H-plane (yz, phi = 90 deg)"))):
        th, d_i, g_i = _cut(iso, fi, phi)
        _, d_e, g_e = _cut(emb, fi, phi)
        a.plot(np.radians(th), g_i, color="#1f5fa8", lw=1.8,
               label="isolated element")
        a.plot(np.radians(th), g_e, color="#c0392b", lw=1.8, ls="--",
               label="embedded (P6 in the 4x4 array)")
        a.set_theta_zero_location("N")
        a.set_theta_direction(-1)
        a.set_thetagrids(range(0, 360, 45))
        lo = min(g_i.min(), g_e.min())
        a.set_rlim(lo - 2, max(g_i.max(), g_e.max()) + 2)
        a.set_title(name, fontsize=10, pad=12)
        a.grid(True, alpha=0.4)
        a.legend(loc="lower center", bbox_to_anchor=(0.5, -0.22),
                 fontsize=8, ncol=2)
    fig.tight_layout()
    out = FIG / "fig08_element_vs_embedded.png"
    fig.savefig(out)
    plt.close(fig)
    print(f"wrote {out}")
    fi_e = f77_index(emb)
    return {"embedded_dmax77": float(emb["Dmax"][fi_e]),
            "embedded_eff77": float(emb["eff"][fi_e])}


def fig_coupling() -> dict:
    """fig09: coupling column S*n6 of the embedded run."""
    f, _ = load_s(RUNS["embedded"], "S_P6")
    s66 = load_s(RUNS["embedded"], "S_P6")[1]
    k = s77_index(f)

    entries = []
    for n in range(1, 17):
        _, s = load_s(RUNS["embedded"], f"S_P{n}")
        xk, yk = am.port_xy(n)
        x6, y6 = am.port_xy(6)
        dist = math.hypot(xk - x6, yk - y6)
        entries.append((n, dist, float(s11_db(np.array([f[k]]), np.array([s[k]]))[0])))
    entries.sort(key=lambda e: e[1])

    fig, ax = plt.subplots(1, 2, figsize=(11.0, 4.3), dpi=140)
    dists = [e[1] for e in entries if e[0] != 6]
    cpls = [e[2] for e in entries if e[0] != 6]
    ax[0].scatter(dists, cpls, s=42, color="#1f5fa8", zorder=3)
    for idx, (n, d, c) in enumerate(entries):
        if n != 6:
            dy = 4 if idx % 2 == 0 else -10
            ax[0].annotate(f"P{n}", (d, c), fontsize=7,
                           xytext=(4, dy), textcoords="offset points")
    ax[0].axhline(s11_db(np.array([f[k]]), np.array([s66[k]]))[0],
                  color="#c0392b", ls="--", lw=1.3,
                  label=f"|S66| = {s11_db(np.array([f[k]]), np.array([s66[k]]))[0]:.1f} dB")
    ax[0].set_xlabel("distance from P6 (mm)")
    ax[0].set_ylabel("|S_i6| (dB)")
    ax[0].set_title("Coupling column at 77 GHz (embedded run)")
    ax[0].grid(True, alpha=0.4)
    ax[0].legend(fontsize=8)

    for n, style in ((5, "-"), (2, "--"), (1, "-."), (10, ":")):
        _, s = load_s(RUNS["embedded"], f"S_P{n}")
        ax[1].plot(f / 1e9, s11_db(f, s), ls=style, lw=1.5,
                   label=f"|S{n}6| @ {math.dist(am.port_xy(n), am.port_xy(6)):.2f} mm")
    ax[1].plot(f / 1e9, s11_db(f, s66), color="#c0392b", lw=1.8,
               label="|S66| (self)")
    ax[1].set_xlabel("Frequency (GHz)")
    ax[1].set_ylabel("|S| (dB)")
    ax[1].set_title("Coupling vs frequency")
    ax[1].grid(True, alpha=0.4)
    ax[1].legend(fontsize=8)

    fig.tight_layout()
    out = FIG / "fig09_coupling.png"
    fig.savefig(out)
    plt.close(fig)
    print(f"wrote {out}")
    near = [c for n, d, c in entries if n != 6 and d < 2.0]
    return {"s66_77_db": float(s11_db(np.array([f[k]]), np.array([s66[k]]))[0]),
            "coupling_max_77_db": float(max(cpls)),
            "coupling_nearest_77_db": float(max(near)) if near else None,
            "coupling_table": [{"port": n, "dist_mm": d, "db": c}
                               for n, d, c in entries if n != 6]}


def _phased_active_gamma(work: str, phases_deg: dict[int, float],
                         f: np.ndarray) -> np.ndarray:
    """Active |Gamma_k| per port at every frequency (see module docstring)."""
    ph1 = phases_deg.get(1, 0.0)
    gam = []
    for n in range(1, 17):
        _, s = load_s(work, f"S_P{n}")
        dphi = math.radians(phases_deg.get(n, 0.0) - ph1)
        gam.append(s * np.exp(-1j * dphi))
    return np.array(gam)          # (16, nfreq)


def fig_array_pattern(work: str, name: str, figid: str,
                      phases_deg: dict[int, float] | None) -> dict:
    """3D pattern + E/H cuts vs the array-factor prediction."""
    nf = load_nf(work)
    fi = f77_index(nf)
    pos, ph30 = element_positions()
    if phases_deg is None:                       # broadside
        exc_ph = np.zeros(len(pos))
        beam_deg = 0.0
    else:
        exc_ph = np.array([phases_deg[n] for n in range(1, 17)])
        beam_deg = am.SCAN30_THETA_DEG

    iso = load_nf(RUNS["element"])
    th_e, _, g_e_iso = _cut(iso, fi, 0.0)
    th_h, _, g_h_iso = _cut(iso, fi, 90.0)
    af_e = array_factor(th_e, np.zeros_like(th_e), pos, exc_ph)
    af_h = array_factor(th_h, np.full_like(th_h, 90.0), pos, exc_ph)
    g_pred_E = g_e_iso + 10 * np.log10(np.clip(af_e, 1e-12, None))
    g_pred_H = g_h_iso + 10 * np.log10(np.clip(af_h, 1e-12, None))

    fig = plt.figure(figsize=(12.8, 4.8), dpi=140)
    ax3 = fig.add_subplot(1, 3, 1, projection="3d")
    peak = _plot_pattern_3d(
        ax3, nf,
        f"{name}\n77 GHz, peak {nf['gain'][fi].max():.2f} dBi", fi)
    m = plt.cm.ScalarMappable(
        cmap="jet", norm=Normalize(nf["gain"][fi].max() + DB_FLOOR,
                                   nf["gain"][fi].max()))
    m.set_array([])
    fig.colorbar(m, ax=ax3, shrink=0.55, pad=0.02, label="gain (dBi)")

    ax2 = fig.add_subplot(1, 3, 2)
    th, d_fw, g_fw = _cut(nf, fi, 0.0)
    floor = float(np.max(g_fw)) - 45.0     # 45 dB display window
    top = float(np.max(g_fw)) + 8.0
    ax2.plot(th, g_fw, color="#1f5fa8", lw=1.9, label="full wave")
    if g_pred_E is not None:
        ax2.plot(th, np.clip(g_pred_E, floor, None), color="#c0392b", lw=1.5,
                 ls="--", label="element x array factor")
        ax2.axvline(beam_deg, color="#27ae60", lw=1.1, ls=":")
        ax2.text(beam_deg + 2.0, top - 4.0,
                 f"design beam {beam_deg:g}$^\\circ$", color="#27ae60",
                 fontsize=8)
    ax2.set_ylim(floor, top)
    pk = int(np.argmax(g_fw))
    ax2.plot(th[pk], g_fw[pk], "o", color="#1f5fa8", ms=6)
    ax2.annotate(f"peak {th[pk]:g}$^\\circ$", (th[pk], g_fw[pk]),
                 xytext=(8, -14), textcoords="offset points", fontsize=8)
    ax2.set_xlabel("theta (deg), phi = 0 (E-plane)")
    ax2.set_ylabel("gain (dBi)")
    ax2.set_title("E-plane cut", fontsize=10)
    ax2.grid(True, alpha=0.4)
    ax2.legend(fontsize=8, loc="lower left")

    ax3b = fig.add_subplot(1, 3, 3)
    th, d_fw, g_fw_h = _cut(nf, fi, 90.0)
    ax3b.plot(th, g_fw_h, color="#1f5fa8", lw=1.9, label="full wave")
    if g_pred_H is not None:
        ax3b.plot(th, np.clip(g_pred_H, floor, None), color="#c0392b", lw=1.5,
                  ls="--", label="element x array factor")
    ax3b.set_ylim(floor, top)
    ax3b.set_xlabel("theta (deg), phi = 90 (H-plane)")
    ax3b.set_title("H-plane cut", fontsize=10)
    ax3b.grid(True, alpha=0.4)
    ax3b.legend(fontsize=8, loc="lower left")

    fig.tight_layout()
    out = FIG / f"{figid}"
    fig.savefig(out)
    plt.close(fig)
    print(f"wrote {out}")

    # beam direction from the 3D pattern
    G = nf["gain"][fi]
    t_idx, p_idx = np.unravel_index(int(np.argmax(G)), G.shape)
    theta_pk = float(nf["theta"][t_idx])
    phi_pk = float(nf["phi"][p_idx])
    return {"dmax77": float(nf["Dmax"][fi]),
            "gain77": float(nf["gain"][fi][t_idx, p_idx]),
            "eff77_raw": float(nf["eff"][fi]),
            "theta_peak_deg": theta_pk, "phi_peak_deg": phi_pk}


def fig_scan_summary(broad: dict, scan: dict) -> dict:
    """fig12: E-plane overlay of both excitations + active |Gamma| bars."""
    nf_b = load_nf(RUNS["broadside"])
    nf_s = load_nf(RUNS["scan30"])
    fi = f77_index(nf_b)
    ph_b = {n: 0.0 for n in range(1, 17)}
    ph_s = am.scan_phases_deg(am.SCAN30_THETA_DEG)

    f_b, _ = load_s(RUNS["broadside"])
    g_b = _phased_active_gamma(RUNS["broadside"], ph_b, f_b)
    g_s = _phased_active_gamma(RUNS["scan30"], ph_s, f_b)
    k = s77_index(f_b)

    fig, ax = plt.subplots(1, 2, figsize=(11.0, 4.3), dpi=140)
    th, _, gb = _cut(nf_b, fi, 0.0)
    _, _, gs = _cut(nf_s, fi, 0.0)
    ax[0].plot(th, gb, color="#1f5fa8", lw=1.9, label="broadside (in phase)")
    ax[0].plot(th, gs, color="#c0392b", lw=1.9, label="scan +30$^\\circ$")
    ax[0].axvline(30, color="#27ae60", lw=1.1, ls=":")
    ax[0].text(30.5, ax[0].get_ylim()[1] * 0.9 if ax[0].get_ylim()[1] else 0,
               "30$^\\circ$", color="#27ae60", fontsize=8)
    ax[0].set_xlabel("theta (deg), phi = 0 (E-plane)")
    ax[0].set_ylabel("gain (dBi)")
    ax[0].set_title("Beam steering at 77 GHz")
    ax[0].grid(True, alpha=0.4)
    ax[0].legend(fontsize=8)

    db_b = 20 * np.log10(np.clip(np.abs(g_b[:, k]), 1e-9, None))
    db_s = 20 * np.log10(np.clip(np.abs(g_s[:, k]), 1e-9, None))
    xs = np.arange(1, 17)
    w = 0.38
    ax[1].bar(xs - w / 2, db_b, w, color="#1f5fa8", label="broadside")
    ax[1].bar(xs + w / 2, db_s, w, color="#c0392b", label="scan +30$^\\circ$")
    ax[1].set_xticks(xs)
    ax[1].set_xlabel("port")
    ax[1].set_ylabel(r"active $|\Gamma_k|$ (dB)")
    ax[1].set_title("Active reflection per port at 77 GHz")
    ax[1].grid(True, axis="y", alpha=0.4)
    ax[1].legend(fontsize=8)

    fig.tight_layout()
    out = FIG / "fig12_scan_summary.png"
    fig.savefig(out)
    plt.close(fig)
    print(f"wrote {out}")
    return {
        "gamma_active_broadside_db": [float(v) for v in db_b],
        "gamma_active_scan30_db": [float(v) for v in db_s],
        "phase_step_deg": 360.0 * am.PITCH_MM
        * math.sin(math.radians(am.SCAN30_THETA_DEG)) / am.LAMBDA0_MM,
    }


# ---------------------------------------------------------------------------
# near-field snapshots
# ---------------------------------------------------------------------------

def _read_dump(path: Path) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Complex E field (3, nz, ny, nx) at 77 GHz plus the three axes."""
    import h5py

    with h5py.File(path, "r") as f:
        re = np.asarray(f["FieldData/FD/f0_real"], dtype=float)
        im = np.asarray(f["FieldData/FD/f0_imag"], dtype=float)
        x = np.asarray(f["Mesh/x"], dtype=float)
        y = np.asarray(f["Mesh/y"], dtype=float)
        z = np.asarray(f["Mesh/z"], dtype=float)
    return re + 1j * im, x, y, z


def fig_fields() -> dict:
    """fig13/14: near-field snapshots of the broadside array at 77 GHz."""
    dump = HERE / RUNS["broadside"] / "simulation"

    # --- E_xz cut (y = 0): |E| over the xz plane -------------------------
    e, x, y, z = _read_dump(dump / "E_xz.h5")
    amp = np.sqrt((np.abs(e) ** 2).sum(axis=0))[:, 0, :]      # (nz, nx)
    db = 20 * np.log10(np.clip(amp / amp.max(), 1e-6, None))
    X, Z = np.meshgrid(x * 1e3, z * 1e3)                      # mm

    fig, ax = plt.subplots(figsize=(8.6, 4.8), dpi=140)
    pc = ax.pcolormesh(X, Z, db, cmap="inferno", vmin=-40, vmax=0,
                       shading="auto")
    ax.axhline(0.0, color="#1f5fa8", lw=1.0, ls="-")
    ax.axhline(am.SUBSTRATE_H_MM, color="#1f5fa8", lw=1.0, ls="-")
    for i in range(am.N_X):
        x0, _, x1, _ = am.patch_rect(i, 1)
        ax.plot([x0, x1], [am.SUBSTRATE_H_MM] * 2, color="#f1c40f", lw=2.6,
                solid_capstyle="butt")
    fig.colorbar(pc, ax=ax, label="|E| (dB, 0 dB = maximum)")
    ax.set_xlabel("x (mm)")
    ax.set_ylabel("z (mm)")
    ax.set_title("Broadside array, 77 GHz - |E| in the xz plane (E plane)")
    ax.set_aspect("equal")
    fig.tight_layout()
    out = FIG / "fig13_nearfield_xz.png"
    fig.savefig(out)
    plt.close(fig)
    print(f"wrote {out}")

    # --- E_near (horizontal plane): |E| over xy ---------------------------
    e, x, y, z = _read_dump(dump / "E_near.h5")
    amp = np.sqrt((np.abs(e) ** 2).sum(axis=0))[0, :, :]      # (ny, nx)
    db = 20 * np.log10(np.clip(amp / amp.max(), 1e-6, None))
    X, Y = np.meshgrid(x * 1e3, y * 1e3)

    fig, ax = plt.subplots(figsize=(7.4, 6.6), dpi=140)
    pc = ax.pcolormesh(X, Y, db, cmap="inferno", vmin=-40, vmax=0,
                       shading="auto")
    for j in range(am.N_Y):
        for i in range(am.N_X):
            x0, y0, x1, y1 = am.patch_rect(i, j)
            ax.add_patch(plt.Rectangle((x0, y0), x1 - x0, y1 - y0,
                                       fill=False, edgecolor="#2ecc71",
                                       linewidth=1.0))
    fig.colorbar(pc, ax=ax, label="|E| (dB, 0 dB = maximum)")
    ax.set_xlabel("x (mm)")
    ax.set_ylabel("y (mm)")
    ax.set_title(f"Broadside array, 77 GHz - |E| at z = {z[0] * 1e3:.2f} mm")
    ax.set_aspect("equal")
    fig.tight_layout()
    out = FIG / "fig14_nearfield_top.png"
    fig.savefig(out)
    plt.close(fig)
    print(f"wrote {out}")
    return {"nearfield_z_mm": float(z[0] * 1e3)}


# ---------------------------------------------------------------------------
# accepted-power correction for the driven arrays
# ---------------------------------------------------------------------------

def array_power_numbers(work: str, nf: dict) -> dict:
    """Corrected eta/gain: Pin = accepted power of ALL driven ports."""
    freqs = nf["freq"]
    p_acc = accepted_power(work, freqs)
    fi = int(np.argmin(np.abs(freqs - F77)))
    eta = nf["Prad"] / p_acc
    gain = nf["D"] + 10.0 * np.log10(np.clip(eta, 1e-9, None))[:, None, None]
    return {"p_acc": p_acc, "eta": eta, "gain": gain, "fi": fi}


# ---------------------------------------------------------------------------
# main
# ---------------------------------------------------------------------------

def main() -> int:
    FIG.mkdir(parents=True, exist_ok=True)
    results: dict = {"design": {
        "f0_ghz": am.F0_HZ / 1e9, "band_ghz": [am.F_START_HZ / 1e9,
                                               am.F_STOP_HZ / 1e9],
        "patch_l_mm": am.PATCH_L_MM, "patch_w_mm": am.PATCH_W_MM,
        "feed_inset_mm": am.FEED_INSET_MM,
        "pitch_mm": am.PITCH_MM, "board_mm": am.BOARD_MM,
        "n_x": am.N_X, "n_y": am.N_Y,
        "substrate": am.SUBSTRATE_NAME, "eps_r": am.SUBSTRATE_EPS_R,
        "tand": am.SUBSTRATE_TAND, "h_mm": am.SUBSTRATE_H_MM,
        "lambda0_mm": am.LAMBDA0_MM,
        "scan_theta_deg": am.SCAN30_THETA_DEG,
    }}

    results["element"] = fig_element_zin()
    results["retune"] = fig_retune_journey()
    results["element"].update(fig_element_pattern())
    results["element"].update(fig_element_vs_embedded())
    results["coupling"] = fig_coupling()

    results["broadside"] = fig_array_pattern(
        RUNS["broadside"], "Broadside", "fig10_broadside.png", None)
    results["scan30"] = fig_array_pattern(
        RUNS["scan30"], f"Scan +{am.SCAN30_THETA_DEG:g}$^\\circ$",
        "fig11_scan30.png", am.scan_phases_deg(am.SCAN30_THETA_DEG))
    results["scan_summary"] = fig_scan_summary(results["broadside"],
                                               results["scan30"])
    results["fields"] = fig_fields()

    # corrected efficiencies for the driven arrays (see module docstring)
    for key in ("broadside", "scan30"):
        nf = load_nf(RUNS[key])
        corr = array_power_numbers(RUNS[key], nf)
        fi = corr["fi"]
        results[key]["eff77_corr"] = float(corr["eta"][fi])
        results[key]["gain77_corr"] = float(corr["gain"][fi][
            int(np.argmin(np.abs(nf["theta"] - results[key]["theta_peak_deg"]))),
            int(np.argmin(np.abs(nf["phi"] - results[key]["phi_peak_deg"])))])
        results[key]["p_acc77"] = float(corr["p_acc"][fi])
        results[key]["prad77"] = float(nf["Prad"][fi])
        print(f"{key}: eta(77) raw {results[key]['eff77_raw']:.3f} "
              f"-> corrected {corr['eta'][fi]:.3f}  (Pin from "
              f"{16} driven ports)")

    # cross-check: single-driven runs must agree with the bundled numbers
    nf_e = load_nf(RUNS["element"])
    p_e = accepted_power(RUNS["element"], nf_e["freq"])
    print("element: bundled eff", [round(v, 3) for v in nf_e["eff"]],
          " recomputed", [round(v, 3) for v in nf_e["Prad"] / p_e])

    # before-design reference numbers
    f_b, s_b = load_s(RUNS["before"])
    k = s77_index(f_b)
    nf_b = load_nf(RUNS["before"])
    results["before"] = {
        "zin77_real": float(z_of(f_b, s_b)[k].real),
        "zin77_imag": float(z_of(f_b, s_b)[k].imag),
        "s11_77_db": float(s11_db(f_b, s_b)[k]),
        "dmax77": float(nf_b["Dmax"][f77_index(nf_b)]),
        "eff77": float(nf_b["eff"][f77_index(nf_b)]),
    }

    # mesh line counts from the validator output (if present)
    meshes: dict[str, list[int]] = {}
    for key, vdir in (
            ("element", "work_validate_patch_element_77ghz"),
            ("embedded", "work_validate_patch_array_4x4_77ghz"),
            ("broadside", "work_validate_patch_array_4x4_broadside"),
            ("scan30", "work_validate_patch_array_4x4_scan30")):
        p = HERE / vdir / "grid_lines.json"
        if p.is_file():
            g = json.loads(p.read_text(encoding="utf-8"))
            meshes[key] = [len(g[k]) for k in ("x", "y", "z")]
    results["mesh"] = meshes

    out = HERE / "docs" / "results.json"
    out.write_text(json.dumps(results, indent=2, ensure_ascii=False) + "\n",
                   encoding="utf-8")
    print(f"wrote {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
