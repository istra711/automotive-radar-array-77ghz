"""Fit the element feed model and predict the (L, inset) that matches at 77 GHz.

Why a model at all
------------------
Three measurements are available:

``work_probe_only_77ghz``  probe with an *open* top  -> Z_open(f)
``work_inset005_77ghz``    patch, probe at 0.05 L    -> Z_total(f)
``work_patch_element_77ghz_wideband``  patch, probe at 0.30 L -> Z_total(f)

The probe is not a lumped inductance: the open-top control measures
-205 ohm (capacitive) at 77 GHz, which is exactly what a short
open-circuited line section looks like (l = 0.254 mm < lambda/4).  So the
feed has to be treated as a transmission-line transformer

    Z_total = Zp * (Z_patch + j*Zp*tan(th)) / (Zp + j*Z_patch*tan(th))

with a characteristic impedance ``Zp`` and electrical length ``th(f)``
fitted from the open-top control.  The patch itself is a parallel RLC
near its fundamental (R_peak at f_r, quality factor Q) - the classic
cavity/transmission-line model of a rectangular patch.

With the probe characterised, the patch curves of both inset runs are
de-embedded and fitted, and the (L, inset) that transforms to 50 ohm at
77 GHz is solved numerically.

Usage::

    python fit_model.py
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
from scipy.optimize import least_squares

from analyze_zin import dft, load_signal

HERE = Path(__file__).resolve().parent
F77 = 77.0e9
L_PROBE_MM = 0.254          # probe length = substrate thickness
PATCH_L_MM = 1.1383         # nominal design length (array_model.PATCH_L_MM)


# ---------------------------------------------------------------------------
# measurements
# ---------------------------------------------------------------------------

def z_from_signals(work: str, fmin=47e9, fmax=107e9, steps=1201):
    sim = HERE / work / "simulation"
    tu, ut = load_signal(sim / "port_ut_1")
    ti, it = load_signal(sim / "port_it_1")
    f = np.linspace(fmin, fmax, steps)
    return f, dft(tu, ut, f) / dft(ti, it, f)


def z_from_sparams(work: str):
    d = json.loads((HERE / work / "sparams" / "sparams.json")
                   .read_text(encoding="utf-8"))
    f = np.array(d["freq"])
    s = np.array(d["S11"]["real"]) + 1j * np.array(d["S11"]["imag"])
    return f, 50 * (1 + s) / (1 - s)


# ---------------------------------------------------------------------------
# probe: open-circuited line section
# ---------------------------------------------------------------------------

def probe_theta(f: np.ndarray, eps_p: float, phi: float) -> np.ndarray:
    """Electrical length incl. an open-end phase correction (rad)."""
    return 2 * np.pi * f * np.sqrt(eps_p) * L_PROBE_MM * 1e-3 / 299_792_458.0 \
        + phi


def fit_probe(f, z) -> tuple[float, float, float]:
    """Z_open = -j*Zp*cot(th)  ->  (Zp, eps_p, phi)."""
    band = (f >= 50e9) & (f <= 105e9)

    def res(p):
        zp, eps_p, phi = p
        return (-zp / np.tan(probe_theta(f[band], eps_p, phi))
                - z.imag[band]) / 100.0

    r = least_squares(res, [120.0, 2.0, 0.0],
                      bounds=([20, 1.0, -1.5], [400, 4.0, 1.5]))
    return float(r.x[0]), float(r.x[1]), float(r.x[2])


def probe_transform(z_load, f, zp, eps_p, phi):
    """Input impedance of the probe line loaded with z_load."""
    t = np.tan(probe_theta(f, eps_p, phi))
    return zp * (z_load + 1j * zp * t) / (zp + 1j * z_load * t)


def probe_deembed(z_in, f, zp, eps_p, phi):
    """Inverse of :func:`probe_transform` (z_load from z_in)."""
    t = np.tan(probe_theta(f, eps_p, phi))
    return zp * (z_in - 1j * zp * t) / (zp - 1j * z_in * t)


# ---------------------------------------------------------------------------
# patch: parallel RLC near the fundamental
# ---------------------------------------------------------------------------

def patch_zlc(f, r_peak, f_r, q):
    w = 2 * np.pi * f
    w0 = 2 * np.pi * f_r
    g = 1.0 / r_peak
    b = (q * g) * (f / f_r - f_r / f)
    return 1.0 / (g + 1j * b)


def fit_patch(f, z_patch, band):
    def res(p):
        r_peak, f_r, q = p
        zz = patch_zlc(f[band], *p)
        return np.concatenate([(zz.real - z_patch[band].real),
                               (zz.imag - z_patch[band].imag)]) / 50.0

    r = least_squares(res, [60.0, 70e9, 6.0],
                      bounds=([5, 40e9, 0.5], [2000, 160e9, 60.0]))
    return float(r.x[0]), float(r.x[1]), float(r.x[2])


# ---------------------------------------------------------------------------
# main
# ---------------------------------------------------------------------------

def main() -> int:
    f_o, z_o = z_from_signals("work_probe_only_77ghz")
    zp, eps_p, phi = fit_probe(f_o, z_o)
    print("=== probe (open top) ===")
    print(f"  Zp = {zp:.1f} ohm, eps_eff = {eps_p:.2f}, "
          f"open-end phase = {np.degrees(phi):.1f} deg")
    print(f"  th(77 GHz) = {np.degrees(probe_theta(np.array([F77]), eps_p, phi)[0]):.1f} deg")
    pred = -zp / np.tan(probe_theta(f_o, eps_p, phi))
    err = np.abs(pred - z_o.imag) / np.abs(z_o.imag)
    print(f"  fit error on Im Z: mean {100*err.mean():.1f} %, "
          f"max {100*err.max():.1f} %")

    runs = [("inset 0.30 L (nominal)", "work_patch_element_77ghz_wideband", 0.30),
            ("inset 0.05 L", "work_inset005_77ghz", 0.05)]

    curves: dict[float, tuple[np.ndarray, np.ndarray]] = {}
    for label, work, inset in runs:
        f, z_tot = z_from_sparams(work)
        z_patch = probe_deembed(z_tot, f, zp, eps_p, phi)
        curves[inset] = (f, z_patch)
        band = (f >= 55e9) & (f <= 100e9)
        r_peak, f_r, q = fit_patch(f, z_patch, band)
        k = int(np.argmin(np.abs(f - F77)))
        # cos^2 inset law -> edge resistance
        r_edge = r_peak / np.cos(np.pi * inset) ** 2
        print(f"\n=== {label}: {work} ===")
        print(f"  Z_total(77) = {z_tot[k].real:6.1f} {z_tot[k].imag:+6.1f}j ohm")
        print(f"  Z_patch(77) = {z_patch[k].real:6.1f} {z_patch[k].imag:+6.1f}j ohm")
        print(f"  RLC fit     : R_peak = {r_peak:.1f} ohm @ {f_r/1e9:.2f} GHz, "
              f"Q = {q:.1f}")
        print(f"  -> R_edge    = {r_edge:.1f} ohm")

    # ------------------------------------------------------------------
    # prediction: Newton step in (inset fraction s, ln L) from the two
    # measured patch curves, targeting Z_port(77) = 50 ohm.
    #
    #   * dz/ds        : secant between the two measured insets at 77 GHz
    #                    (measured, no cos^2 assumption)
    #   * dz/dlnL      : f * dZ/df of the inset005 curve - changing L by
    #                    dlnL shifts the whole resonance curve to
    #                    f*(1 + dlnL), a standard detuning argument
    # ------------------------------------------------------------------
    z_need = probe_deembed(np.array([50.0 + 0j]), np.array([F77]),
                           zp, eps_p, phi)[0]
    print("\n=== prediction: Z_port(77 GHz) = 50 ohm ===")
    print(f"  required Z_patch(feed) = {z_need.real:6.1f} "
          f"{z_need.imag:+6.1f}j ohm")

    (f1, zp1), (f2, zp2) = curves[0.05], curves[0.30]

    def at_77(f, z):
        return (float(np.interp(F77, f, z.real))
                + 1j * float(np.interp(F77, f, z.imag)))

    z_a, z_b = at_77(f1, zp1), at_77(f2, zp2)
    dz_ds = (z_b - z_a) / (0.30 - 0.05)             # per unit inset fraction

    band = (f1 >= 74e9) & (f1 <= 80e9)
    fr_ghz = f1[band] / 1e9
    pr = np.polyfit(fr_ghz, zp1[band].real, 1)
    pq = np.polyfit(fr_ghz, zp1[band].imag, 1)
    dz_dlnl = F77 * (pr[0] + 1j * pq[0]) * 1e-9     # dz/dlnL = f * dz/df

    print(f"  measured Z_patch(77): s=0.05 -> {z_a.real:6.1f} {z_a.imag:+6.1f}j"
          f"   s=0.30 -> {z_b.real:6.1f} {z_b.imag:+6.1f}j")
    print(f"  dZ/ds (secant)      = {dz_ds.real:7.1f} {dz_ds.imag:+7.1f}j")
    print(f"  dZ/dlnL (= f dZ/df) = {dz_dlnl.real:7.1f} {dz_dlnl.imag:+7.1f}j")

    mat = np.array([[dz_ds.real, dz_dlnl.real],
                    [dz_ds.imag, dz_dlnl.imag]])
    rhs = np.array([z_need.real - z_a.real, z_need.imag - z_a.imag])
    ds, dlnl = np.linalg.solve(mat, rhs)
    s_new = 0.05 + ds
    l_new = PATCH_L_MM * float(np.exp(dlnl))
    inset_new = s_new * l_new

    print(f"\n  solution: L = {l_new:.4f} mm ({100*dlnl:+.2f} %), "
          f"inset = {s_new:.3f} L = {inset_new:.4f} mm")
    if not (0.03 <= s_new <= 0.45) or abs(dlnl) > 0.10:
        print("  !! outside the validated range - treat as first guess only")
    print("  -> next step: simulate this candidate (make_candidate.py) and"
          " re-linearise")

    # predicted |S11|(f) of the candidate over the production band
    def z_patch_model(f_hz, s, l_scale):
        """Measured curves, L-scaled by f -> f*l_scale, interpolated in s."""
        out = []
        for ff in f_hz:
            fsh = ff * l_scale                       # base-curve frequency
            v1 = (float(np.interp(fsh, f1, zp1.real))
                  + 1j * float(np.interp(fsh, f1, zp1.imag)))
            v2 = (float(np.interp(fsh, f2, zp2.real))
                  + 1j * float(np.interp(fsh, f2, zp2.imag)))
            w = (s - 0.05) / (0.30 - 0.05)
            out.append(v1 + (v2 - v1) * w)
        return np.array(out)

    f_pred = np.linspace(74e9, 80e9, 61)
    z_pred = probe_transform(z_patch_model(f_pred, s_new,
                                           l_new / PATCH_L_MM),
                             f_pred, zp, eps_p, phi)
    s11_pred = 20 * np.log10(np.abs((z_pred - 50) / (z_pred + 50)))
    kmin = int(np.argmin(s11_pred))
    print(f"\n  predicted candidate: Z(77) = {z_pred[30].real:5.1f} "
          f"{z_pred[30].imag:+5.1f}j ohm, "
          f"best |S11| = {s11_pred[kmin]:.1f} dB @ {f_pred[kmin]/1e9:.2f} GHz")
    print(f"  candidate file      : cand_L{l_new:.4f}_s{s_new:.3f}.fcsim.json")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
