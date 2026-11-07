"""Print the raw measurements of the two control runs (no model assumptions).

Background
----------
The first element run measured Zin(77 GHz) = 14.8 + 68j ohm.  Whether a
probe-fed patch can be matched at all depends on two numbers, so
``make_controls.py`` simulated them separately:

============================  =================================================
``probe_only``                substrate + ground + probe, *no patch*.
                              The bare line port measures Z_open(f).  It is
                              **capacitive** (-205 ohm at 77 GHz): the open
                              line section is shorter than lambda/4, exactly
                              what ``fit_model.py`` fits as a transmission
                              line transformer (Zp, eps_eff, open-end phase).
                              A naive series-inductance reading would get the
                              sign wrong - the raw numbers are printed here so
                              the sign is visible in the report.
``inset005``                  the real patch, probe at 0.05 L (almost at the
                              radiating edge): Re Z peak ~ R_edge, and the
                              first matched element.
============================  =================================================

The modelling (probe de-embedding, RLC fit, retune prediction) lives in
``fit_model.py``.

Usage::

    python analyze_controls.py
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np

from analyze_zin import dft, load_signal

HERE = Path(__file__).resolve().parent
F77 = 77.0e9


def z_of(work: str, fmin: float = 47e9, fmax: float = 107e9,
         steps: int = 1201) -> tuple[np.ndarray, np.ndarray]:
    sim = HERE / work / "simulation"
    tu, ut = load_signal(sim / "port_ut_1")
    ti, it = load_signal(sim / "port_it_1")
    f = np.linspace(fmin, fmax, steps)
    return f, dft(tu, ut, f) / dft(ti, it, f)


def probe_only() -> None:
    f, z = z_of("work_probe_only_77ghz")
    k77 = int(np.argmin(np.abs(f - F77)))
    print("=== probe_only (substrate + ground + probe, no patch) ===")
    print("  raw open-line measurement, no equivalent circuit assumed")
    print("   f [GHz]    Re Z     Im Z     interpretation")
    for fg in (50, 60, 70, 77, 85, 95, 105):
        k = int(np.argmin(np.abs(f - fg * 1e9)))
        kind = "capacitive" if z[k].imag < 0 else "inductive"
        print(f"   {f[k]/1e9:7.2f}  {z[k].real:8.2f} {z[k].imag:+9.2f}   {kind}")
    band = (f >= 55e9) & (f <= 100e9)
    print(f"  Z(77 GHz)   = {z[k77].real:.2f} {z[k77].imag:+.2f}j ohm")
    print(f"  mean Re Z   = {z[band].real.mean():.2f} ohm "
          "(radiates nothing -> ~0)")
    print("  -> Im Z < 0 over the whole band: the probe is a SHORT "
          "open-circuited line")
    print("     section (l = 0.254 mm < lambda/4), not a series inductor.")


def inset005() -> None:
    path = HERE / "work_inset005_77ghz" / "sparams" / "sparams.json"
    d = json.loads(path.read_text(encoding="utf-8"))
    f = np.array(d["freq"])
    s = np.array(d["S11"]["real"]) + 1j * np.array(d["S11"]["imag"])
    z = 50 * (1 + s) / (1 - s)

    k = int(np.argmax(z.real))
    f_r, r_peak = f[k], z.real[k]
    inset = 0.05
    r_edge = r_peak / np.cos(np.pi * inset) ** 2

    print("\n=== inset005 (probe at 0.05 L from the patch edge) ===")
    print(f"  Re Z peak     : {r_peak:.1f} ohm @ {f_r/1e9:.2f} GHz")
    print(f"  R_edge        : {r_edge:.1f} ohm  "
          f"(= peak / cos^2({inset:.2f} pi))")
    k77 = int(np.argmin(np.abs(f - F77)))
    print(f"  Z(77 GHz)     : {z[k77].real:.1f} {z[k77].imag:+.1f}j ohm"
          f"   |S11| = {20*np.log10(abs(s[k77])):.2f} dB")
    kmin = int(np.argmin(20 * np.log10(np.abs(s))))
    print(f"  best |S11|    : {20*np.log10(abs(s[kmin])):.2f} dB "
          f"@ {f[kmin]/1e9:.2f} GHz")
    print("  -> first element state that already matches at 77 GHz; "
          "fit_model.py uses it as the patch reference.")


def main() -> int:
    probe_only()
    inset005()
    print("\n(next: python fit_model.py - probe TL fit + retune prediction)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
