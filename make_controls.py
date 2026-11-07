"""Control simulations that isolate the two unknowns of the element feed.

Background
----------
The first element run gave Zin(77 GHz) = 14.8 + 68j ohm - strongly
inductive, and the wideband sweep (47..107 GHz) shows the input reactance
staying positive over the whole band while Re Z peaks at only 42 ohm
(@ 69 GHz).  Two numbers decide whether a probe-fed patch can be matched
at all, and neither is known yet:

``L_probe``   series inductance of the 0.254 mm long, ~0.1 mm wide probe.
              X(77 GHz) = 2*pi*f*L is 70..90 ohm for L = 0.15..0.18 nH -
              larger than the reactance the patch can produce.

``R_edge``    resistance of the patch at its radiating edge.  The classic
              inset law R(s) = R_edge*cos^2(pi*s/L) turns the measured
              42 ohm (at s = 0.30 L) into R_edge = 60 ohm, but that is an
              inference; theory predicts ~150 ohm.  The value decides
              everything (an edge-fed patch can only cancel the probe if
              R_edge is large).

Two runs measure them directly:

===========================  ==========================================
``probe_only``               substrate + ground + probe, *no patch*.
                             Zin = j*w*L_probe (+ stray tip capacitance),
                             so Im Z is a straight line through the origin
                             with slope 2*pi*L_probe.
``inset005``                 same patch, probe moved to s = 0.05*L.
                             Re Z peak then is 0.976 * R_edge, i.e. the
                             edge resistance itself.
===========================  ==========================================

Both use the wideband 47..107 GHz excitation of ``make_wideband.py``.

Usage::

    python make_controls.py
    python run_headless.py probe_only_77ghz.fcsim.json      # ~2 min
    python run_headless.py inset005_77ghz.fcsim.json        # ~3 min
"""

from __future__ import annotations

import json
from pathlib import Path

import array_model as am

HERE = Path(__file__).resolve().parent

#: wideband excitation (kept identical to make_wideband.py)
WB_START_HZ = 47.0e9
WB_STOP_HZ = 107.0e9


def _wideband(data: dict) -> dict:
    data["mesh"]["start_frequency"] = WB_START_HZ
    data["mesh"]["stop_frequency"] = WB_STOP_HZ
    return data


def _write(stem: str, data: dict) -> Path:
    path = HERE / f"{stem}.fcsim.json"
    path.write_text(json.dumps(data, indent=2, ensure_ascii=False) + "\n",
                    encoding="utf-8")
    print(f"wrote {path}  ({path.stat().st_size} bytes)")
    return path


def build_probe_only() -> Path:
    """Substrate + ground + probe, no patch -> measures L_probe."""
    data = am.project_dict(name="probe_only_77ghz", variant="element")
    data["primitives"] = [p for p in data["primitives"]
                          if not p["name"].startswith("patch")]
    # nothing radiates: the NF2FF box would only transform noise
    data["nf2ff"]["enabled"] = False
    data["sim"]["nf2ff_enabled"] = False
    return _write("probe_only_77ghz", _wideband(data))


def build_inset005() -> Path:
    """Probe at 5 % of L from the patch edge -> Re Z peak = R_edge."""
    saved = am.FEED_INSET_MM
    am.FEED_INSET_MM = 0.05 * am.PATCH_L_MM
    try:
        data = am.project_dict(name="inset005_77ghz", variant="element")
    finally:
        am.FEED_INSET_MM = saved
    return _write("inset005_77ghz", _wideband(data))


def main() -> int:
    print(f"probe: inset was {am.FEED_INSET_MM:.4f} mm "
          f"({am.FEED_INSET_MM / am.PATCH_L_MM:.2f} L)")
    build_probe_only()
    build_inset005()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
