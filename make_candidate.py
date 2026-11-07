"""Build a candidate element project for one (L, inset) retune step.

Retune loop
-----------
``fit_model.py`` linearises the two measured element runs (inset 0.05 L and
0.30 L) around the nominal design and predicts the ``(L, inset)`` that puts
Zin on 50 ohm at 77 GHz.  Every prediction is checked by one cheap wideband
run (~2.5 min), then the model is re-linearised with the new measurement::

    python fit_model.py                       # predict (L, inset)
    python make_candidate.py --L 1.1610 --s 0.05
    python run_headless.py cand_L1p1610_s0p050.fcsim.json
    python fit_model.py                       # re-linearise, verify

The candidate is written with the wideband 47..107 GHz excitation of
``make_wideband.py`` (mesh parameters untouched, validated grid kept) unless
``--narrowband`` selects the production 74..80 GHz band.

Usage::

    python make_candidate.py --L 1.1610 --s 0.05
    python make_candidate.py --L 1.1683 --s 0.02 --narrowband
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import array_model as am

HERE = Path(__file__).resolve().parent

#: wideband excitation (kept identical to make_wideband.py)
WB_START_HZ = 47.0e9
WB_STOP_HZ = 107.0e9


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--L", type=float, required=True,
                    help="patch length in mm (resonant dimension, x)")
    ap.add_argument("--s", type=float, required=True,
                    help="feed inset as a fraction of L (0 = patch edge)")
    ap.add_argument("--name", default="",
                    help="project name (default: cand_L<s>_s<s>)")
    ap.add_argument("--narrowband", action="store_true",
                    help="production band 74..80 GHz instead of 47..107 GHz")
    args = ap.parse_args()

    if not (0.0 <= args.s <= 0.45):
        raise SystemExit("--s must lie in [0, 0.45]")
    if not (0.8 <= args.L <= 1.6):
        raise SystemExit("--L must lie in [0.8, 1.6] mm (77 GHz patch)")

    # override the design in the single source of truth for this build only
    am.PATCH_L_MM = args.L
    am.FEED_INSET_MM = args.s * args.L

    name = args.name or (f"cand_L{args.L:.4f}_s{args.s:.3f}"
                         .replace(".", "p"))
    data = am.project_dict(name=name, variant="element")
    if not args.narrowband:
        data["mesh"]["start_frequency"] = WB_START_HZ
        data["mesh"]["stop_frequency"] = WB_STOP_HZ

    path = HERE / f"{name}.fcsim.json"
    path.write_text(json.dumps(data, indent=2, ensure_ascii=False) + "\n",
                    encoding="utf-8")
    band = ("47..107 GHz" if not args.narrowband else "74..80 GHz")
    print(f"wrote {path}  ({path.stat().st_size} bytes)")
    print(f"  patch   : L = {args.L:.4f} mm x {am.PATCH_W_MM:.4f} mm")
    print(f"  feed    : inset {am.FEED_INSET_MM:.4f} mm ({args.s:.3f} L) "
          f"from the -x edge")
    print(f"  excitation band: {band}")
    print(f"  run     : python run_headless.py {path.name}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
