"""Write a wideband diagnostic element project (47..107 GHz excitation).

Why this exists
---------------
The production band is 74..80 GHz, which makes the Gaussian excitation only
``fc = 3 GHz`` wide: the port signal carries useful energy in a 6 GHz window
alone, so ``analyze_zin.py`` cannot see the resonance if it lies outside it.
The first element run measured Zin(77 GHz) = 14.7 + 68j (far from a match),
and 74..80 GHz showed no dip at all - so the real resonance may sit anywhere.

Widening the mesh frequency range grows the excitation bandwidth while
keeping everything else identical:

=========================  ==========================================
``f0`` (excitation centre)  ``target_frequency`` = 77 GHz  (unchanged)
``fc`` (-20 dB corners)     0.5 * (stop - start) = 30 GHz (was 3 GHz)
S-parameter sweep          linspace(47, 107, 201), 300 MHz step
mesh                       unchanged: lambda_air(107 GHz)/10 = 0.28 mm
                           > ``max_cellsize`` 0.25 mm, so the validated
                           grid (1.73 M cells) is kept
=========================  ==========================================

Usage::

    python make_wideband.py
    python run_headless.py patch_element_77ghz_wideband.fcsim.json
    python analyze_zin.py work_element_wideband --print --fmin 45e9 --fmax 110e9
"""

from __future__ import annotations

import json
from pathlib import Path

import array_model as am

HERE = Path(__file__).resolve().parent

STEM = "patch_element_77ghz_wideband"
WB_START_HZ = 47.0e9
WB_STOP_HZ = 107.0e9


def main() -> int:
    data = am.project_dict(name=STEM, variant="element")
    data["mesh"]["start_frequency"] = WB_START_HZ
    data["mesh"]["stop_frequency"] = WB_STOP_HZ

    path = HERE / f"{STEM}.fcsim.json"
    path.write_text(json.dumps(data, indent=2, ensure_ascii=False) + "\n",
                    encoding="utf-8")
    print(f"wrote {path}  ({path.stat().st_size} bytes)")
    print(f"  excitation f0 = {data['mesh']['target_frequency']/1e9:.0f} GHz"
          f"   fc = {0.5*(WB_STOP_HZ - WB_START_HZ)/1e9:.0f} GHz")
    print(f"  sweep {WB_START_HZ/1e9:.0f}..{WB_STOP_HZ/1e9:.0f} GHz,"
          " mesh parameters untouched (validated grid is kept)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
