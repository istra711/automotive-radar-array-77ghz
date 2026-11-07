"""Write one VisualopenEMS project per beam-steering angle.

The projects are the input of the full-wave beam-steering sweep that feeds
``make_animation.py``.  Everything except the port phases is identical to the
production array, so the runs are directly comparable.

Usage::

    python make_sweep.py
    python make_sweep.py --angles -30 0 30
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import array_model as am

HERE = Path(__file__).resolve().parent


def project_path(theta_deg: float) -> Path:
    return HERE / f"scan_{am.sweep_tag(theta_deg)}_77ghz.fcsim.json"


def write_project(theta_deg: float) -> Path:
    tag = am.sweep_tag(theta_deg)
    data = am.project_dict(
        f"scan_{tag}_77ghz",
        "scan",
        theta_deg=float(theta_deg),
    )
    path = project_path(theta_deg)
    path.write_text(json.dumps(data, indent=2, ensure_ascii=False) + "\n",
                    encoding="utf-8")
    return path


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--angles", type=float, nargs="*",
                    default=list(am.SWEEP_ANGLES_DEG),
                    help="steering angles in degrees (default: the full sweep)")
    args = ap.parse_args()

    print(am.summary())
    for theta in args.angles:
        p = write_project(theta)
        ph = am.scan_phases_deg(theta)
        step = ph[2] - ph[1]
        print(f"  theta {theta:+6.1f} deg -> phase step {step:+7.2f} deg/px  "
              f"({p.name})")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
