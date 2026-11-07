"""Run the whole study end to end (projects, validation, FDTD, report).

Every step is skipped when its output already exists, so the script can be
re-run after an interruption.  Use ``--force`` to recompute everything.

Usage::

    python run_all.py
    python run_all.py --force --no-sim     # only figures + report
"""

from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
PY = sys.executable

PROJECTS = [
    "patch_element_77ghz.fcsim.json",
    "patch_array_4x4_77ghz.fcsim.json",
    "patch_array_4x4_broadside.fcsim.json",
    "patch_array_4x4_scan30.fcsim.json",
]

#: project stem -> work directory written by run_headless.py
WORK = {
    "patch_element_77ghz": "work_element",
    "patch_array_4x4_77ghz": "work_array_embedded",
    "patch_array_4x4_broadside": "work_array_broadside",
    "patch_array_4x4_scan30": "work_array_scan30",
}


def run(step: str, *args: str) -> None:
    print(f"\n=== {step} ===", flush=True)
    subprocess.run([PY, *args], cwd=str(HERE), check=True)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--force", action="store_true",
                    help="recompute even if the output exists")
    ap.add_argument("--no-sim", action="store_true",
                    help="skip the FDTD runs")
    args = ap.parse_args()

    if args.force or not (HERE / "docs/figures/fig04_array3d.png").is_file():
        run("projects", "make_project.py", "--preview")
    if args.force or not (HERE / "work_validate_patch_element_77ghz").is_dir():
        run("validation", "validate_project.py")

    if not args.no_sim:
        for proj in PROJECTS:
            stem = proj.split(".")[0]
            work = WORK[stem]
            if args.force or not (HERE / work / "nf2ff" / "nf2ff.json").is_file():
                run(f"FDTD {proj}", "run_headless.py", proj)

    if args.force or not (HERE / "docs/results.json").is_file():
        run("figures", "plot_results.py")
    run("report", "build_docs.py")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
