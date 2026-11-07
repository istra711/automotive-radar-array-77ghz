"""4x4 patch antenna array at 77 GHz on RT/duroid 5880 - single source of truth.

Everything the project needs is derived here:

* the element design (Hammerstad patch equations),
* the 4x4 array layout with lambda0/2 pitch,
* the VisualopenEMS project dictionaries (``make_project.py``),
* the excitation variants (element / embedded / broadside / scan +30 deg).

Drawing unit is **mm**, frequencies are **Hz**.

Geometry (all variants share the same board)
-------------------------------------------
The board is a 12 x 12 mm RT/duroid 5880 slab (eps_r = 2.2, tan_d = 0.0009,
h = 0.254 mm) with a full copper ground on the bottom face (z = 0).  The
patches are zero-thickness PEC sheets on the top face (z = h): the resonant
length **L runs along x**, the radiating width **W along y**, so the far field
is x-polarised and the E-plane is the xz-plane.  Each patch is probe fed: a
vertical 50 ohm lumped port from the ground plane to the patch, inset by
``FEED_INSET_MM`` from the -x edge (the standard inset match).

Array variants
--------------
===========  ===========================================================
element      one patch at the board centre, P1 excited
embedded     16 patches, only P6 (an inner element) excited, the other
             15 ports are 50 ohm terminations -> embedded element pattern
             + the first column of the S matrix (mutual coupling)
broadside    all 16 ports excited in phase -> true 4x4 broadside pattern
scan30       all 16 ports excited with a progressive phase for a beam
             steered to theta = +30 deg in the E-plane (xz)
===========  ===========================================================
"""

from __future__ import annotations

import math
from typing import Any

# ---------------------------------------------------------------------------
# Physical constants
# ---------------------------------------------------------------------------

C0 = 299_792_458.0
EPS0 = 8.854_187_812_8e-12

# ---------------------------------------------------------------------------
# Design parameters
# ---------------------------------------------------------------------------

#: Automotive radar centre frequency (ETSI 76-77 GHz band)
F0_HZ = 77.0e9
#: Simulated band (S-parameters are evaluated here)
F_START_HZ = 74.0e9
F_STOP_HZ = 80.0e9
#: 20 dB corner frequency of the Gaussian excitation
FC_HZ = 0.5 * (F_STOP_HZ - F_START_HZ)

#: Substrate: Rogers RT/duroid 5880 - the low-loss reference for 77 GHz radar
SUBSTRATE_NAME = "rt_duroid_5880"
SUBSTRATE_EPS_R = 2.2
SUBSTRATE_TAND = 0.0009
SUBSTRATE_H_MM = 0.254          # 10 mil
COPPER_UM = 35.0                # modelled as ideal PEC sheets

#: Array layout
N_X = 4                         # elements along x (E-plane)
N_Y = 4                         # elements along y (H-plane)
BOARD_MM = 12.0                 # square ground plane / substrate edge

#: Air margin between structure and PML: [xmin, xmax, ymin, ymax, zmin, zmax]
AIR_MARGIN_MM = [4.0, 4.0, 4.0, 4.0, 4.0, 6.0]

#: NF2FF integration box (manual, inside the domain, outside the board)
NF2FF_START_MM = [-8.0, -8.0, -2.5]
NF2FF_STOP_MM = [8.0, 8.0, 5.5]
NF2FF_CENTER_MM = [0.0, 0.0, SUBSTRATE_H_MM / 2.0]
NF2FF_FREQS_HZ = [74e9, 75e9, 76e9, 77e9, 78e9, 79e9, 80e9]
NF2FF_THETA_STEPS = 91           # 0..180 deg in 2 deg steps
NF2FF_PHI_STEPS = 181            # 0..360 deg in 2 deg steps

#: Solver
END_CRITERIA = 1e-4              # -40 dB energy criterion
MAX_TIMESTEPS = 40_000           # safety cap (the criterion stops earlier)

#: E-field dumps (frequency domain at f0) used for the field-slice figures
DUMPS = [
    {
        "name": "E_near",
        "quantity": "E",
        "start": [-8.0, -8.0, 1.0],
        "stop": [8.0, 8.0, 1.0],
    },
    {
        "name": "E_xz",
        "quantity": "E",
        "start": [-9.5, 0.0, -3.5],
        "stop": [9.5, 0.0, 6.0],
    },
]

#: Excitation variant of each generated project file
VARIANTS = ("element", "embedded", "broadside", "scan30")

#: Scan angle of the ``scan30`` variant, degrees from broadside in the E-plane
SCAN30_THETA_DEG = 30.0


# ---------------------------------------------------------------------------
# Derived element design (Hammerstad / standard rectangular patch)
# ---------------------------------------------------------------------------

LAMBDA0_MM = C0 / F0_HZ * 1e3


def _w_for_free_space() -> float:
    """Radiating width W: half-wavelength in the quasi-static medium."""
    return LAMBDA0_MM / 2.0 * math.sqrt(2.0 / (SUBSTRATE_EPS_R + 1.0))


PATCH_W_MM = _w_for_free_space()          # along y


def patch_eps_eff(w_mm: float = PATCH_W_MM) -> float:
    """Effective permittivity (Hammerstad) for the given width."""
    h = SUBSTRATE_H_MM
    return ((SUBSTRATE_EPS_R + 1.0) / 2.0
            + (SUBSTRATE_EPS_R - 1.0) / 2.0 * (1.0 + 12.0 * h / w_mm) ** -0.5)


def _delta_l(w_mm: float = PATCH_W_MM) -> float:
    """Fringing length extension dL (mm)."""
    h = SUBSTRATE_H_MM
    eps_eff = patch_eps_eff(w_mm)
    wh = w_mm / h
    return (0.412 * h * (eps_eff + 0.3) * (wh + 0.264)
            / ((eps_eff - 0.258) * (wh + 0.8)))


def patch_length(eps_eff: float | None = None) -> float:
    """Resonant length L (mm): lambda_eff/2 minus two fringing extensions."""
    if eps_eff is None:
        eps_eff = patch_eps_eff()
    return LAMBDA0_MM / (2.0 * math.sqrt(eps_eff)) - 2.0 * _delta_l()


#: Analytic Hammerstad start (never used directly, kept for the report)
PATCH_L_HAMMERSTAD_MM = patch_length()
#: Retuned resonant length: the Hammerstad value 1.1383 mm resonated at
#: ~73 GHz and left Zin(77 GHz) = 66.6 + 10.8j ohm at the near-edge feed;
#: ``fit_model.py`` linearised the measured runs and ``make_candidate.py``
#: verified L = 1.1610 mm -> Zin(77) = 53.8 + 14.6j ohm, |S11| = -16.9 dB.
PATCH_L_MM = 1.1610

#: Element pitch = lambda0/2 (no grating lobes up to end-fire scan)
PITCH_MM = LAMBDA0_MM / 2.0

#: Probe inset from the -x patch edge.  The feed probe is an open
#: transmission-line section (Zp = 188 ohm, theta = 43 deg at 77 GHz, fitted
#: from the probe_only control run), so the classic 50 ohm inset at
#: 0.30 L transforms to 14.8 + 68j ohm - unusable.  The retune keeps the
#: probe almost at the radiating edge (0.05 L) and moves the resonance with
#: L instead (see fit_model.py / make_candidate.py).
FEED_INSET_MM = 0.05 * PATCH_L_MM

#: Port resistance
PORT_R_OHM = 50.0


def substrate_kappa() -> float:
    """Equivalent conductivity (S/m) for tan_d at f0.

    VisualopenEMS writes openEMS ``AddMaterial(epsilon=..., kappa=...)`` -
    openEMS has no loss-tangent argument, so the dielectric loss is converted
    at the design frequency (same convention as the GUI's preset builders).
    """
    return 2.0 * math.pi * F0_HZ * EPS0 * SUBSTRATE_EPS_R * SUBSTRATE_TAND


def retune_inset(r_measured_ohm: float, inset_old_mm: float) -> float:
    """Correct the probe inset from one measured input resistance.

    For an inset feed ``R_in(y0) = R_edge * cos^2(pi*y0/L)`` with ``y0``
    measured from the patch edge, so one measurement at ``inset_old`` gives
    ``R_edge`` and from it the inset that lands on 50 ohm::

        R_edge = R_meas / cos^2(pi*y0_old/L)
        y0_new = (L/pi) * arccos(sqrt(50/R_edge))

    Returns the new inset in mm (clamped to 5 % .. 45 % of L).
    """
    if not (0.0 < r_measured_ohm < 1e6):
        raise ValueError(f"implausible measured resistance {r_measured_ohm}")
    c_old = math.cos(math.pi * inset_old_mm / PATCH_L_MM)
    r_edge = r_measured_ohm / max(c_old * c_old, 1e-12)
    ratio = PORT_R_OHM / r_edge
    if ratio >= 1.0:
        # edge resistance already <= 50 ohm: keep the feed close to the edge
        y0 = 0.05 * PATCH_L_MM
    else:
        y0 = PATCH_L_MM / math.pi * math.acos(math.sqrt(ratio))
    return min(max(y0, 0.05 * PATCH_L_MM), 0.45 * PATCH_L_MM)


# ---------------------------------------------------------------------------
# Layout helpers
# ---------------------------------------------------------------------------

def element_xy(i: int, j: int) -> tuple[float, float]:
    """Centre of element (i, j) with i = x index, j = y index (0-based)."""
    x = (i - (N_X - 1) / 2.0) * PITCH_MM
    y = (j - (N_Y - 1) / 2.0) * PITCH_MM
    return x, y


def port_number(i: int, j: int) -> int:
    """1-based port number of element (i, j) - row-major, x fastest."""
    return j * N_X + i + 1


def port_xy(n: int) -> tuple[float, float]:
    """Probe position of port ``n`` (1-based)."""
    i = (n - 1) % N_X
    j = (n - 1) // N_X
    x, y = element_xy(i, j)
    # probe inset from the -x edge of the patch
    return x - PATCH_L_MM / 2.0 + FEED_INSET_MM, y


def patch_rect(i: int, j: int) -> tuple[float, float, float, float]:
    """(x_min, y_min, x_max, y_max) of element (i, j)."""
    xc, yc = element_xy(i, j)
    return (xc - PATCH_L_MM / 2.0, yc - PATCH_W_MM / 2.0,
            xc + PATCH_L_MM / 2.0, yc + PATCH_W_MM / 2.0)


def array_extent_mm() -> tuple[float, float]:
    """Bounding box of all patches: (x extent, y extent)."""
    return ((N_X - 1) * PITCH_MM + PATCH_L_MM,
            (N_Y - 1) * PITCH_MM + PATCH_W_MM)


#: Beam steering angles (deg, E-plane) of the full-wave sweep animation.
SWEEP_ANGLES_DEG = [-50, -40, -30, -20, -10, 0, 10, 20, 30, 40, 50]


def sweep_tag(theta_deg: float) -> str:
    """Short tag for a sweep angle: -50 -> 'm50', 0 -> '000', 30 -> 'p30'."""
    a = int(round(theta_deg))
    if a < 0:
        return f"m{abs(a):02d}"
    if a > 0:
        return f"p{a:02d}"
    return "000"


def scan_phases_deg(theta_deg: float, phi_deg: float = 0.0) -> dict[int, float]:
    """Excitation phase (deg) per port for a beam at (theta, phi).

    Pattern multiplication: the array far field is maximum where the element
    phase ``k0 * r_n . r_hat`` is compensated.  With ``r_n = (x_n, y_n, 0)``::

        phase_n = -360 deg * (x_n*sin(theta)*cos(phi) + y_n*sin(theta)*sin(phi)) / lambda0

    The **minus** sign follows openEMS's ``exp(+j*omega*t)`` phasor
    convention: a plane wave travelling in ``+x`` carries ``exp(-j*k*x)``, so
    the excitation of the element at ``x_n`` must lead by ``-k0*x_n``.  This
    was verified by simulation: with the earlier plus sign the "scan +30 deg"
    run put its maximum at ``phi = 180 deg`` (i.e. -x); the flipped sign puts
    it at ``phi = 0 deg``.  ``phi = 0`` and ``phi = 90`` stay broadside planes
    either way.
    """
    th = math.radians(theta_deg)
    ph = math.radians(phi_deg)
    phases: dict[int, float] = {}
    for j in range(N_Y):
        for i in range(N_X):
            x, y = element_xy(i, j)
            phase = (-360.0 / LAMBDA0_MM
                     * (x * math.sin(th) * math.cos(ph)
                        + y * math.sin(th) * math.sin(ph)))
            phases[port_number(i, j)] = phase
    return phases


# ---------------------------------------------------------------------------
# VisualopenEMS project (.fcsim.json)
# ---------------------------------------------------------------------------

def _mesh_dict() -> dict[str, Any]:
    """FDTD mesh for the 77 GHz structure (drawing unit mm).

    ``max_cellsize`` = lambda_air(80 GHz)/19 keeps every cell below the
    lambda/10 rule in air and below lambda/10 *inside* the substrate
    (2.52 mm / 10 = 0.25 mm); ``min_cellsize_z`` = h/4 resolves the 0.254 mm
    slab with four cells.  ``refined_cellsize`` is the target at metal edges.
    """
    return {
        "preset": "high",
        "drawing_unit": "mm",
        "start_frequency": F_START_HZ,
        "stop_frequency": F_STOP_HZ,
        "target_frequency": F0_HZ,
        "refined_cellsize": 0.10,
        "min_cellsize": 0.07,
        "max_cellsize": 0.25,
        "min_cellsize_z": 0.06,
        "cells_per_wavelength": 20.0,
        "max_grading_ratio": 1.6,
        "smooth_metal_edge": "one_third_two_thirds",
        "use_circle_detection": False,
        "use_arc_detection": False,
        "handle_closely_placed_edges": True,
        "num_lines": 5,
        "dirs": "xyz",
        "align_coordinates": True,
        "diagnose_cell_sources": True,
        "manual_mesh_lines": None,
        "manual_mesh_lines_removed": None,
        "mesh_repairs": None,
        "enforce_min_cellsize": False,
        "accepted_cells": None,
        "add_graded_mesh_lines_at_material_transitions": True,
        "coord_system": "cartesian",
        "excitation": "gauss",
        "excitation_file": "",
    }


def _port_dict(n: int, x: float, y: float, excite: bool, phase_deg: float,
               name: str | None = None) -> dict[str, Any]:
    return {
        "name": name or f"P{n}",
        "type": "lumped",
        "start": [round(x, 6), round(y, 6), 0.0],
        "stop": [round(x, 6), round(y, 6), SUBSTRATE_H_MM],
        "resistance": PORT_R_OHM,
        "inductance": 0.0,
        "capacitance": 0.0,
        "direction": "z",
        "excite_port": bool(excite),
        "excite_phase_deg": float(phase_deg),
        "active": True,
        "prop_direction": "",
        "exc_direction": "",
        "mesh_dirs": "xy",
        "metal_edge_res": None,
        "lumped_width": "line",
        "metal": "",
        "wg_a": 0.0,
        "wg_b": 0.0,
        "wg_mode": "TE10",
        "reverse": False,
    }


def _dump_dict(d: dict[str, str]) -> dict[str, Any]:
    return {
        "name": d["name"],
        "quantity": d["quantity"],
        "start": list(d["start"]),
        "stop": list(d["stop"]),
        "frequency_domain": True,
        "frequency": F0_HZ,
        "sub_sampling": [1, 1, 1],
        "enabled": True,
        "file_type": 1,
        "dump_mode": 1,
        "frequencies": None,
    }


def _primitives(variant: str) -> list[dict[str, Any]]:
    """Substrate + ground + the patches of the selected variant."""
    half = BOARD_MM / 2.0
    prims: list[dict[str, Any]] = [
        {
            "kind": "box",
            "name": "substrate",
            "material": SUBSTRATE_NAME,
            "priority": 10,
            "start": [-half, -half, 0.0],
            "stop": [half, half, SUBSTRATE_H_MM],
            "edges2grid": "",
        },
        {
            "kind": "box",
            "name": "ground",
            "material": "copper",
            "priority": 100,
            "start": [-half, -half, 0.0],
            "stop": [half, half, 0.0],
            "edges2grid": "xy",
        },
    ]

    if variant == "element":
        patches = [(0, 0)]
    else:
        patches = [(i, j) for j in range(N_Y) for i in range(N_X)]

    for i, j in patches:
        if variant == "element":
            # single patch in the centre of the same 12 x 12 mm board
            x0, y0 = -PATCH_L_MM / 2.0, -PATCH_W_MM / 2.0
            x1, y1 = PATCH_L_MM / 2.0, PATCH_W_MM / 2.0
        else:
            x0, y0, x1, y1 = patch_rect(i, j)
        tag = "patch" if variant == "element" else f"patch_{j * N_X + i + 1:02d}"
        prims.append({
            "kind": "box",
            "name": tag,
            "material": "copper",
            "priority": 100,
            "start": [round(x0, 6), round(y0, 6), SUBSTRATE_H_MM],
            "stop": [round(x1, 6), round(y1, 6), SUBSTRATE_H_MM],
            "edges2grid": "xy",
        })
    return prims


def _ports(variant: str, theta_deg: float | None = None) -> list[dict[str, Any]]:
    if variant == "element":
        # probe of the single, board-centred patch
        return [_port_dict(1, -PATCH_L_MM / 2.0 + FEED_INSET_MM, 0.0,
                           True, 0.0, name="P1")]

    if variant == "scan":
        if theta_deg is None:
            raise ValueError("variant 'scan' needs theta_deg")
        phases = scan_phases_deg(theta_deg)
    elif variant == "scan30":
        phases = scan_phases_deg(SCAN30_THETA_DEG)
    else:
        phases = {}

    excited_only = variant == "embedded"
    ports = []
    for j in range(N_Y):
        for i in range(N_X):
            n = port_number(i, j)
            x, y = port_xy(n)
            if variant in ("broadside", "scan30", "scan"):
                excite, phase = True, phases.get(n, 0.0)
            elif excited_only:
                # P6 = inner element (i=1, j=1): embedded element pattern and
                # the coupling column S*n6 of the S matrix in one run
                excite, phase = (n == 6), 0.0
            else:  # pragma: no cover - defensive
                excite, phase = False, 0.0
            ports.append(_port_dict(n, x, y, excite, round(phase, 6)))
    return ports


def _notes(variant: str, theta_deg: float | None = None) -> str:
    common = (
        f"4x4 patch array, 77 GHz automotive radar, RT/duroid 5880 "
        f"(h={SUBSTRATE_H_MM} mm, eps_r={SUBSTRATE_EPS_R}, "
        f"tan_d={SUBSTRATE_TAND}); patch {PATCH_L_MM:.3f} x {PATCH_W_MM:.3f} mm "
        f"(L along x), pitch {PITCH_MM:.3f} mm = lambda0/2, board "
        f"{BOARD_MM:g} x {BOARD_MM:g} mm, probe feed {FEED_INSET_MM:.3f} mm "
        f"inset, 50 ohm lumped ports."
    )
    extra = {
        "element": " Single patch at the board centre (isolated reference).",
        "embedded": (" 16 patches, only P6 excited - embedded element pattern "
                     "and coupling column S*n6."),
        "broadside": (" 16 ports excited in phase - true 4x4 broadside "
                      "pattern including mutual coupling."),
        "scan30": (f" 16 ports with progressive phase for a beam steered to "
                   f"theta=+{SCAN30_THETA_DEG:g} deg (E-plane)."),
        "scan": (f" 16 ports with progressive phase for a beam steered to "
                 f"theta={theta_deg:g} deg (E-plane)."),
    }[variant]
    return common + extra


def project_dict(name: str, variant: str,
                 theta_deg: float | None = None) -> dict[str, Any]:
    """Complete VisualopenEMS project dictionary for one excitation variant.

    ``variant="scan"`` additionally needs ``theta_deg`` (beam steering angle
    used for the progressive element phases).
    """
    allowed = VARIANTS + ("scan",)
    if variant not in allowed:
        raise ValueError(f"unknown variant {variant!r}, expected one of {allowed}")
    if variant == "scan" and theta_deg is None:
        raise ValueError("variant 'scan' needs theta_deg")

    half = BOARD_MM / 2.0
    return {
        "name": name,
        "version": 1,
        "geometry_files": [],
        "materials": [
            {
                "name": "copper",
                "kind": "pec",
                "epsilon": 1.0,
                "permeability": 1.0,
                "kappa": 0.0,
                "tand": 0.0,
                "color": [0.85, 0.55, 0.2, 1.0],
            },
            {
                "name": SUBSTRATE_NAME,
                "kind": "dielectric",
                "epsilon": SUBSTRATE_EPS_R,
                "permeability": 1.0,
                "kappa": round(substrate_kappa(), 12),
                "tand": SUBSTRATE_TAND,
                "color": [0.2, 0.5, 0.85, 0.6],
            },
        ],
        "assignments": [],
        "ports": _ports(variant, theta_deg),
        "dumps": [_dump_dict(d) for d in DUMPS],
        "nf2ff": {
            "enabled": True,
            "box_mode": "manual",
            "start": list(NF2FF_START_MM),
            "stop": list(NF2FF_STOP_MM),
            "center": list(NF2FF_CENTER_MM),
            "radius": 1.0,
            "frequencies": list(NF2FF_FREQS_HZ),
            "freq_mode": "list",
            "freq_count": len(NF2FF_FREQS_HZ),
            "theta_steps": NF2FF_THETA_STEPS,
            "phi_steps": NF2FF_PHI_STEPS,
        },
        "mesh": _mesh_dict(),
        "boundary": {
            "xmin": "PML",
            "xmax": "PML",
            "ymin": "PML",
            "ymax": "PML",
            "zmin": "PML",
            "zmax": "PML",
            "pml_size": 8,
            "air_margin": list(AIR_MARGIN_MM),
        },
        "sim": {
            "end_criteria": END_CRITERIA,
            "max_timesteps": MAX_TIMESTEPS,
            "max_time": None,
            "nf2ff_enabled": True,
            "field_dumps_enabled": True,
            "num_threads": 0,
            "output_dir": "",
            "antenna_mode": False,
            "antenna_pml_lambda": 0.5,
            "antenna_nf2ff_lambda": 0.5,
            "antenna_cells_per_wavelength": 15.0,
            "verbose": 2,
            "engine": "",
            "setup_only": False,
            "cleanup": True,
            "disable_dumps": False,
            "dump_statistics": False,
            "debug_material": False,
            "debug_pec": False,
            "debug_boxes": False,
            "debug_operator": False,
            "debug_csx": False,
            "show_probe_discretization": False,
            "legacy_hdf5_dumps": False,
            "over_sampling": None,
            "time_step_factor": None,
            "time_step_method": None,
            "cell_constant_material": False,
        },
        "mesh2freecad": {
            "export_step": False,
            "fcstd_object_type": "mesh",
            "voxel_grid": True,
            "voxel_grid_mode": "all",
            "include_sheets": True,
            "export_ports": True,
            "keep_temp": False,
        },
        "gerber": _gerber_block(half),
        "kicad": _kicad_block(half),
        "cst": {
            "source": "",
            "include_materials": True,
            "include_ports": True,
            "include_frequency": True,
            "include_boundaries": True,
            "include_sheets": True,
            "flip_y": False,
            "port_names": [],
        },
        "freecad_cmd": "",
        "notes": _notes(variant, theta_deg),
        "patch_params": None,
        "primitives": _primitives(variant),
        "object_display": {},
        "results": None,
    }


def _gerber_block(half: float) -> dict[str, Any]:
    return {
        "source": "",
        "board_thickness_mm": SUBSTRATE_H_MM,
        "copper_thickness_mm": COPPER_UM / 1000.0,
        "copper_as_sheet": True,
        "via_plating_radius_mm": 0.025,
        "via_max_diameter_mm": 0.8,
        "include_npth": False,
        "arc_error_mm": 0.01,
        "simplify_mm": 0.0,
        "min_area_mm2": 0.005,
        "clip_to_outline": True,
        "create_substrate": True,
        "substrate_material": SUBSTRATE_NAME,
        "substrate_epsilon": SUBSTRATE_EPS_R,
        "substrate_tand": SUBSTRATE_TAND,
        "include_mask": False,
        "solder_mask_thickness_mm": 0.02,
        "solder_mask_epsilon": 3.3,
        "solder_mask_tand": 0.02,
        "solder_mask_material": "solder_mask",
        "layer_z": {},
        "via_spans": {},
        "suggest_ports": True,
    }


def _kicad_block(half: float) -> dict[str, Any]:
    return {
        "source": "",
        "use_stackup": True,
        "board_thickness_mm": SUBSTRATE_H_MM,
        "copper_thickness_mm": COPPER_UM / 1000.0,
        "copper_as_sheet": True,
        "copper_mode": "",
        "substrate_material": SUBSTRATE_NAME,
        "substrate_epsilon": SUBSTRATE_EPS_R,
        "substrate_tand": SUBSTRATE_TAND,
        "via_plating_radius_mm": 0.025,
        "include_npth": False,
        "include_zones": True,
        "include_unfilled_zones": False,
        "include_pad_barrels": True,
        "include_tracks": True,
        "include_smd_pads": True,
        "include_vias": True,
        "flip_y": True,
        "suggest_ports": True,
        "arc_error_mm": 0.01,
        "simplify_mm": 0.0,
        "min_area_mm2": 0.005,
        "clip_to_outline": True,
        "create_substrate": True,
    }


# ---------------------------------------------------------------------------
# Human readable summary
# ---------------------------------------------------------------------------

def summary() -> str:
    """Short text block for the README / console."""
    ext_x, ext_y = array_extent_mm()
    eps_eff = patch_eps_eff()
    return (
        f"77 GHz automotive radar, 4x4 patch array on RT/duroid 5880\n"
        f"  substrate : h={SUBSTRATE_H_MM} mm (10 mil), eps_r={SUBSTRATE_EPS_R}, "
        f"tan_d={SUBSTRATE_TAND}, kappa={substrate_kappa():.4f} S/m\n"
        f"  lambda0   : {LAMBDA0_MM:.4f} mm  ->  pitch d = lambda0/2 = "
        f"{PITCH_MM:.4f} mm\n"
        f"  patch     : L = {PATCH_L_MM:.4f} mm (x, resonant), "
        f"W = {PATCH_W_MM:.4f} mm (y), eps_eff = {eps_eff:.4f}\n"
        f"  feed      : 50 ohm probe, inset {FEED_INSET_MM:.4f} mm from the "
        f"-x edge\n"
        f"  array     : {N_X}x{N_Y}, extent {ext_x:.2f} x {ext_y:.2f} mm on a "
        f"{BOARD_MM:g} x {BOARD_MM:g} mm board\n"
        f"  band      : {F_START_HZ/1e9:.0f}..{F_STOP_HZ/1e9:.0f} GHz "
        f"(gauss f0={F0_HZ/1e9:.0f} GHz, fc={FC_HZ/1e9:.1f} GHz)\n"
        f"  scan demo : theta = +{SCAN30_THETA_DEG:g} deg -> phase step "
        f"{360.0 * PITCH_MM * math.sin(math.radians(SCAN30_THETA_DEG)) / LAMBDA0_MM:.1f} deg "
        f"per element in x"
    )


if __name__ == "__main__":
    print(summary())
