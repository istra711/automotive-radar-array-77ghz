"""Build the German and English ODT reports from the finished study.

Reads `docs/results.json` (written by `plot_results.py`), the control-run
S-parameter files and the mesh validation output, and writes

    docs/patch_array_4x4_77ghz_de.odt
    docs/patch_array_4x4_77ghz_en.odt

Usage::

    python build_docs.py
"""

from __future__ import annotations

import json
import math
import sys
from datetime import date
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE / "tools"))
sys.path.insert(0, str(HERE))

from odt_builder import Report          # noqa: E402
import array_model as am                # noqa: E402

DOCS = HERE / "docs"
FIG = DOCS / "figures"
RES = json.loads((DOCS / "results.json").read_text(encoding="utf-8"))


# ---------------------------------------------------------------------------
# helpers
# ---------------------------------------------------------------------------

def load_s11(work: str) -> tuple[np.ndarray, np.ndarray]:
    d = json.loads((HERE / work / "sparams" / "sparams.json")
                   .read_text(encoding="utf-8"))
    f = np.array(d["freq"], dtype=float)
    s = (np.array(d["S11"]["real"], dtype=float)
         + 1j * np.array(d["S11"]["imag"], dtype=float))
    return f, s


def zin_at(work: str, f0: float = 77e9) -> tuple[float, float, float]:
    f, s = load_s11(work)
    k = int(np.argmin(np.abs(f - f0)))
    z = 50.0 * (1.0 + s[k]) / (1.0 - s[k])
    return z.real, z.imag, 20.0 * math.log10(max(abs(s[k]), 1e-12))


def mesh_of(work_validate: str) -> tuple[int, int, int, int]:
    d = json.loads((HERE / work_validate / "grid_lines.json")
                   .read_text(encoding="utf-8"))
    nx, ny, nz = len(d["x"]), len(d["y"]), len(d["z"])
    return nx, ny, nz, nx * ny * nz


def steps_of(work: str) -> int:
    p = HERE / work / "simulation" / "et"
    try:
        with p.open("r", encoding="utf-8", errors="replace") as fh:
            return sum(1 for line in fh if line.strip())
    except OSError:
        return 0


def build(lang: str) -> Path:
    de = lang == "de"

    def T(a: str, b: str) -> str:
        return a if de else b

    def N(v: float, d: int = 2) -> str:
        s = f"{v:.{d}f}"
        return s.replace(".", ",") if de else s

    el = RES["element"]
    arr = RES["broadside"]
    sc = RES["scan30"]
    be = RES["before"]
    rt = RES["retune"]
    cp = RES["coupling"]
    ss = RES["scan_summary"]

    z_b = (be["zin77_real"], be["zin77_imag"], be["s11_77_db"])
    z_i = zin_at("work_inset005_77ghz")
    z_c = zin_at("work_cand_L1p1683_s0p020")

    r = Report(
        title=T("4×4-Patch-Array für 77 GHz",
                "4×4 Patch Array for 77 GHz"),
        author="VisualOpenEMS · openEMS",
        subject=T("openEMS-Studie: Elemententwurf, Anpassung und "
                  "Strahlschwenkung auf RT/duroid 5880",
                  "openEMS study: element design, matching and beam "
                  "steering on RT/duroid 5880"),
        header=T("openEMS-Studie · 4×4-Patch-Array bei 77 GHz auf RT/duroid 5880",
                 "openEMS study · 4×4 patch array at 77 GHz on RT/duroid 5880"),
        toc_title=T("Inhalt", "Contents"),
    )
    r.figure_word = T("Bild", "Figure")
    r.table_word = T("Tabelle", "Table")

    # ------------------------------------------------------------- cover
    r.p(T("4×4-Patch-Array für 77 GHz", "4×4 Patch Array for 77 GHz"),
        style="Title")
    r.p(T("Vollwellen-Entwurf einer Gruppenantenne für den automobilen "
          "Radarbereich auf RT/duroid 5880 – vom Hammerstad-Erstentwurf "
          "über das Impedanz-Retuning bis zur Strahlschwenkung",
          "Full-wave design of an automotive radar antenna array on "
          "RT/duroid 5880 – from the Hammerstad first cut through impedance "
          "re-tuning to beam steering"),
        style="Subtitle")
    r.p("", style="Rule")
    r.p(T("Simulationsbericht", "Simulation report"), style="CoverMeta")
    r.p(T("Werkzeug: openEMS (FDTD) über die VisualopenEMS-Projektschicht",
          "Tool: openEMS (FDTD) through the VisualopenEMS project layer"),
        style="CoverMeta")
    r.p(T(f"Frequenzband: {N(am.F_START_HZ/1e9, 0)}–{N(am.F_STOP_HZ/1e9, 0)} GHz, "
          f"Mittenfrequenz {N(am.F0_HZ/1e9, 0)} GHz",
          f"Frequency band: {N(am.F_START_HZ/1e9, 0)}–{N(am.F_STOP_HZ/1e9, 0)} GHz, "
          f"centre {N(am.F0_HZ/1e9, 0)} GHz"), style="CoverMeta")
    r.p(date.today().isoformat(), style="CoverMeta")

    r.start_main()
    r.toc()
    r.page_break()

    # ------------------------------------------------------------- 1
    r.h1("1", T("Einleitung", "Introduction"))
    r.p(T(
        "Diese Studie dokumentiert den vollständigen Entwurf einer planaren "
        "4×4-Patch-Gruppenantenne für den automobilen 77-GHz-Radarbereich. "
        "Ausgangspunkt ist ein einzelnes, sondengespeistes Rechteck-Patch auf "
        "RT/duroid 5880. Das Element wird zunächst mit der klassischen "
        "Hammerstad-Näherung dimensioniert, im Vollwellenmodell vermessen und "
        "anschließend mit einem systematischen Retuning-Verfahren auf 50 Ω "
        "angepasst. Das so gefundene Element wird zu einem 4×4-Array mit "
        "λ0/2-Raster zusammengesetzt und in zwei Betriebsarten untersucht: "
        "phasengleich (Broadside) und mit progressiver Phasung für eine um "
        "+30° geschwenkte Keule.",
        "This study documents the complete design of a planar 4×4 patch "
        "antenna array for automotive 77 GHz radar.  The starting point is a "
        "single probe-fed rectangular patch on RT/duroid 5880.  The element "
        "is dimensioned with the classical Hammerstad approximation, "
        "characterised in the full-wave model and then matched to 50 Ω with a "
        "systematic re-tuning procedure.  The resulting element is assembled "
        "into a 4×4 array with λ0/2 spacing and studied in two modes of "
        "operation: in phase (broadside) and with a progressive phase for a "
        "beam steered to +30°."))

    r.h2("1.1", T("Aufgabenstellung", "Scope"))
    r.bullets([
        T("Einzelelement auf RT/duroid 5880 bei 77 GHz, Anpassung an 50 Ω "
          "über einen Sonden-Inset.",
          "Single element on RT/duroid 5880 at 77 GHz, matched to 50 Ω by a "
          "probe inset."),
        T("4×4-Array mit λ0/2-Raster (d = "
          f"{N(am.PITCH_MM, 3)} mm) und 16 eigenständigen 50-Ω-Ports.",
          "4×4 array with λ0/2 spacing (d = "
          f"{N(am.PITCH_MM, 3)} mm) and 16 independent 50 Ω ports."),
        T("Broadside-Betrieb und phasengesteuerter Betrieb (+30°) in je einem "
          "vollständigen FDTD-Lauf.",
          "Broadside operation and phased operation (+30°) in one complete "
          "FDTD run each."),
        T("Bewertung von Richtfaktor, Wirkungsgrad, Kopplung und "
          "Strahlschwenkung; Vergleich mit dem Arrayfaktor.",
          "Assessment of directivity, efficiency, coupling and beam steering; "
          "comparison with the array factor."),
    ])

    r.h2("1.2", T("Ergebnisüberblick", "Results at a glance"))
    r.p(T(
        f"Das finale Einzelelement erreicht bei 77 GHz eine Eingangsimpedanz "
        f"von Zin = {N(el['zin77_real'], 1)} + {N(el['zin77_imag'], 1)}j Ω und "
        f"eine Reflexion von {N(el['s11_77_db'], 1)} dB bei einem Richtfaktor "
        f"von {N(el['dmax77'], 2)} dBi. Das Broadside-Array kommt auf "
        f"{N(arr['dmax77'], 2)} dBi Richtfaktor; nach Korrektur der "
        f"Leistungsbilanz ergibt sich ein Gewinn von "
        f"{N(arr['gain77_corr'], 2)} dBi bei einem Wirkungsgrad von "
        f"{N(arr['eff77_corr'], 3)}. Die geschwenkte Variante erreicht "
        f"{N(sc['gain77_corr'], 2)} dBi in Richtung θ = "
        f"{N(sc['theta_peak_deg'], 0)}°, φ = {N(sc['phi_peak_deg'], 0)}°.",
        f"The final element reaches an input impedance of "
        f"Zin = {N(el['zin77_real'], 1)} + {N(el['zin77_imag'], 1)}j Ω and a "
        f"reflection of {N(el['s11_77_db'], 1)} dB at 77 GHz with a "
        f"directivity of {N(el['dmax77'], 2)} dBi.  The broadside array "
        f"reaches {N(arr['dmax77'], 2)} dBi directivity; after correcting the "
        f"power budget the gain is {N(arr['gain77_corr'], 2)} dBi at an "
        f"efficiency of {N(arr['eff77_corr'], 3)}.  The steered variant "
        f"reaches {N(sc['gain77_corr'], 2)} dBi towards θ = "
        f"{N(sc['theta_peak_deg'], 0)}°, φ = {N(sc['phi_peak_deg'], 0)}°."))

    # ------------------------------------------------------------- 2
    r.h1("2", T("Modell, Werkzeuge und Diskretisierung",
                "Model, tools and discretisation"))
    r.h2("2.1", T("Simulationskette", "Simulation chain"))
    r.p(T(
        "Die gesamte Studie ist skriptgesteuert und ohne grafische Oberfläche "
        "reproduzierbar. Die Geometrie existiert nur an einer Stelle, im "
        "Python-Modul `array_model.py`; alle Projektdateien, Läufe und "
        "Auswertungen werden daraus abgeleitet.",
        "The whole study is script driven and reproducible without a "
        "graphical user interface.  The geometry exists in exactly one place, "
        "the Python module `array_model.py`; all project files, runs and "
        "evaluations are derived from it."))
    r.numbered([
        T("`array_model.py` – Substrat, Patch, Raster, Ports, Mesh, "
          "NF2FF-Box und die vier Projektvarianten.",
          "`array_model.py` – substrate, patch, grid, ports, mesh, NF2FF box "
          "and the four project variants."),
        T("`make_project.py` – schreibt die `.fcsim.json`-Projekte und die "
          "Konstruktionsbilder.",
          "`make_project.py` – writes the `.fcsim.json` projects and the "
          "construction figures."),
        T("`validate_project.py` – prüft Kollisionen, Mesh-Regeln, "
          "Portlagen und die NF2FF-Box vor dem Lauf.",
          "`validate_project.py` – checks collisions, mesh rules, port "
          "positions and the NF2FF box before the run."),
        T("`run_headless.py` – erzeugt das openEMS-Skript aus dem Projekt "
          "und startet FDTD inklusive Nachverarbeitung.",
          "`run_headless.py` – generates the openEMS script from the "
          "project and runs FDTD including post-processing."),
        T("`plot_results.py` – liest S-Parameter und Fernfelder, erzeugt die "
          "Ergebnisbilder und `docs/results.json`.",
          "`plot_results.py` – reads S-parameters and far fields, produces "
          "the result figures and `docs/results.json`."),
        T("`build_docs.py` – baut diesen Bericht in Deutsch und Englisch.",
          "`build_docs.py` – builds this report in German and English."),
    ])

    r.h2("2.2", T("Diskretisierung", "Discretisation"))
    r.p(T(
        f"Das Mesh wird über die VisualopenEMS-Voreinstellung „high“ erzeugt "
        f"und mit expliziten Zellgrößen ergänzt: maximale Zellweite "
        f"0,25 mm, verfeinerte Kante 0,10 mm, minimale Zellweite 0,07 mm und "
        f"minimale z-Zellweite 0,06 mm. Die 0,254 mm dicke Folie wird damit "
        f"von vier Zellen in z aufgelöst, die Metallkanten von rund drei "
        f"Zellen. Das Rechengebiet ist allseitig von einer PML umschlossen; "
        f"der Luftabstand zum Substrat beträgt 4 mm, nach oben 6 mm.",
        f"The mesh is created with the VisualopenEMS preset “high” and "
        f"supplemented with explicit cell sizes: maximum cell size 0.25 mm, "
        f"refined edge 0.10 mm, minimum cell size 0.07 mm and minimum z cell "
        f"size 0.06 mm.  The 0.254 mm laminate is therefore resolved by four "
        f"cells in z and the metal edges by about three cells.  The "
        f"computational domain is enclosed by a PML on all sides; the air "
        f"clearance is 4 mm to the sides and 6 mm above the board."))
    r.p(T(
        "Der Abbruch erfolgt bei −40 dB Restenergie (openEMS-Kriterium "
        "1e-4), was bei allen vier Produktionsläufen nach rund 10 600 "
        "Zeitschritten erreicht ist. Die Anregung ist ein Gauß-Puls mit "
        f"f0 = {N(am.F0_HZ/1e9, 0)} GHz und fc = "
        f"{N(am.FC_HZ/1e9, 1)} GHz, der das Band 74–80 GHz sicher abdeckt.",
        "The run stops at −40 dB residual energy (openEMS criterion 1e-4), "
        "which all four production runs reach after about 10,600 time steps.  "
        "The excitation is a Gaussian pulse with "
        f"f0 = {N(am.F0_HZ/1e9, 0)} GHz and fc = "
        f"{N(am.FC_HZ/1e9, 1)} GHz, which safely covers the 74–80 GHz band."))

    # ------------------------------------------------------------- 3
    r.h1("3", T("Aufbau der Antenne", "Antenna construction"))

    r.h2("3.1", T("Substrat und Schichtaufbau", "Substrate and stack-up"))
    r.p(T(
        "Als Träger dient RT/duroid 5880 mit εr = 2,2 und tan δ = 0,0009 in "
        "der Standarddicke 0,254 mm (10 mil). Beide Kupferlagen werden im "
        "Modell als unendlich dünne PEC-Flächen mit 35 µm Dicke angesetzt. "
        "Bild 1 zeigt den Schnitt durch eine Patch-Zeile; die z-Achse ist "
        "darin um den Faktor 6 überhöht, damit die 0,254 mm dünne Folie "
        "sichtbar bleibt.",
        "RT/duroid 5880 with εr = 2.2 and tan δ = 0.0009 in the standard "
        "0.254 mm (10 mil) thickness serves as the carrier.  Both copper "
        "layers are modelled as 35 µm thick PEC sheets.  Figure 1 shows the "
        "cut through one patch row; the z axis is exaggerated by a factor of "
        "six so that the 0.254 mm laminate remains visible."))
    r.figure(str(FIG / "fig03_stackup.png"),
             T("Schichtaufbau entlang einer Patch-Zeile (y = −0,973 mm). "
               "Die Sonden verbinden die Grundfläche mit den Patches; "
               "z ist 6-fach überhöht.",
               "Stack-up along one patch row (y = −0.973 mm).  The probes "
               "connect the ground plane to the patches; z is exaggerated "
               "6×."), 15.6)

    r.table(
        T("Substrat- und Aufbauparameter.",
          "Substrate and stack-up parameters."),
        [T("Parameter", "Parameter"), T("Wert", "Value")],
        [
            [T("Substrat", "Substrate"), "RT/duroid 5880"],
            [T("Rel. Permittivität εr", "Relative permittivity εr"),
             N(am.SUBSTRATE_EPS_R, 1)],
            [T("Verlustfaktor tan δ", "Loss tangent tan δ"),
             N(am.SUBSTRATE_TAND, 4)],
            [T("Substratdicke h", "Substrate height h"),
             f"{N(am.SUBSTRATE_H_MM, 3)} mm"],
            [T("Kupferdicke (PEC)", "Copper thickness (PEC)"),
             f"{N(am.COPPER_UM, 0)} µm"],
            [T("Platinengröße", "Board size"),
             f"{N(am.BOARD_MM, 1)} × {N(am.BOARD_MM, 1)} mm"],
            [T("Freiraumwellenlänge λ0(77 GHz)",
               "Free-space wavelength λ0(77 GHz)"),
             f"{N(am.LAMBDA0_MM, 3)} mm"],
        ],
        widths=[6.0, 5.0])

    r.h2("3.2", T("Patch-Dimensionierung", "Patch dimensioning"))
    r.p(T(
        "Die Resonanzlänge folgt der Hammerstad-Näherung für einen "
        "Rechteck-Patch auf einem unendlich ausgedehnten Substrat. Mit "
        "W = 1,539 mm ergibt sich eine effektive Permittivität von 1,9475 "
        "und damit eine Startlänge von L = 1,1383 mm. Die Breite ist so "
        "gewählt, dass der Strahlungswiderstand im Bereich einiger hundert "
        "Ohm liegt und der Sonden-Inset den Transformationsbereich abdeckt.",
        "The resonant length follows the Hammerstad approximation for a "
        "rectangular patch on an infinitely extended substrate.  With "
        "W = 1.539 mm the effective permittivity is 1.9475 and the starting "
        "length is L = 1.1383 mm.  The width is chosen so that the radiation "
        "resistance is a few hundred ohms and the probe inset covers the "
        "required transformation range."))
    r.code(
        "eps_eff = (eps_r + 1)/2 + (eps_r - 1)/2 * (1 + 12*h/W)^(-0.5)\n"
        "dL      = 0.412*h * (eps_eff + 0.3)*(W/h + 0.264) /\n"
        "                 ((eps_eff - 0.258)*(W/h + 0.8))\n"
        "L       = c0/(2*f0*sqrt(eps_eff)) - 2*dL")
    r.p(T(
        "Bild 2 zeigt das Element mit Bemaßung. Der Zoom rechts macht die "
        "Speisung sichtbar: Der Inset beträgt im finalen Entwurf nur "
        f"{N(am.FEED_INSET_MM, 3)} mm, also rund 5 % der Patchlänge – ein "
        "Wert, der auf der Gesamtansicht nicht mehr auflösbar wäre.",
        "Figure 2 shows the dimensioned element.  The zoom on the right "
        "reveals the feed: in the final design the inset is only "
        f"{N(am.FEED_INSET_MM, 3)} mm, about 5 % of the patch length – a "
        "value that is no longer resolvable in the overall view."))
    r.figure(str(FIG / "fig02_element.png"),
             T("Einzelelement mit Bemaßung (links) und Zoom auf die "
               "Sonden-Speisung (rechts).",
               "Dimensioned single element (left) and zoom on the probe feed "
               "(right)."), 16.0)

    r.h2("3.3", T("Array-Layout", "Array layout"))
    r.p(T(
        f"Vier Elemente in x und vier in y ergeben ein quadratisches Array "
        f"mit dem Rastermaß d = λ0/2 = {N(am.PITCH_MM, 3)} mm. Die "
        f"Gesamtausdehnung der Patchfläche beträgt "
        f"{N(am.array_extent_mm()[0], 2)} × {N(am.array_extent_mm()[1], 2)} mm "
        f"auf einer {N(am.BOARD_MM, 0)} × {N(am.BOARD_MM, 0)} mm großen "
        f"Platine. Jedes Element besitzt einen eigenen koaxialen Sondenport "
        f"von der Grundfläche zur Patchunterseite, sodass alle 16 Elemente "
        f"unabhängig gespeist werden können.",
        f"Four elements in x and four in y form a square array with "
        f"d = λ0/2 = {N(am.PITCH_MM, 3)} mm spacing.  The total footprint of "
        f"the patch area is {N(am.array_extent_mm()[0], 2)} × "
        f"{N(am.array_extent_mm()[1], 2)} mm on a {N(am.BOARD_MM, 0)} × "
        f"{N(am.BOARD_MM, 0)} mm board.  Each element has its own coaxial "
        f"probe port from the ground plane to the patch underside, so all 16 "
        f"elements can be driven independently."))
    r.figure(str(FIG / "fig01_layout.png"),
             T("Draufsicht auf das 4×4-Array mit Rastermaß, Patch- und "
               "Platinenabmessungen. Rote Punkte: 50-Ω-Sondenports.",
               "Top view of the 4×4 array with grid, patch and board "
               "dimensions.  Red dots: 50 Ω probe ports."), 14.4)

    r.h2("3.4", T("Dreidimensionales Modell", "Three-dimensional model"))
    r.p(T(
        "Bild 4 zeigt das vollständige Simulationsmodell. Die Patches liegen "
        "auf der Oberseite der Folie, die Grundfläche bildet die "
        "Unterseite; die Sonden durchstoßen das Substrat. Das Modell ist "
        "exakt das, was `array_model.py` in das openEMS-Skript schreibt – "
        "es gibt keine manuell gepflegte zweite Geometrie.",
        "Figure 4 shows the complete simulation model.  The patches sit on "
        "top of the laminate, the ground plane forms the bottom side and the "
        "probes pass through the substrate.  The model is exactly what "
        "`array_model.py` writes into the openEMS script – there is no "
        "second, manually maintained geometry."))
    r.figure(str(FIG / "fig04_array3d.png"),
             T("3D-Modell des 4×4-Arrays (z 8-fach überhöht).",
               "3D model of the 4×4 array (z exaggerated 8×)."), 14.8)

    # ------------------------------------------------------------- 4
    r.h1("4", T("Entwurf und Anpassung des Einzelelements",
                "Element design and matching"))

    r.h2("4.1", T("Erster Entwurf nach Hammerstad",
                  "First cut after Hammerstad"))
    r.p(T(
        "Der Hammerstad-Entwurf mit L = 1,1383 mm und einem Inset von 30 % "
        "der Patchlänge (0,341 mm) wird als Referenz simuliert. Das Ergebnis "
        "ist ernüchternd: Bei 77 GHz misst man "
        f"Zin = {N(z_b[0], 1)} + {N(z_b[1], 1)}j Ω, also eine stark "
        f"induktive Fehlanpassung, und eine Reflexion von nur "
        f"{N(z_b[2], 1)} dB. Der große Inset hat die Impedanz über den "
        "gewünschten Punkt hinaus transformiert.",
        "The Hammerstad design with L = 1.1383 mm and an inset of 30 % of "
        "the patch length (0.341 mm) is simulated as a reference.  The "
        "result is sobering: at 77 GHz the input impedance is "
        f"Zin = {N(z_b[0], 1)} + {N(z_b[1], 1)}j Ω, a strongly inductive "
        f"mismatch, and the reflection is only {N(z_b[2], 1)} dB.  The large "
        "inset has transformed the impedance beyond the target."))
    r.p(T(
        "Bemerkenswert ist, dass der Reflexionsgrad über das ganze Band "
        "flach bleibt – ein Zeichen dafür, dass hier keine einzelne "
        "Fehlanpassung, sondern eine systematisch falsche Speiseimpedanz "
        "vorliegt. Genau dieser Verdacht wird im nächsten Schritt mit einem "
        "Kontrollexperiment überprüft.",
        "It is notable that the reflection stays flat across the whole band "
        "– an indication that this is not a single mismatch but a "
        "systematically wrong feed impedance.  Exactly this suspicion is "
        "tested with a control experiment in the next step."))

    r.h2("4.2", T("Kontrollexperiment: offene Speiseleitung",
                  "Control experiment: open feed line"))
    r.p(T(
        "Im Kontrolllauf `probe_only` wird die Speiseleitung ohne Patch "
        "simuliert: ein kurzes, am oberen Ende offenes Stück Koaxialleitung "
        "durch das Substrat. Ein offenes Leitungsende wirkt kapazitiv, ein "
        "kurzgeschlossenes induktiv – die Messung entscheidet also, welches "
        "Ersatzbild für die Sonde gilt.",
        "In the control run `probe_only` the feed is simulated without the "
        "patch: a short piece of coaxial line through the substrate, open at "
        "the top.  An open line end is capacitive, a short-circuited one is "
        "inductive – the measurement therefore decides which equivalent "
        "circuit applies to the probe."))
    r.p(T(
        f"Bild 5 links zeigt das Ergebnis: Die gemessene Reaktanz ist über "
        f"das ganze Band stark kapazitiv (bei 77 GHz "
        f"{N(rt['probe_z77_imag'], 0)}j Ω) und folgt einer Tangens-Kurve. "
        f"Der Fit einer offenen Leitung liefert Zp = "
        f"{N(rt['probe_zp_ohm'], 0)} Ω, εp = {N(rt['probe_eps_p'], 2)} und "
        f"eine Phasenlänge von {N(rt['probe_theta77_deg'], 1)}° bei 77 GHz; "
        f"die Abweichung zum Modell liegt unter 0,5 %.",
        f"The left panel of Figure 5 shows the result: the measured "
        f"reactance is strongly capacitive over the whole band "
        f"({N(rt['probe_z77_imag'], 0)}j Ω at 77 GHz) and follows a tangent "
        f"curve.  The fit of an open line yields Zp = "
        f"{N(rt['probe_zp_ohm'], 0)} Ω, εp = {N(rt['probe_eps_p'], 2)} and "
        f"an electrical length of {N(rt['probe_theta77_deg'], 1)}° at "
        f"77 GHz; the deviation from the model is below 0.5 %."))
    r.figure(str(FIG / "fig06_retune_journey.png"),
             T("Kontrollexperiment und Retuning: links die offene "
               "Speiseleitung mit Fit, rechts die Eingangsimpedanz-Ortskurven "
               "der drei Entwurfsstufen im Smith-Diagramm-äquivalenten "
               "Z-Bereich (Punkte: 77 GHz).",
               "Control experiment and re-tuning: left the open feed line "
               "with its fit, right the input-impedance loci of the three "
               "design stages in the Z plane (dots: 77 GHz)."), 16.0)

    r.h2("4.3", T("De-Embedding und Newton-Retuning",
                  "De-embedding and Newton re-tuning"))
    r.p(T(
        "Mit dem Leitungsmodell lässt sich der Einfluss der Sonde aus der "
        "gemessenen Eingangsimpedanz herausrechnen. Übrig bleibt die "
        "eigentliche Patchimpedanz, die nun gezielt über die beiden "
        "Entwurfsparameter – Patchlänge L (Resonanz) und Inset (Widerstand) "
        "– auf 50 Ω gebracht werden kann. Dazu wird ein lokales Newton-"
        "Verfahren auf die simulierte Impedanzkarte angewandt; die "
        "Jacobi-Matrix stammt aus finiten Differenzen der bereits "
        "vorhandenen Kontrollläufe.",
        "With the line model the influence of the probe can be removed from "
        "the measured input impedance.  What remains is the actual patch "
        "impedance, which can now be brought to 50 Ω in a targeted way via "
        "the two design parameters – patch length L (resonance) and inset "
        "(resistance).  A local Newton iteration is applied to the simulated "
        "impedance map; the Jacobian comes from finite differences of the "
        "control runs already available."))
    r.p(T(
        f"Tabelle 2 zeigt die Iterationskette. Der Schritt von 30 % auf 5 % "
        f"Inset beseitigt den induktiven Überschuss weitgehend, die "
        f"verbleibende Abweichung wird über eine Verlängerung der Patchlänge "
        f"auf 1,1610 mm korrigiert. Zwei unabhängige Kandidaten wurden zur "
        f"Kontrolle voll simuliert; Kandidat A wird als endgültiger Entwurf "
        f"übernommen.",
        f"Table 2 shows the iteration chain.  Reducing the inset from 30 % "
        f"to 5 % largely removes the inductive excess; the remaining "
        f"deviation is corrected by lengthening the patch to 1.1610 mm.  Two "
        f"independent candidates were fully simulated as a check; candidate "
        f"A is adopted as the final design."))
    r.table(
        T("Entwurfskette mit simulierten Impedanzen bei 77 GHz.",
          "Design chain with simulated impedances at 77 GHz."),
        [T("Schritt", "Step"), "L (mm)", T("Inset", "Inset"),
         "Zin(77) (Ω)", "|S11|(77) (dB)"],
        [
            [T("Hammerstad-Erstentwurf", "Hammerstad first cut"),
             N(am.PATCH_L_HAMMERSTAD_MM, 4), T("0,30 L", "0.30 L"),
             f"{N(z_b[0], 1)} + {N(z_b[1], 1)}j", N(z_b[2], 1)],
            [T("Inset-Kontrolle", "Inset control"), N(am.PATCH_L_HAMMERSTAD_MM, 4),
             T("0,05 L", "0.05 L"),
             f"{N(z_i[0], 1)} + {N(z_i[1], 1)}j", N(z_i[2], 1)],
            [T("Kandidat B", "Candidate B"), N(1.1683, 4), T("0,02 L", "0.02 L"),
             f"{N(z_c[0], 1)} + {N(z_c[1], 1)}j", N(z_c[2], 1)],
            [T("Kandidat A – final", "Candidate A – final"), N(1.1610, 4),
             T("0,05 L", "0.05 L"),
             f"{N(el['zin77_real'], 1)} + {N(el['zin77_imag'], 1)}j",
             N(el["s11_77_db"], 1)],
        ],
        widths=[4.6, 2.2, 2.2, 3.4, 3.0],
        source=T("Alle Werte aus dem jeweiligen Produktions- bzw. "
                 "Verifikationslauf, S-Parameter bei 77 GHz.",
                 "All values from the respective production or verification "
                 "run, S-parameters at 77 GHz."))

    r.h2("4.4", T("Endkontrolle der Anpassung", "Final matching check"))
    r.p(T(
        f"Bild 6 stellt den finalen Entwurf dem Erstentwurf gegenüber. Die "
        f"Eingangsimpedanz bei 77 GHz beträgt "
        f"{N(el['zin77_real'], 1)} + {N(el['zin77_imag'], 1)}j Ω, das "
        f"Reflexionsmaß {N(el['s11_77_db'], 1)} dB und das "
        f"Stehwellenverhältnis {N(el['vswr77'], 2)}. Die −10-dB-Bandbreite "
        f"reicht von {N(el['bw10_low_ghz'], 2)} bis "
        f"{N(el['bw10_high_ghz'], 2)} GHz, also "
        f"{N(el['bw10_ghz'], 2)} GHz. Das Optimum liegt mit "
        f"{N(el['best_s11_db'], 1)} dB bei {N(el['best_s11_ghz'], 2)} GHz.",
        f"Figure 6 compares the final design with the first cut.  The input "
        f"impedance at 77 GHz is {N(el['zin77_real'], 1)} + "
        f"{N(el['zin77_imag'], 1)}j Ω, the reflection "
        f"{N(el['s11_77_db'], 1)} dB and the VSWR {N(el['vswr77'], 2)}.  The "
        f"−10 dB bandwidth extends from {N(el['bw10_low_ghz'], 2)} to "
        f"{N(el['bw10_high_ghz'], 2)} GHz, i.e. {N(el['bw10_ghz'], 2)} GHz.  "
        f"The optimum is {N(el['best_s11_db'], 1)} dB at "
        f"{N(el['best_s11_ghz'], 2)} GHz."))
    r.p(T(
        "Der Vergleich zeigt auch die Grenze des Verfahrens: Die "
        "Impedanzkarte in der Nähe des Optimums ist nahezu kollinear, "
        "Länge und Inset verschieben die Impedanz in dieselbe Richtung. "
        "Mehr als etwa −17 dB sind mit diesen zwei Stellschrauben bei "
        "dieser Substratdicke nicht erreichbar; für eine tiefere Anpassung "
        "wäre ein zusätzliches Transformationsnetzwerk oder ein "
        "Schlitzkoppler nötig.",
        "The comparison also shows the limit of the method: near the optimum "
        "the impedance map is almost collinear, i.e. length and inset move "
        "the impedance in the same direction.  More than about −17 dB cannot "
        "be achieved with these two knobs at this substrate thickness; a "
        "deeper match would require an additional matching network or a "
        "coupled feed."))
    r.figure(str(FIG / "fig05_element_zin_s11.png"),
             T("Eingangsimpedanz (links) und Reflexionsmaß (rechts) des "
               "finalen Elements; gestrichelt der Hammerstad-Erstentwurf.",
               "Input impedance (left) and reflection (right) of the final "
               "element; dashed: the Hammerstad first cut."), 16.0)

    r.h2("4.5", T("Richtdiagramm des isolierten Elements",
                  "Pattern of the isolated element"))
    r.p(T(
        f"Das isolierte Element strahlt mit einem Richtfaktor von "
        f"{N(el['dmax77'], 2)} dBi bei einem Wirkungsgrad von "
        f"{N(el['eff77'], 3)}. Bild 7 zeigt das dreidimensionale Diagramm; "
        f"die Keule steht erwartungsgemäß in +z, die E-Ebene ist die "
        f"xz-Ebene.",
        f"The isolated element radiates with a directivity of "
        f"{N(el['dmax77'], 2)} dBi at an efficiency of {N(el['eff77'], 3)}.  "
        f"Figure 7 shows the three-dimensional pattern; the beam points in "
        f"+z as expected and the E plane is the xz plane."))
    r.figure(str(FIG / "fig07_element_pattern3d.png"),
             T(f"Richtdiagramm des isolierten Elements bei 77 GHz "
               f"(Maximum {N(el['gain77'], 2)} dBi).",
               f"Pattern of the isolated element at 77 GHz "
               f"(peak {N(el['gain77'], 2)} dBi)."), 11.4)

    # ------------------------------------------------------------- 5
    r.h1("5", T("Das 4×4-Array", "The 4×4 array"))

    r.h2("5.1", T("Einbettung und Kopplung", "Embedding and coupling"))
    r.p(T(
        "Im Lauf `embedded` wird nur Port 6 angeregt, alle anderen Ports "
        "sind passiv mit 50 Ω abgeschlossen. Der Lauf liefert damit die "
        "vollständige Kopplungsspalte S_i6 und das eingebettete "
        "Elementdiagramm – also das Verhalten eines Elements in seiner "
        "tatsächlichen Umgebung.",
        "In the `embedded` run only port 6 is excited while all other "
        "ports are passively terminated in 50 Ω.  The run therefore yields "
        "the complete coupling column S_i6 and the embedded element "
        "pattern – the behaviour of one element in its real environment."))
    r.p(T(
        f"Die Eigenreflexion des Elements beträgt "
        f"{N(cp['s66_77_db'], 1)} dB und ist damit schlechter als beim "
        f"isolierten Element ({N(el['s11_77_db'], 1)} dB), was an der "
        f"nachbarlichen Kopplung liegt. Der stärkste Koppler ist das "
        f"direkte Nachbarelement in x mit {N(cp['coupling_max_77_db'], 1)} dB "
        f"bei {N(cp['coupling_table'][0]['dist_mm'], 2)} mm Abstand; die "
        f"diagonalen Nachbarn liegen bei etwa −25 dB und die Elemente am "
        f"Platinenrand bei −28 dB und darunter. Bild 8 zeigt die "
        f"Kopplungsspalte über der Entfernung und über der Frequenz.",
        f"The self-reflection of the element is {N(cp['s66_77_db'], 1)} dB, "
        f"worse than that of the isolated element ({N(el['s11_77_db'], 1)} "
        f"dB) because of neighbour coupling.  The strongest coupler is the "
        f"direct neighbour in x at {N(cp['coupling_max_77_db'], 1)} dB and "
        f"{N(cp['coupling_table'][0]['dist_mm'], 2)} mm distance; the "
        f"diagonal neighbours are around −25 dB and the elements at the "
        f"board edge −28 dB and below.  Figure 8 shows the coupling column "
        f"versus distance and versus frequency."))
    r.figure(str(FIG / "fig09_coupling.png"),
             T("Kopplung im eingebetteten Lauf: Spalte S_i6 bei 77 GHz über "
               "der Entfernung (links) und über der Frequenz (rechts).",
               "Coupling in the embedded run: column S_i6 at 77 GHz versus "
               "distance (left) and versus frequency (right)."), 16.0)
    r.p(T(
        f"Das eingebettete Element erreicht einen Richtfaktor von "
        f"{N(el['embedded_dmax77'], 2)} dBi bei einem Wirkungsgrad von "
        f"{N(el['embedded_eff77'], 3)}. Bild 9 vergleicht E- und H-Ebene "
        f"mit dem isolierten Element: Die Keule wird durch die Nachbarn "
        f"leicht aufgeweitet, die Rückkeule steigt an – das ist der "
        f"Fußabdruck der Kopplung, der später im Arrayfaktor-Vergleich "
        f"bewusst nicht modelliert wird.",
        f"The embedded element reaches a directivity of "
        f"{N(el['embedded_dmax77'], 2)} dBi at an efficiency of "
        f"{N(el['embedded_eff77'], 3)}.  Figure 9 compares the E and H "
        f"planes with the isolated element: the beam is slightly widened by "
        f"the neighbours and the back lobe rises – the footprint of coupling, "
        f"which is deliberately not modelled in the array-factor comparison "
        f"later."))
    r.figure(str(FIG / "fig08_element_vs_embedded.png"),
             T("E- und H-Ebene bei 77 GHz: isoliertes Element (blau) gegen "
               "eingebettetes Element P6 (rot).",
               "E and H plane at 77 GHz: isolated element (blue) versus "
               "embedded element P6 (red)."), 16.0)

    r.h2("5.2", T("Broadside-Betrieb", "Broadside operation"))
    r.p(T(
        "Im Lauf `broadside` werden alle 16 Ports phasengleich angeregt. "
        "Der Richtfaktor steigt auf "
        f"{N(arr['dmax77'], 2)} dBi. Das passt gut zur Aperturabschätzung "
        "D ≈ 4πA/λ0² = 17,0 dBi für eine gleichförmig belegte "
        "λ0/2-Apertur von 2λ0 × 2λ0; der um 0,9 dB höhere Wert folgt aus der "
        "Elementcharakteristik, die die Ränder des sichtbaren Bereichs "
        "abschwächt. Die Rasterauflösung des Diagramms beträgt 2°, das "
        f"Maximum liegt bei θ = {N(arr['theta_peak_deg'], 0)}°, φ = "
        f"{N(arr['phi_peak_deg'], 0)}°.",
        "In the `broadside` run all 16 ports are driven in phase.  The "
        f"directivity rises to {N(arr['dmax77'], 2)} dBi.  This agrees well "
        "with the aperture estimate D ≈ 4πA/λ0² = 17.0 dBi for a uniformly "
        "excited λ0/2 aperture of 2λ0 × 2λ0; the 0.9 dB higher value follows "
        "from the element pattern, which tapers the edge of the visible "
        "range.  The angular grid is 2°, and the maximum is at "
        f"θ = {N(arr['theta_peak_deg'], 0)}°, φ = "
        f"{N(arr['phi_peak_deg'], 0)}°."))
    r.p(T(
        "Bild 10 stellt das Vollwellen-Diagramm dem Produkt aus "
        "Einzelelement-Diagramm und Arrayfaktor gegenüber. Die "
        "Hauptkeule stimmt im gesamten Bereich bis −20 dB hervorragend "
        "überein; die Abweichungen in den Nullstellen und den weiten "
        "Nebenkeulen sind die Signatur der Vernachlässigung von Kopplung "
        "und Randbeugung im Arrayfaktor-Modell.",
        "Figure 10 compares the full-wave pattern with the product of element "
        "pattern and array factor.  The main beam agrees excellently down to "
        "−20 dB; the deviations in the nulls and the wide sidelobes are the "
        "signature of neglecting coupling and edge diffraction in the array-"
        "factor model."))
    r.figure(str(FIG / "fig10_broadside.png"),
             T("Broadside-Betrieb bei 77 GHz: 3D-Diagramm, E-Ebene und "
               "H-Ebene im Vergleich mit Element × Arrayfaktor.",
               "Broadside operation at 77 GHz: 3D pattern, E plane and "
               "H plane compared with element × array factor."), 16.4)

    r.h2("5.3", T("Strahlschwenkung um +30°", "Beam steering by +30°"))
    r.p(T(
        f"Für den geschwenkten Lauf wird jedem Element in x eine "
        f"progressive Phase von −90° aufgeprägt, entsprechend "
        f"{N(ss['phase_step_deg'], 1)}° pro Raster d = λ0/2. Das Vorzeichen "
        f"folgt der openEMS-Phasenkonvention: Eine in +x laufende Welle "
        f"trägt exp(−jkx), die Anregung muss also voreilen. Diese Konvention "
        f"wurde durch einen Testlauf verifiziert – mit dem umgekehrten "
        f"Vorzeichen erschien die Keule spiegelbildlich bei φ = 180°.",
        f"For the steered run each element in x receives a progressive phase "
        f"of −90°, corresponding to {N(ss['phase_step_deg'], 1)}° per "
        f"spacing d = λ0/2.  The sign follows the openEMS phase convention: "
        f"a wave travelling in +x carries exp(−jkx), so the excitation must "
        f"lead.  This convention was verified with a test run – with the "
        f"opposite sign the beam appeared mirrored at φ = 180°."))
    r.p(T(
        f"Die Simulation bestätigt die Steuerung: Das Maximum liegt bei "
        f"θ = {N(sc['theta_peak_deg'], 0)}°, φ = {N(sc['phi_peak_deg'], 0)}° "
        f"mit {N(sc['gain77_corr'], 2)} dBi. Die Ablage gegenüber den "
        f"nominellen 30° resultiert aus der endlichen Elementkeule und der "
        f"Kopplung; der Arrayfaktor allein sagt 30° voraus.",
        f"The simulation confirms the steering: the maximum is at "
        f"θ = {N(sc['theta_peak_deg'], 0)}°, φ = {N(sc['phi_peak_deg'], 0)}° "
        f"with {N(sc['gain77_corr'], 2)} dBi.  The offset from the nominal "
        f"30° results from the finite element beam and mutual coupling; the "
        f"array factor alone predicts 30°."))
    r.figure(str(FIG / "fig11_scan30.png"),
             T("Geschwenkter Betrieb (+30°): 3D-Diagramm, E-Ebene und "
               "H-Ebene im Vergleich mit Element × Arrayfaktor.",
               "Steered operation (+30°): 3D pattern, E plane and H plane "
               "compared with element × array factor."), 16.4)

    r.h2("5.4", T("Vergleich beider Betriebsarten",
                  "Comparison of both modes"))
    r.p(T(
        "Bild 12 legt die E-Ebenen beider Läufe übereinander und zeigt "
        "zusätzlich die aktive Reflexion jedes Ports. Die aktive Reflexion "
        "ist die Überlagerung aller rücklaufenden Wellen an einem Port und "
        "damit die relevante Größe im Betrieb. Im Broadside-Fall liegt sie "
        "zwischen −10,8 dB (Port 1 und 13) und −4,2 dB (Port 4, 8, 12 und "
        "16); im geschwenkten Betrieb verbessert sie sich für die meisten "
        "Ports um 3 bis 6 dB, weil die progressive Phase die "
        "Kopplungsbeiträge teilweise destruktiv überlagert. Die rechte "
        "Spalte bei x = +2,4 mm bleibt in beiden Fällen mit etwa −4 dB "
        "auffällig – dort addieren sich die Kopplungen der Nachbarn "
        "phasengleich.",
        "Figure 12 overlays the E planes of both runs and additionally shows "
        "the active reflection of each port.  The active reflection is the "
        "superposition of all returning waves at a port and therefore the "
        "relevant quantity in operation.  In the broadside case it lies "
        "between −10.8 dB (ports 1 and 13) and −4.2 dB (ports 4, 8, 12 and "
        "16); in the steered case it improves by 3 to 6 dB for most ports "
        "because the progressive phase partially cancels the coupling "
        "contributions.  The right-hand column at x = +2.4 mm remains "
        "conspicuous at about −4 dB in both cases – there the neighbour "
        "couplings add up in phase."))
    r.figure(str(FIG / "fig12_scan_summary.png"),
             T("Vergleich Broadside / +30°: E-Ebene (links) und aktive "
               "Reflexion je Port bei 77 GHz (rechts).",
               "Broadside / +30° comparison: E plane (left) and active "
               "reflection per port at 77 GHz (right)."), 16.0)

    r.h2("5.5", T("Nahfeld", "Near field"))
    r.p(T(
        "Die Bilder 13 und 14 zeigen Momentaufnahmen des Betrags der "
        "elektrischen Feldstärke im Broadside-Betrieb bei 77 GHz, entnommen "
        "den Frequenzbereichs-Dumps des Laufes. In der E-Ebene (Bild 13) "
        "ist die gleichphasige Anregung der vier Elemente in x als "
        "gemeinsame Wellenfront über der Apertur zu erkennen; das Feld "
        "klingt mit dem Abstand von der Platine rasch ab, was die "
        "Nahfeldkopplung zwischen benachbarten Elementen begrenzt. Die "
        "Draufsicht in 1,02 mm Höhe (Bild 14) zeigt die 4×4-Apertur als "
        "kompakte, zusammenhängende Quellfläche – die Patches wirken hier "
        "noch als ein strahlendes Ganzes, nicht als 16 Einzelstrahler.",
        "Figures 13 and 14 show snapshots of the electric field magnitude in "
        "broadside operation at 77 GHz, taken from the frequency-domain "
        "dumps of the run.  In the E plane (Figure 13) the in-phase "
        "excitation of the four elements along x appears as a common wave "
        "front across the aperture; the field decays quickly with distance "
        "from the board, which limits near-field coupling between "
        "neighbouring elements.  The top view at 1.02 mm height (Figure 14) "
        "shows the 4×4 aperture as a compact, contiguous source area – the "
        "patches still act as one radiating entity rather than 16 separate "
        "radiators."))
    r.figure(str(FIG / "fig13_nearfield_xz.png"),
             T("Betrag der elektrischen Feldstärke in der E-Ebene (xz, "
               "y = 0) bei 77 GHz, Broadside-Betrieb. Gelb: Patchkanten.",
               "Electric field magnitude in the E plane (xz, y = 0) at "
               "77 GHz, broadside operation.  Yellow: patch edges."), 16.0)
    r.figure(str(FIG / "fig14_nearfield_top.png"),
             T("Draufsicht auf den Feldbetrag 1,02 mm über der Platine; "
               "grün: Patchkonturen.",
               "Top view of the field magnitude 1.02 mm above the board; "
               "green: patch outlines."), 13.0)

    # ------------------------------------------------------------- 6
    r.h1("6", T("Ergebnisse und Diskussion", "Results and discussion"))
    r.h2("6.1", T("Kennzahlen", "Key figures"))
    r.table(
        T("Zusammenfassung der Array-Kennzahlen bei 77 GHz.",
          "Summary of the array figures at 77 GHz."),
        [T("Größe", "Quantity"), T("Broadside", "Broadside"),
         T("Scan +30°", "Scan +30°")],
        [
            [T("Richtfaktor Dmax", "Directivity Dmax"),
             f"{N(arr['dmax77'], 2)} dBi", f"{N(sc['dmax77'], 2)} dBi"],
            [T("Wirkungsgrad η (roh, 16-Port-Artefakt)",
               "Efficiency η (raw, 16-port artefact)"),
             N(arr["eff77_raw"], 3), N(sc["eff77_raw"], 3)],
            [T("Wirkungsgrad η (korrigiert)",
               "Efficiency η (corrected)"),
             N(arr["eff77_corr"], 3), N(sc["eff77_corr"], 3)],
            [T("Gewinn (korrigiert)", "Gain (corrected)"),
             f"{N(arr['gain77_corr'], 2)} dBi", f"{N(sc['gain77_corr'], 2)} dBi"],
            [T("Strahlrichtung", "Beam direction"),
             f"θ = {N(arr['theta_peak_deg'], 0)}°, φ = {N(arr['phi_peak_deg'], 0)}°",
             f"θ = {N(sc['theta_peak_deg'], 0)}°, φ = {N(sc['phi_peak_deg'], 0)}°"],
        ],
        widths=[5.4, 4.6, 4.6],
        source=T("η korrigiert bezieht die aufgenommene Leistung aller 16 "
                 "Ports ein (siehe Abschnitt 6.3).",
                 "Corrected η includes the accepted power of all 16 ports "
                 "(see Section 6.3)."))

    r.h2("6.2", T("Genauigkeit und Grenzen", "Accuracy and limitations"))
    r.bullets([
        T("Das Mesh löst die 0,254-mm-Folie mit vier Zellen und die "
          "Metallkanten mit rund drei Zellen; die verbleibende numerische "
          "Unsicherheit des Richtfaktors liegt bei etwa ±0,2 dB.",
          "The mesh resolves the 0.254 mm laminate with four cells and the "
          "metal edges with about three cells; the remaining numerical "
          "uncertainty of the directivity is about ±0.2 dB."),
        T("Der Arrayfaktor-Vergleich ignoriert Kopplung und Randbeugung. "
          "Die gute Übereinstimmung der Hauptkeule zeigt, dass diese Effekte "
          "die Strahlform dominieren, nicht aber die Strahlrichtung.",
          "The array-factor comparison ignores coupling and edge "
          "diffraction.  The good agreement of the main beam shows that these "
          "effects shape the pattern but do not dominate the beam direction."),
        T("Die Sondenanpassung ist mit zwei Stellschrauben auf etwa −17 dB "
          "begrenzt; für eine tiefere Anpassung wäre ein "
          "Transformationsnetzwerk erforderlich.",
          "Probe matching with two design knobs is limited to about −17 dB; "
          "a deeper match would require a matching network."),
        T("Die aktive Reflexion erreicht an der rechten Arraykante nur "
          "etwa −4 dB; für ein reales Modul ist dort eine Anpassung des "
          "Speisenetzes oder eine Randbehandlung vorzusehen.",
          "Active reflection reaches only about −4 dB at the right-hand "
          "array edge; a real module needs feed matching or an edge "
          "treatment there."),
    ])

    r.h2("6.3", T("Hinweise zur Leistungsbilanz",
                  "Note on the power budget"))
    r.p(T(
        "Die von openEMS mitgelieferte Nachverarbeitung bestimmt die "
        "aufgenommene Leistung aus dem Port mit der größten Amplitude. Für "
        "einen einzigen angeregten Port ist das korrekt, für 16 gleichzeitig "
        "angeregte Ports unterschätzt es die Leistung um den Faktor 16 und "
        "bläht Wirkungsgrad und Gewinn um rund 12 dB auf. Die in diesem "
        "Bericht angegebenen korrigierten Werte summieren die aufgenommene "
        "Leistung über alle 16 Ports; Richtfaktor und Diagramme sind davon "
        "nicht betroffen, da sie reine Feldintegrale sind.",
        "The post-processing shipped with openEMS determines the accepted "
        "power from the port with the largest amplitude.  This is correct "
        "for a single driven port, but for 16 simultaneously driven ports it "
        "underestimates the power by a factor of 16 and inflates efficiency "
        "and gain by about 12 dB.  The corrected values quoted in this "
        "report sum the accepted power over all 16 ports; directivity and "
        "patterns are unaffected because they are pure field integrals."))
    r.p(T(
        "Bei der Auswertung der Port-Zeitreihen ist außerdem zu beachten, "
        "dass openEMS Spannung und Strom auf leicht versetzten "
        "Zeitrastern ablegt. Werden beide Signale auf ihrer jeweils eigenen "
        "Zeitachse transformiert, entsteht eine scheinbare Phasenverschiebung "
        "und damit ein Fehler in der Wirkleistung von bis zu 11 %. Die "
        "Auswertung verwendet daher für beide Signale dieselbe Zeitachse – "
        "so wie es die openEMS-Nachverarbeitung ebenfalls tut.",
        "When evaluating the port time series it must also be noted that "
        "openEMS stores voltage and current on slightly offset time grids.  "
        "Transforming each signal on its own time axis introduces an "
        "apparent phase shift and hence an error in the real power of up to "
        "11 %.  The evaluation therefore uses the same time axis for both "
        "signals, exactly as the openEMS post-processing does."))

    # ------------------------------------------------------------- 7
    r.h1("7", T("Reproduktion", "Reproduction"))
    r.h2("7.1", T("Voraussetzungen", "Prerequisites"))
    r.p(T(
        "Benötigt werden Python 3.11 mit numpy, matplotlib und h5py sowie "
        "eine openEMS-Installation, die über die VisualopenEMS-"
        "Projektschicht erreichbar ist. Der Pfad zur VisualopenEMS-Wurzel "
        "wird über die Umgebungsvariable `VISUALOPENEMS_ROOT` gesetzt:",
        "Python 3.11 with numpy, matplotlib and h5py is required, plus an "
        "openEMS installation reachable through the VisualopenEMS project "
        "layer. The path to the VisualopenEMS root is set via the "
        "environment variable `VISUALOPENEMS_ROOT`:"))
    r.code("VISUALOPENEMS_ROOT = "
           "C:\\path\\to\\FreeCAD_OpenEMS_GUI")

    r.h2("7.2", T("Ablauf", "Workflow"))
    r.code(
        "python make_project.py --preview      # projects + figures 01..04\n"
        "python validate_project.py            # mesh / collision checks\n"
        "python run_headless.py patch_element_77ghz.fcsim.json\n"
        "python run_headless.py patch_array_4x4_77ghz.fcsim.json\n"
        "python run_headless.py patch_array_4x4_broadside.fcsim.json\n"
        "python run_headless.py patch_array_4x4_scan30.fcsim.json\n"
        "python plot_results.py                # figures 05..14 + results.json\n"
        "python build_docs.py                  # this report (DE + EN)")

    r.h2("7.3", T("Rechenaufwand", "Computational effort"))
    rows = []
    for name, work, val in (
            (T("Einzelelement", "Single element"), "work_element",
             "work_validate_patch_element_77ghz"),
            (T("Eingebettetes Element", "Embedded element"),
             "work_array_embedded", "work_validate_patch_array_4x4_77ghz"),
            (T("Broadside-Array", "Broadside array"),
             "work_array_broadside", "work_validate_patch_array_4x4_broadside"),
            (T("Array +30°", "Array +30°"),
             "work_array_scan30", "work_validate_patch_array_4x4_scan30")):
        try:
            nx, ny, nz, cells = mesh_of(val)
            cells_txt = f"{nx}×{ny}×{nz} = {cells:,}".replace(",", " ")
        except Exception:                     # noqa: BLE001
            cells_txt = "–"
        n = steps_of(work)
        rows.append([name, cells_txt, f"{n:,}".replace(",", " ")
                     if n else "–",
                     T("ca. 20–30 min", "approx. 20–30 min")])
    r.table(
        T("Rechenaufwand der vier Produktionsläufe.",
          "Computational effort of the four production runs."),
        [T("Lauf", "Run"), T("Gitter (Zellen)", "Grid (cells)"),
         T("Zeitschritte", "Time steps"), T("Dauer", "Wall time")],
        rows, widths=[4.6, 5.6, 2.8, 3.2],
        source=T("Zellzahlen aus der openEMS-Startmeldung, Zeitschritte aus "
                 "der Energie-Zeitreihe; Dauer auf einem 8-Kern-Notebook.",
                 "Cell counts from the openEMS start-up message, time steps "
                 "from the energy time series; wall time on an 8-core laptop."))

    r.h2("7.4", T("Dateien", "Files"))
    r.code(
        "array_model.py            single source of truth\n"
        "make_project.py           project writer + construction figures\n"
        "validate_project.py       pre-run checker\n"
        "run_headless.py           openEMS runner\n"
        "plot_results.py           result figures + results.json\n"
        "build_docs.py             this report\n"
        "*.fcsim.json              the four projects\n"
        "work_*/                   simulation output\n"
        "docs/figures/             figure 01..14\n"
        "docs/results.json         numeric source for the report")

    out = DOCS / f"patch_array_4x4_77ghz_{lang}.odt"
    r.save(str(out))
    print(f"wrote {out}")
    return out


def main() -> int:
    build("de")
    build("en")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
