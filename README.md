<h1 align="center">4×4-Patch-Array für 77 GHz auf RT/duroid 5880</h1>

<p align="center">
  <em>openEMS-Vollwellenstudie einer planaren Gruppenantenne für den automobilen
  Radarbereich – vom Hammerstad-Erstentwurf über das Impedanz-Retuning bis zur
  phasengesteuerten Strahlschwenkung</em>
</p>

<p align="center">
  <img alt="solver" src="https://img.shields.io/badge/solver-openEMS%20(FDTD)-1f5fa8">
  <img alt="band" src="https://img.shields.io/badge/band-74%E2%80%9380%20GHz-0b7285">
  <img alt="python" src="https://img.shields.io/badge/python-3.11-3776AB?logo=python&logoColor=white">
  <img alt="status" src="https://img.shields.io/badge/status-abgeschlossen-success">
</p>

<p align="center">
  <img src="docs/figures/fig04_array3d.png" alt="3D-Modell des 4×4-Arrays" width="760">
</p>

---

## Kurzbeschreibung

Dieses Repository enthält eine vollständig skriptgesteuerte openEMS-Studie
einer 4×4-Patch-Gruppenantenne auf **RT/duroid 5880** für **77 GHz**. Die
Geometrie existiert nur an einer Stelle (`array_model.py`); Projekte, Läufe,
Auswertungen und der PDF/ODT-Bericht werden daraus abgeleitet. Untersucht
werden vier Betriebsfälle:

| Projekt | Inhalt |
|---|---|
| `patch_element_77ghz` | Einzelelement, Produktionsband 74–80 GHz |
| `patch_array_4x4_77ghz` | 4×4-Array, nur Port P6 angeregt (eingebettetes Element, Kopplung) |
| `patch_array_4x4_broadside` | alle 16 Ports phasengleich |
| `patch_array_4x4_scan30` | alle 16 Ports, progressive Phase, Keule bei +30° |

Dazu kommen vier Kontroll- und Verifikationsläufe, die den Anpassungsprozess
dokumentieren (offene Speiseleitung, Inset-Variation, zwei Kandidaten).

## Ergebnisse auf einen Blick

| Größe | Element | Broadside | Scan +30° |
|---|---|---|---|
| Eingangsimpedanz bei 77 GHz | 53,8 + 14,5j Ω | – | – |
| Reflexion \|S11\|(77) | −16,9 dB | – | – |
| −10-dB-Bandbreite | 74,06–79,10 GHz (5,04 GHz) | – | – |
| Richtfaktor D<sub>max</sub> | 9,38 dBi | 17,91 dBi | 17,88 dBi |
| Wirkungsgrad η | 0,939 | 0,949 | 0,970 |
| Gewinn | 9,11 dBi | 17,69 dBi | 17,75 dBi |
| Keulenrichtung | +z | θ = 2°, φ = 0° | θ = 28°, φ = 0° |
| Aktive Reflexion je Port | – | −10,8 … −4,2 dB | −14,4 … −3,9 dB |

> **Hinweis zur Leistungsbilanz:** Die mitgelieferte openEMS-Nachverarbeitung
> bestimmt die aufgenommene Leistung nur aus dem amplitudenstärksten Port. Bei
> 16 gleichzeitig gespeisten Ports unterschätzt das die Leistung um den Faktor
> 16 und bläht η/Gewinn um ~12 dB auf. Die Tabelle enthält die **korrigierten**
> Werte (Summe über alle 16 Ports); Richtfaktor und Diagramme sind davon
> unabhängig. Details in Abschnitt [Hinweise zur Auswertung](#hinweise-zur-auswertung).

---

## 1. Aufbau der Antenne

### 1.1 Layout

Vier Elemente in x und vier in y, Rastermaß d = λ<sub>0</sub>/2 = 1,947 mm,
Platine 12 × 12 mm. Jedes Element hat einen eigenen koaxialen Sondenport von
der Grundfläche zur Patchunterseite, sodass alle 16 Elemente unabhängig
gespeist werden können.

<p align="center">
  <img src="docs/figures/fig01_layout.png" alt="Array-Layout" width="620">
</p>

### 1.2 Einzelelement und Speisung

L = 1,161 mm (resonant in x), W = 1,539 mm. Der Sonden-Inset beträgt im
finalen Entwurf nur **0,058 mm** – rund 5 % der Patchlänge. Rechts im Bild ist
deshalb ein Zoom auf die Speisung dargestellt; auf der Gesamtansicht wäre der
Unterschied nicht auflösbar.

<p align="center">
  <img src="docs/figures/fig02_element.png" alt="Einzelelement mit Zoom auf die Speisung" width="860">
</p>

### 1.3 Schichtaufbau

RT/duroid 5880 mit ε<sub>r</sub> = 2,2, tan δ = 0,0009 und h = 0,254 mm
(10 mil); beide Kupferlagen als 35 µm dicke PEC-Schichten. Die z-Achse ist im
Schnitt 6-fach überhöht.

<p align="center">
  <img src="docs/figures/fig03_stackup.png" alt="Schichtaufbau" width="780">
</p>

---

## 2. Anpassung des Einzelelements

Der Hammerstad-Erstentwurf (L = 1,1383 mm, Inset 0,30 L = 0,341 mm) liefert
bei 77 GHz nur **14,8 + 68,0j Ω** und **−1,77 dB** – eine stark induktive
Fehlanpassung. Ein Kontrollexperiment mit offener Speiseleitung ohne Patch
zeigt, dass die Sonde selbst eine kurze, offene Leitung mit
**Z<sub>p</sub> = 188 Ω, ε<sub>p</sub> = 3,40** ist (−205j Ω bei 77 GHz) und
keineswegs als Induktivität wirkt. Mit diesem Leitungsmodell lässt sich der
Sondeneinfluss herausrechnen; ein lokales Newton-Verfahren führt anschließend
über die beiden Stellschrauben Patchlänge und Inset auf 50 Ω.

<p align="center">
  <img src="docs/figures/fig06_retune_journey.png" alt="Kontrollexperiment und Zin-Ortskurven" width="880">
</p>

**Entwurfskette** (simulierte Werte bei 77 GHz):

| Schritt | L (mm) | Inset | Z<sub>in</sub>(77) | \|S11\|(77) |
|---|---|---|---|---|
| Hammerstad-Erstentwurf | 1,1383 | 0,30 L | 14,8 + 68,0j Ω | −1,8 dB |
| Inset-Kontrolle | 1,1383 | 0,05 L | 66,6 + 10,8j Ω | −15,4 dB |
| Kandidat B | 1,1683 | 0,02 L | 38,7 + 20,7j Ω | −11,7 dB |
| **Kandidat A – final** | **1,1610** | **0,05 L** | **53,8 + 14,5j Ω** | **−16,9 dB** |

Das Optimum liegt mit −17,4 dB bei 76,58 GHz. Mehr als etwa −17 dB sind mit
diesen zwei Stellschrauben bei dieser Substratdicke nicht erreichbar: Die
Impedanzkarte ist nahe am Optimum fast kollinear, Länge und Inset verschieben
die Impedanz in dieselbe Richtung.

<p align="center">
  <img src="docs/figures/fig05_element_zin_s11.png" alt="Eingangsimpedanz und Reflexion" width="880">
</p>

---

## 3. Strahlungsverhalten

### 3.1 Isoliertes Element

Richtfaktor 9,38 dBi, Gewinn 9,11 dBi, Wirkungsgrad 0,939. Die Keule steht
erwartungsgemäß in +z, die E-Ebene ist die xz-Ebene.

<p align="center">
  <img src="docs/figures/fig07_element_pattern3d.png" alt="3D-Richtdiagramm des Elements" width="560">
</p>

### 3.2 Eingebettetes Element und Kopplung

Im Lauf `embedded` wird nur Port P6 angeregt, alle anderen Ports sind passiv
mit 50 Ω abgeschlossen. Die Eigenreflexion steigt auf **−18,0 dB**, der
stärkste Koppler ist das direkte Nachbarelement in x mit **−14,3 dB** bei
1,95 mm; die Diagonalen liegen bei etwa −25 dB, die Randelemente unter −28 dB.
Das eingebettete Element erreicht 9,16 dBi bei η = 0,814 – die Differenz zum
isolierten Element ist die in die Nachbarn abgegebene Leistung, die im
voll gespeisten Array wieder abgestrahlt wird.

<p align="center">
  <img src="docs/figures/fig09_coupling.png" alt="Kopplung im eingebetteten Lauf" width="880">
</p>

<p align="center">
  <img src="docs/figures/fig08_element_vs_embedded.png" alt="Isoliertes vs. eingebettetes Element" width="880">
</p>

---

## 4. Das 4×4-Array

### 4.1 Broadside

Alle 16 Ports phasengleich: **17,91 dBi** Richtfaktor. Das passt gut zur
Aperturabschätzung D ≈ 4πA/λ<sub>0</sub>² = 17,0 dBi für eine gleichförmig
belegte λ<sub>0</sub>/2-Apertur; der um 0,9 dB höhere Wert folgt aus der
Elementcharakteristik. Der Vergleich mit Element × Arrayfaktor stimmt bis
−20 dB hervorragend überein – Abweichungen in den Nullstellen und weiten
Nebenkeulen zeigen die vernachlässigte Kopplung und Randbeugung.

<p align="center">
  <img src="docs/figures/fig10_broadside.png" alt="Broadside-Diagramm" width="900">
</p>

### 4.2 Strahlschwenkung um +30°

Progressiv −90° je Element in x (entspricht 90° pro Raster d = λ<sub>0</sub>/2).
Das Vorzeichen folgt der openEMS-Phasenkonvention (eine in +x laufende Welle
trägt exp(−jkx), die Anregung muss voreilen) und wurde mit einem Testlauf
verifiziert: Mit umgekehrtem Vorzeichen erscheint die Keule spiegelbildlich
bei φ = 180°.

Das Maximum liegt bei **θ = 28°, φ = 0°** mit 17,75 dBi; die Ablage gegenüber
den nominellen 30° folgt aus der endlichen Elementkeule und der Kopplung.

<p align="center">
  <img src="docs/figures/fig11_scan30.png" alt="Geschwenkte Keule" width="900">
</p>

### 4.3 Vergleich der Betriebsarten

Die aktive Reflexion liegt im Broadside-Fall zwischen −10,8 dB (Port 1/13) und
−4,2 dB (Port 4/8/12/16). Im geschwenkten Betrieb verbessert sie sich für die
meisten Ports um 3–6 dB, weil die progressive Phase die Kopplungsbeiträge
teilweise destruktiv überlagert; die rechte Spalte bei x = +2,4 mm bleibt mit
etwa −4 dB in beiden Fällen auffällig.

<p align="center">
  <img src="docs/figures/fig12_scan_summary.png" alt="Vergleich Broadside / Scan" width="880">
</p>

---

## 5. Nahfeld

Frequenzbereichs-Snapshots bei 77 GHz aus dem Broadside-Lauf. In der E-Ebene
ist die gleichphasige Anregung der vier x-Elemente als gemeinsame Wellenfront
zu erkennen; die Draufsicht 1,02 mm über der Platine zeigt die 4×4-Apertur als
zusammenhängende Quellfläche.

<p align="center">
  <img src="docs/figures/fig13_nearfield_xz.png" alt="Nahfeld E-Ebene" width="760">
</p>

<p align="center">
  <img src="docs/figures/fig14_nearfield_top.png" alt="Nahfeld Draufsicht" width="600">
</p>

---

## 6. Diskretisierung und Rechenaufwand

Mesh-Preset „high“ mit expliziten Zellgrößen: max. 0,25 mm, verfeinerte Kante
0,10 mm, min. 0,07 mm, min. z 0,06 mm. Die 0,254-mm-Folie wird von vier Zellen
in z aufgelöst, Metallkanten von rund drei Zellen. Abbruch bei −40 dB
Restenergie.

| Lauf | Gitter (Zellen) | Zeitschritte | Dauer |
|---|---|---|---|
| Einzelelement | 154×155×75 = 1 790 250 | 10 506 | ca. 20–30 min |
| Eingebettetes Element | 163×165×75 = 2 017 125 | 10 598 | ca. 20–30 min |
| Broadside-Array | 163×165×75 = 2 017 125 | 10 598 | ca. 20–30 min |
| Array +30° | 163×165×75 = 2 017 125 | 10 598 | ca. 20–30 min |

---

## 7. Reproduktion

**Voraussetzungen:** Python 3.11 mit numpy, matplotlib, h5py und Pillow sowie
eine openEMS-Installation, die über die VisualopenEMS-Projektschicht
erreichbar ist:

```powershell
$env:VISUALOPENEMS_ROOT = "C:\path\to\FreeCAD_OpenEMS_GUI"
```

**Kompletter Ablauf:**

```powershell
python make_project.py --preview      # Projekte + Konstruktionsbilder 01..04
python validate_project.py            # Mesh-/Kollisionsprüfung
python run_headless.py patch_element_77ghz.fcsim.json
python run_headless.py patch_array_4x4_77ghz.fcsim.json
python run_headless.py patch_array_4x4_broadside.fcsim.json
python run_headless.py patch_array_4x4_scan30.fcsim.json
python plot_results.py                # Bilder 05..14 + docs/results.json
python build_docs.py                  # Bericht DE + EN
```

Oder alles auf einmal (überspringt bereits vorhandene Ergebnisse):

```powershell
python run_all.py
```

---

## 8. Projektstruktur

```
array_model.py            einzige Quelle der Geometrie, Ports, Mesh, NF2FF-Box
make_project.py           schreibt die .fcsim.json-Projekte + Bilder 01..04
validate_project.py       prüft Kollisionen, Mesh-Regeln, NF2FF-Box
run_headless.py           erzeugt das openEMS-Skript und startet den Lauf
plot_results.py           wertet S-Parameter und Fernfelder aus (Bilder 05..14)
build_docs.py             baut den deutschen und englischen ODT-Bericht
run_all.py                Gesamtlauf mit Skip-Logik
analyze_zin.py            Zin(f) schnell auswerten
fit_model.py              Leitungsmodell der Sonde + Newton-Retuning
make_candidate.py         Kandidatenprojekte erzeugen
analyze_controls.py       Kontrollläufe vergleichen
make_controls.py          Kontrollprojekte erzeugen (probe_only, inset005)
make_wideband.py          Breitband-Referenzlauf
*.fcsim.json              die Projekte (Produktion + Kontrollen)
work_*/                   Simulationsausgaben (nicht versioniert)
docs/figures/             Bilder 01..14
docs/results.json         Zahlenquelle des Berichts
docs/*.odt                Bericht Deutsch + Englisch
```

---

## 9. Dokumentation

| Datei | Inhalt |
|---|---|
| [`docs/patch_array_4x4_77ghz_de.odt`](docs/patch_array_4x4_77ghz_de.odt) | 15-seitiger Bericht, 14 Bilder, 4 Tabellen |
| [`docs/patch_array_4x4_77ghz_en.odt`](docs/patch_array_4x4_77ghz_en.odt) | englische Fassung |
| [`docs/results.json`](docs/results.json) | alle Kennzahlen als JSON |

## Hinweise zur Auswertung

* **Wirkungsgrad/Gewinn der 16-Port-Läufe.** Die gebündelte
  openEMS-Nachverarbeitung nutzt nur den amplitudenstärksten Port
  (`postproc.nf2ff._load_pin_spectral`). Bei 16 gleichzeitig gespeisten Ports
  ist Pin dadurch ~16× zu klein, η und Gewinn sind ~12 dB zu hoch.
  `plot_results.py` summiert die aufgenommene Leistung über alle Ports;
  Richtfaktor und Diagramme sind reine Feldintegrale und bleiben unberührt.
* **Zeitachsen der Portsignale.** openEMS legt Spannung und Strom auf leicht
  versetzten Zeitrastern ab (≈0,03–0,04 dt). Werden beide Signale auf ihrer
  eigenen Achse transformiert, entsteht eine scheinbare Phasenverschiebung und
  ein Fehler in der Wirkleistung von bis zu 11 %. Die Auswertung verwendet
  deshalb – wie die openEMS-Nachverarbeitung – dieselbe Zeitachse für beide
  Signale.
* **Phasenvorzeichen.** „+30°“ bedeutet Strahlrichtung +x (φ = 0°). Die
  Konvention wurde mit einem Testlauf verifiziert.

## Genauigkeit und Grenzen

* Numerische Unsicherheit des Richtfaktors ≈ ±0,2 dB.
* Der Arrayfaktor-Vergleich ignoriert Kopplung und Randbeugung – er beschreibt
  die Strahlform gut, die Strahlrichtung nur näherungsweise.
* Die Sondenanpassung ist mit zwei Stellschrauben auf ≈ −17 dB begrenzt.
* Die rechte Arraykante zeigt eine aktive Reflexion von nur ≈ −4 dB; ein
  reales Modul bräuchte dort eine Anpassung des Speisenetzes.

---

# English summary

A fully script-driven **openEMS FDTD** study of a planar **4×4 patch array at
77 GHz on RT/duroid 5880** (ε<sub>r</sub> = 2.2, tan δ = 0.0009, h = 0.254 mm).

* **Element** (L = 1.161 mm, W = 1.539 mm, probe inset 0.058 mm):
  Z<sub>in</sub>(77) = 53.8 + 14.5j Ω, |S11| = −16.9 dB, −10 dB bandwidth
  74.06–79.10 GHz, D<sub>max</sub> = 9.38 dBi, η = 0.939.
* **Retuning story:** the Hammerstad first cut gives 14.8 + 68.0j Ω
  (−1.8 dB). A control run with an open feed line shows that the probe is a
  short open line (Z<sub>p</sub> = 188 Ω, ε<sub>p</sub> = 3.40), not an
  inductance. After de-embedding it, a Newton iteration on patch length and
  inset reaches the final design.
* **Broadside array:** D<sub>max</sub> = 17.91 dBi, η = 0.949, gain
  17.69 dBi, beam at θ = 2°.
* **Phased array (+30°):** progressive −90° per element, beam at
  θ = 28°, φ = 0°, gain 17.75 dBi, η = 0.970.
* **Reports:** `docs/patch_array_4x4_77ghz_de.odt` and `..._en.odt`
  (15 pages, 14 figures, 4 tables each).

Run `python run_all.py` to regenerate everything; simulation output
(`work_*/`) is not versioned.

---

<p align="center"><sub>VisualopenEMS · openEMS 2026</sub></p>
