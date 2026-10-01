# Project AEROLAB

A low-cost instrumented wind tunnel built to measure airfoil lift and validate
results against XFOIL and thin-airfoil theory. Designed and built over one month
as a first-year mechanical engineering portfolio project.

---

## What it is

An open-circuit desktop wind tunnel (80 × 80 mm test section) driven by a 120 mm
USB fan, with:
- A **Raspberry Pi + HX711 load cell** force balance measuring lift in real time
- A **pitot tube + U-tube manometer** measuring airspeed
- Three test airfoils: flat plate, NACA 0012, NACA 4412
- Angles of attack swept from −4° to +16° in 2° steps

Measured lift curves are compared against XFOIL polars and thin-airfoil theory
(2π/rad slope), with low-Reynolds-number effects (Re ≈ 25,000) discussed.

---

## Results

<!-- Add your validation plot here once you have it -->
<!-- Drag an image file into this editor on GitHub and it will insert automatically -->

| Airfoil | Measured slope (/deg) | XFOIL slope (/deg) | Stall angle |
|---|---|---|---|
| Flat plate | | | |
| NACA 0012 | | | |
| NACA 4412 | | | |

---

## Repo structure

```
├── 06-Code/
│   ├── daq_lift_logger.py       # Raspberry Pi DAQ — reads HX711, logs to CSV
│   ├── analyse_lift_curve.py    # Python analysis — C_L curves, slope fit, XFOIL overlay
│   └── analyse_lift_curve.m     # MATLAB cross-check
├── 07-CAD/                      # SolidWorks/Fusion files + 1:1 print templates
├── 08-Data/                     # XFOIL polars + experimental lift_runs.csv
└── report/                      # Final write-up
```

---

## Skills demonstrated

CAD · Fabrication · Python (DAQ + analysis) · MATLAB · Aerodynamics ·
Sensors & electronics · Experimental methods · Data analysis · Git

---

## Build cost

~£43 — full bill of materials in [`04-BOM-and-Budget/`](04-BOM-and-Budget/Bill-of-Materials.md)

---

*Oran · MEng Mechanical Engineering · University of Lincoln*
