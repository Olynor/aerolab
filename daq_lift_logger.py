#!/usr/bin/env python3
"""
daq_lift_logger.py  --  Project AEROLAB

Reads a load cell through an HX711 amplifier on a Raspberry Pi and logs the
LIFT force on a wind-tunnel airfoil against angle of attack.

Three modes:
    raw        -- print live raw counts + converted force (sanity-check wiring)
    calibrate  -- find the scale factor from known masses (do this once)
    log        -- run an angle-of-attack sweep and append rows to a CSV

Use --simulate to run on a laptop with NO hardware (fake readings) so you can
practise the whole workflow before your parts arrive.

Examples
--------
    python3 daq_lift_logger.py raw
    python3 daq_lift_logger.py calibrate
    python3 daq_lift_logger.py log
    python  daq_lift_logger.py --simulate log     # no Pi needed

This file is intentionally heavily commented -- it doubles as Python practice
for the "Python for Engineering" certification.
"""

import argparse
import csv
import os
import random
import statistics
import sys
import time
from datetime import datetime

# ----------------------------------------------------------------------
# CONFIG  -- edit these for your build
# ----------------------------------------------------------------------
DT_PIN = 5            # HX711 DT  -> Pi GPIO5  (BCM numbering)
SCK_PIN = 6           # HX711 SCK -> Pi GPIO6
SAMPLES = 15          # readings averaged per data point
GAIN = 128            # HX711 channel-A gain

# From calibration (see calibrate mode and 03-Notes/Calibration-Records.md):
SCALE_FACTOR = 1000.0   # raw counts per NEWTON   <-- REPLACE after calibrating
OFFSET = 0.0            # raw counts at zero load (tare) <-- REPLACE

# Where the data goes (relative to this script, i.e. ../08-Data/)
HERE = os.path.dirname(os.path.abspath(__file__))
CSV_PATH = os.path.normpath(os.path.join(HERE, "..", "08-Data", "lift_runs.csv"))

G = 9.81  # m/s^2, for converting masses to forces during calibration


# ----------------------------------------------------------------------
# Reader abstraction: real hardware OR a simulator
# ----------------------------------------------------------------------
class SimReader:
    """Fake load cell so you can test everything without a Pi.

    Pretends there is a tiny bit of lift plus random noise so the numbers
    look believable. Not physically accurate -- just for practising the code.
    """

    def __init__(self):
        self._base = 0.0

    def read_raw(self):
        # wander slowly + add noise so it behaves like a real noisy sensor
        self._base += random.uniform(-50, 50)
        return OFFSET + self._base + random.gauss(0, 300)

    def cleanup(self):
        pass


class HX711Reader:
    """Thin wrapper around an HX711 Python library on the Raspberry Pi.

    NOTE: different HX711 libraries use slightly different method names.
    This is written for the common `hx711` API. If yours differs, the only
    line you should need to change is marked  <<< ADJUST  below.
    """

    def __init__(self, dt=DT_PIN, sck=SCK_PIN, gain=GAIN):
        try:
            from hx711 import HX711  # <<< ADJUST import to match your library
        except Exception as exc:  # pragma: no cover - hardware only
            print("Could not import an HX711 library:", exc)
            print("Install one (see 06-Code/README.md) or run with --simulate.")
            sys.exit(1)
        self._hx = HX711(dout_pin=dt, pd_sck_pin=sck, gain=gain)
        try:
            self._hx.reset()
        except Exception:
            pass

    def read_raw(self):
        # <<< ADJUST this call to match your library's "read N and average" method
        value = self._hx.get_raw_data_mean(readings=1)
        if value is False or value is None:
            raise RuntimeError("HX711 returned no data -- check wiring/power.")
        return float(value)

    def cleanup(self):
        try:
            import RPi.GPIO as GPIO
            GPIO.cleanup()
        except Exception:
            pass


def get_reader(simulate):
    return SimReader() if simulate else HX711Reader()


# ----------------------------------------------------------------------
# Helpers
# ----------------------------------------------------------------------
def averaged_reading(reader, n=SAMPLES):
    """Return (mean, std) of n raw readings."""
    vals = []
    for _ in range(max(1, n)):
        vals.append(reader.read_raw())
        time.sleep(0.02)
    mean = statistics.fmean(vals)
    std = statistics.pstdev(vals) if len(vals) > 1 else 0.0
    return mean, std


def raw_to_newtons(raw):
    return (raw - OFFSET) / SCALE_FACTOR


def ensure_csv_header(path):
    new = not os.path.exists(path)
    os.makedirs(os.path.dirname(path), exist_ok=True)
    if new:
        with open(path, "w", newline="") as f:
            csv.writer(f).writerow([
                "timestamp", "airfoil", "airspeed_ms", "alpha_deg",
                "lift_g", "lift_N", "raw_mean", "raw_std", "samples", "notes",
            ])


# ----------------------------------------------------------------------
# Modes
# ----------------------------------------------------------------------
def mode_raw(reader):
    print("Live readings (Ctrl-C to stop). Press up on the airfoil to test sign.\n")
    print(f"{'raw':>12} {'force (N)':>12} {'force (g)':>12}")
    try:
        while True:
            raw = reader.read_raw()
            n = raw_to_newtons(raw)
            print(f"{raw:12.0f} {n:12.4f} {n / G * 1000:12.2f}")
            time.sleep(0.2)
    except KeyboardInterrupt:
        print("\nStopped.")


def mode_calibrate(reader):
    """Interactive calibration: tare, then add known masses, then linear fit."""
    print("=== CALIBRATION ===")
    print("Mount the airfoil + sting. Keep the fan OFF.\n")
    input("Step 1: remove all load, then press Enter to TARE...")
    offset, _ = averaged_reading(reader, n=SAMPLES * 2)
    print(f"  tare (offset) = {offset:.1f} counts\n")

    masses_g, raws = [], []
    print("Step 2: add known masses one at a time (e.g. coins; verify on a")
    print("kitchen scale). Enter the TOTAL mass in grams each time, or blank to finish.\n")
    while True:
        s = input("  total mass on sting (g), or Enter to finish: ").strip()
        if s == "":
            break
        try:
            m = float(s)
        except ValueError:
            print("    not a number, try again")
            continue
        raw, std = averaged_reading(reader, n=SAMPLES * 2)
        masses_g.append(m)
        raws.append(raw)
        print(f"    raw = {raw:.1f}  (std {std:.1f})")

    if len(masses_g) < 2:
        print("\nNeed at least two masses for a fit. Aborting.")
        return

    # Linear fit: raw = k * force_N + offset  ->  scale_factor k = counts per N
    forces_N = [m / 1000.0 * G for m in masses_g]
    try:
        import numpy as np
        k, c = np.polyfit(forces_N, [r - offset for r in raws], 1)
        # R^2
        pred = [k * f + c for f in forces_N]
        ss_res = sum((a - b) ** 2 for a, b in zip([r - offset for r in raws], pred))
        mean_y = statistics.fmean([r - offset for r in raws])
        ss_tot = sum((y - mean_y) ** 2 for y in [r - offset for r in raws])
        r2 = 1 - ss_res / ss_tot if ss_tot else float("nan")
    except ImportError:
        # fallback: simple two-point slope through the extremes
        k = (raws[-1] - raws[0]) / (forces_N[-1] - forces_N[0])
        r2 = float("nan")

    print("\n=== RESULT ===")
    print(f"  SCALE_FACTOR = {k:.2f}   (raw counts per newton)")
    print(f"  OFFSET       = {offset:.2f}")
    if r2 == r2:  # not NaN
        print(f"  linearity R^2 = {r2:.4f}   (want > 0.99)")
    print("\nCopy SCALE_FACTOR and OFFSET into the CONFIG block at the top of")
    print("this file, and record them in 03-Notes/Calibration-Records.md.")


def mode_log(reader):
    """Run an angle-of-attack sweep, appending one row per angle to the CSV."""
    print("=== LOG RUN ===")
    airfoil = input("Airfoil name (e.g. naca4412): ").strip() or "unknown"
    airspeed = input("Test airspeed in m/s (from the manometer): ").strip()
    try:
        airspeed_ms = float(airspeed)
    except ValueError:
        airspeed_ms = float("nan")
        print("  (couldn't parse airspeed; logged as blank -- add it later)")

    ensure_csv_header(CSV_PATH)
    print("\nFan ON, airspeed steady. Enter each angle of attack; blank to finish.")
    print("Tip: tare with the fan OFF first if you bumped anything.\n")

    while True:
        s = input("  angle of attack (deg), or Enter to finish: ").strip()
        if s == "":
            break
        try:
            alpha = float(s)
        except ValueError:
            print("    not a number, try again")
            continue
        note = input("    note (optional): ").strip()
        print(f"    averaging {SAMPLES} samples...")
        raw, std = averaged_reading(reader, n=SAMPLES)
        lift_N = raw_to_newtons(raw)
        lift_g = lift_N / G * 1000.0
        row = [
            datetime.now().isoformat(timespec="seconds"),
            airfoil, airspeed_ms, alpha,
            round(lift_g, 3), round(lift_N, 5),
            round(raw, 1), round(std, 1), SAMPLES, note,
        ]
        with open(CSV_PATH, "a", newline="") as f:
            csv.writer(f).writerow(row)
        print(f"    logged: alpha={alpha} deg  lift={lift_g:.2f} g  ({lift_N:.4f} N)\n")

    print(f"Saved to {CSV_PATH}")


# ----------------------------------------------------------------------
# Entry point
# ----------------------------------------------------------------------
def main():
    p = argparse.ArgumentParser(description="AEROLAB load-cell DAQ / logger")
    p.add_argument("mode", choices=["raw", "calibrate", "log"],
                   help="raw=live readings, calibrate=find scale factor, log=sweep")
    p.add_argument("--simulate", action="store_true",
                   help="use fake readings (no Raspberry Pi / hardware needed)")
    args = p.parse_args()

    if args.simulate:
        print(">>> SIMULATE MODE: fake readings, no hardware.\n")
    reader = get_reader(args.simulate)
    try:
        {"raw": mode_raw, "calibrate": mode_calibrate, "log": mode_log}[args.mode](reader)
    finally:
        reader.cleanup()


if __name__ == "__main__":
    main()
