#!/usr/bin/env python3
"""
analyse_lift_curve.py  --  Project AEROLAB

Turns the logged data (../08-Data/lift_runs.csv) into lift-coefficient vs
angle-of-attack curves, fits the lift-curve slope, finds the stall point,
and (optionally) overlays XFOIL polars for comparison.

    python analyse_lift_curve.py            # analyse your real data
    python analyse_lift_curve.py --demo     # synthesise fake data + plot (no CSV needed)

Needs: numpy, pandas, matplotlib   (pip install numpy pandas matplotlib)

Heavily commented on purpose -- good practice for the Python for Engineering cert.
"""

import argparse
import os
import sys

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# ----------------------------------------------------------------------
# CONFIG -- match your airfoil + air properties
# ----------------------------------------------------------------------
CHORD = 0.060      # m
SPAN = 0.078       # m
RHO_AIR = 1.225    # kg/m^3
S = CHORD * SPAN   # reference area, m^2
MU = 1.81e-5       # Pa.s, air viscosity (for Reynolds number)

LINEAR_REGION_DEG = (-2.0, 8.0)  # alpha range used to fit the lift-curve slope

HERE = os.path.dirname(os.path.abspath(__file__))
CSV_PATH = os.path.normpath(os.path.join(HERE, "..", "08-Data", "lift_runs.csv"))
DATA_DIR = os.path.dirname(CSV_PATH)
FIG_PATH = os.path.join(DATA_DIR, "lift_curve.png")


# ----------------------------------------------------------------------
# Demo data (only used with --demo)
# ----------------------------------------------------------------------
def make_demo_dataframe():
    """Synthesise believable (NOT real) data for three airfoils."""
    rng = np.random.default_rng(0)
    airspeed = 7.0
    alphas = np.arange(-4, 17, 2.0)
    # (name, zero-lift angle deg, slope per deg, stall angle deg, CLmax)
    specs = [
        ("flat_plate", 0.0, 0.075, 10, 0.75),
        ("naca0012", 0.0, 0.090, 12, 0.95),
        ("naca4412", -3.5, 0.095, 11, 1.15),
    ]
    rows = []
    for name, a0, slope, stall, clmax in specs:
        for a in alphas:
            cl = slope * (a - a0)
            if a > stall:  # crude post-stall drop
                cl = clmax - 0.05 * (a - stall)
            cl = min(cl, clmax)
            q = 0.5 * RHO_AIR * airspeed ** 2
            lift_N = cl * q * S
            for _ in range(2):  # two repeats per angle
                noise = rng.normal(0, 0.06 * abs(lift_N) + 0.0008)
                L = lift_N + noise
                rows.append({
                    "airfoil": name, "airspeed_ms": airspeed, "alpha_deg": a,
                    "lift_N": L, "lift_g": L / 9.81 * 1000,
                })
    return pd.DataFrame(rows)


# ----------------------------------------------------------------------
# XFOIL polar parsing (optional overlay)
# ----------------------------------------------------------------------
def load_xfoil_polar(path):
    """Robustly read an XFOIL polar file -> (alpha[], CL[]).

    XFOIL files have a header block, then a row of column names containing
    'alpha' and 'CL', then numeric rows. We find that header row, map the
    columns, and read the numbers.
    """
    with open(path) as f:
        lines = f.readlines()
    hdr_idx = None
    for i, line in enumerate(lines):
        low = line.lower()
        if "alpha" in low and "cl" in low:
            hdr_idx = i
            break
    if hdr_idx is None:
        return None
    cols = lines[hdr_idx].split()
    cols_low = [c.lower() for c in cols]
    try:
        ia, icl = cols_low.index("alpha"), cols_low.index("cl")
    except ValueError:
        ia, icl = 0, 1
    alpha, cl = [], []
    for line in lines[hdr_idx + 1:]:
        parts = line.split()
        if len(parts) <= max(ia, icl):
            continue
        try:
            alpha.append(float(parts[ia]))
            cl.append(float(parts[icl]))
        except ValueError:
            continue
    if not alpha:
        return None
    order = np.argsort(alpha)
    return np.array(alpha)[order], np.array(cl)[order]


def find_polar_for(airfoil):
    """Look for ../08-Data/polar_<airfoil>.txt (case-insensitive-ish)."""
    candidates = [
        os.path.join(DATA_DIR, f"polar_{airfoil}.txt"),
        os.path.join(DATA_DIR, f"polar_{airfoil.lower()}.txt"),
        os.path.join(DATA_DIR, f"{airfoil}_polar.txt"),
    ]
    for c in candidates:
        if os.path.exists(c):
            return c
    return None


# ----------------------------------------------------------------------
# Analysis
# ----------------------------------------------------------------------
def compute_cl(df):
    """Add a C_L column from lift force and airspeed (per row)."""
    if "lift_N" not in df.columns:
        if "lift_g" in df.columns:
            df["lift_N"] = df["lift_g"] / 1000.0 * 9.81
        else:
            sys.exit("CSV needs a 'lift_N' or 'lift_g' column.")
    df = df.copy()
    df["airspeed_ms"] = pd.to_numeric(df["airspeed_ms"], errors="coerce")
    if df["airspeed_ms"].isna().all():
        sys.exit("No airspeed in the data -- can't form a coefficient. "
                 "Add airspeed_ms (m/s) to the CSV.")
    q = 0.5 * RHO_AIR * df["airspeed_ms"] ** 2
    df["CL"] = df["lift_N"] / (q * S)
    return df


def summarise(df):
    """Average repeats at each (airfoil, alpha); return mean + std of C_L."""
    g = df.groupby(["airfoil", "alpha_deg"])["CL"]
    out = g.agg(["mean", "std", "count"]).reset_index()
    out["std"] = out["std"].fillna(0.0)
    return out


def fit_slope(alpha, cl):
    """Fit C_L = slope*alpha + intercept over the linear region. Returns (slope_per_deg, intercept, alpha_L0)."""
    lo, hi = LINEAR_REGION_DEG
    mask = (alpha >= lo) & (alpha <= hi)
    if mask.sum() < 2:
        return None
    slope, intercept = np.polyfit(alpha[mask], cl[mask], 1)
    alpha_L0 = -intercept / slope if slope != 0 else float("nan")
    return slope, intercept, alpha_L0


def reynolds(airspeed):
    return RHO_AIR * airspeed * CHORD / MU


# ----------------------------------------------------------------------
# Plot + report
# ----------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser(description="AEROLAB lift-curve analysis")
    ap.add_argument("--demo", action="store_true",
                    help="use synthetic data if no real CSV is available")
    ap.add_argument("--csv", default=CSV_PATH, help="path to lift_runs.csv")
    args = ap.parse_args()

    synthetic = False
    if args.demo or not os.path.exists(args.csv):
        if not args.demo:
            print(f"No data file at {args.csv} -- showing a DEMO with synthetic data.")
        df = make_demo_dataframe()
        synthetic = True
    else:
        df = pd.read_csv(args.csv)

    df = compute_cl(df)
    summary = summarise(df)

    airspeed = float(np.nanmean(pd.to_numeric(df["airspeed_ms"], errors="coerce")))
    print("\n=== Lift-curve analysis ===")
    print(f"Reference area S = {S:.5f} m^2 (chord {CHORD*1000:.0f} mm x span {SPAN*1000:.0f} mm)")
    print(f"Mean test airspeed ~ {airspeed:.2f} m/s  ->  Re ~ {reynolds(airspeed):,.0f}")
    if synthetic:
        print("** SYNTHETIC DEMO DATA -- not real measurements **")
    print()

    plt.figure(figsize=(8, 6))
    colours = plt.cm.viridis(np.linspace(0, 0.85, summary["airfoil"].nunique()))

    for (name, sub), col in zip(summary.groupby("airfoil"), colours):
        sub = sub.sort_values("alpha_deg")
        a = sub["alpha_deg"].to_numpy()
        cl = sub["mean"].to_numpy()
        err = sub["std"].to_numpy()

        plt.errorbar(a, cl, yerr=err, fmt="o-", color=col, capsize=3,
                     label=f"{name} (measured)")

        fit = fit_slope(a, cl)
        if fit:
            slope, intercept, aL0 = fit
            lo, hi = LINEAR_REGION_DEG
            xs = np.linspace(lo, hi, 50)
            plt.plot(xs, slope * xs + intercept, "--", color=col, alpha=0.6)
            i_max = int(np.argmax(cl))
            print(f"{name}:")
            print(f"   lift-curve slope = {slope:.4f} /deg  ({slope*180/np.pi:.3f} /rad)")
            print(f"   (ideal thin-airfoil  = 0.1097 /deg  = 2*pi /rad)")
            print(f"   zero-lift angle  ~ {aL0:.2f} deg")
            print(f"   max C_L          ~ {cl[i_max]:.3f} at alpha = {a[i_max]:.1f} deg (stall)")
            print()

        # optional XFOIL overlay
        polar = find_polar_for(name)
        if polar:
            res = load_xfoil_polar(polar)
            if res is not None:
                xa, xcl = res
                plt.plot(xa, xcl, ":", color=col, linewidth=2,
                         label=f"{name} (XFOIL)")

    # ideal thin-airfoil reference line through origin
    xs = np.linspace(*LINEAR_REGION_DEG, 50)
    plt.plot(xs, 0.1097 * xs, "k-", alpha=0.3, linewidth=1,
             label="thin-airfoil ideal (2π/rad)")

    plt.axhline(0, color="grey", linewidth=0.6)
    plt.axvline(0, color="grey", linewidth=0.6)
    plt.xlabel("Angle of attack  α  (deg)")
    plt.ylabel("Lift coefficient  $C_L$")
    title = "AEROLAB: lift curves"
    if synthetic:
        title += "  [SYNTHETIC DEMO]"
    plt.title(title)
    plt.legend(fontsize=8)
    plt.grid(True, alpha=0.3)
    plt.tight_layout()

    os.makedirs(DATA_DIR, exist_ok=True)
    plt.savefig(FIG_PATH, dpi=150)
    print(f"Figure saved to {FIG_PATH}")
    try:
        plt.show()
    except Exception:
        pass


if __name__ == "__main__":
    main()
