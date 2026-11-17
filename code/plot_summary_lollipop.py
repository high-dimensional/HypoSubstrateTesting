#########################################################################
#                                      ###     ###    #######   ###     #
#                                      ###     ###   ###        ###     #
#                                      ###     ###   ###        ###     #
#                                      ###     ###   ###        ###     #
#                                       #########     #######   #########
#                                                                       #
#                                                                       #
#                      HypoSubstrateTesting                             #
# Testing hypothesized neural substrates with lesions                   #
#                                                                       #
# This script reads the *_summary.csv files produced by                 #
# run_rldr_simulation.m and draws the radial lollipop figure: four      #
# circles with the percentage of iterations that are significant (p <   #
# 0.05) per noise bin for four compartments each.                       #
# Rows: cortical networks (top), white matter tracts (bottom).          #
# Columns: false positive studies (left, global generator) and true     #
# positive studies (right, substrate as behaviour generator).           #
# Blue circles: lesion volume-unadjusted,                               #
# green triangles: lesion volume-adjusted                               #
# usage: python plot_summary_lollipop.py --summaryDir DIR --out FILE    #
#                                                                       #
#                                                                       #
# Author: Patrik Bey, patrik.bey@ucl.ac.uk                              #
#                                                                       #
#                                                                       #
# last update: 2026/11/10.                                              #
#                                                                       #
#                                                                       #
#########################################################################


import argparse
import glob
import os

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

UNADJ = "#4aa3df"
ADJ = "#2e8b57"
ALPHA_OUT, ALPHA_IN = 0.20, 0.45     # shading opacity: outward (H) / inward (~H)
GAP = 8.0                            # radius at which a spoke at the inward limit ends

# (Compartment, quadrant start angle in degrees, counter-clockwise from east)
ROWS = [[("SM", 0), ("MDN", 90), ("VAN", 180), ("DAN", 270)],
        [("AF", 0), ("SLF", 90), ("ILF", 180), ("IFOF", 270)]]
# (generator, outward limit %, inward limit %); None: generator == substrate
COLUMNS = [("global", 65, 20), (None, 105, 20)]


def load(summary_dir):
    files = sorted(glob.glob(os.path.join(summary_dir, "*_summary.csv")))
    if not files:
        raise SystemExit(f"no *_summary.csv files in {summary_dir}")
    return pd.concat([pd.read_csv(f) for f in files], ignore_index=True)


def spoke_angle(quadrant_start, k, n_bins):
    """Angle (radians) of noise bin k, running clockwise inside the quadrant."""
    return np.deg2rad(quadrant_start + 90 - 90 * (k + 0.5) / n_bins)


def shade(ax, ang, rU, rA, r0, sign, alpha, n=240):
    order = np.argsort(ang)
    a = np.asarray(ang)[order]
    dU = np.abs(np.asarray(rU)[order] - r0)
    dA = np.abs(np.asarray(rA)[order] - r0)
    th = np.linspace(a[0], a[-1], n)
    dU, dA = np.interp(th, a, dU), np.interp(th, a, dA)
    lo, hi = np.minimum(dU, dA), np.maximum(dU, dA)
    adj_lower, unadj_higher = dA <= dU, dU >= dA
    kw = dict(alpha=alpha, linewidth=0, interpolate=True, zorder=2)
    y0, ylo, yhi = np.full(n, r0), r0 + sign * lo, r0 + sign * hi
    ax.fill_between(th, y0, ylo, where=adj_lower, color=ADJ, **kw)
    ax.fill_between(th, y0, ylo, where=~adj_lower, color=UNADJ, **kw)
    ax.fill_between(th, ylo, yhi, where=unadj_higher, color=UNADJ, **kw)
    ax.fill_between(th, ylo, yhi, where=~unadj_higher, color=ADJ, **kw)
    ax.plot(th, r0 + sign * dU, color=UNADJ, linewidth=0.8, zorder=3)
    ax.plot(th, r0 + sign * dA, color=ADJ, linewidth=0.8, zorder=3)


def draw(ax, df, compartments, generator, rmax, rin):
    r0 = rin + GAP
    ax.set_ylim(0, r0 + rmax)
    ticks = list(range(-int(rin), 0, 10)) + list(range(0, int(rmax) + 1, 20))
    ax.set_yticks([r0 + v for v in ticks])
    ax.set_yticklabels([])
    ax.set_xticks([])
    ax.grid(color="0.88", linewidth=0.5)
    ax.spines["polar"].set_visible(False)
    th = np.linspace(0, 2 * np.pi, 400)
    ax.plot(th, np.full_like(th, r0), color="k", linestyle="--", linewidth=0.9, zorder=2)
    for name, start in compartments:
        sub = df[(df.hypothesized == name) & (df.behaviourGenerator == (generator or name))].sort_values("noiseBinStart")
        n = len(sub)
        tips = {key: ([], []) for key in ("U+", "U-", "A+", "A-")}
        for k, row in enumerate(sub.itertuples()):
            a = spoke_angle(start, k, n)
            tip = r0 + max(row.pctTargetFavouring_unadjusted, row.pctTargetFavouring_adjusted)
            ax.plot([a, a], [tip, r0 + 1.055 * rmax], color="0.55", linestyle=":", linewidth=1.2, zorder=1.5)
            for out, inn, colour in ((row.pctTargetFavouring_unadjusted, row.pctNonTargetFavouring_unadjusted, UNADJ),
                                     (row.pctTargetFavouring_adjusted, row.pctNonTargetFavouring_adjusted, ADJ)):
                for value, sign in ((out, 1), (inn, -1)):
                    r = r0 + sign * value
                    ax.plot([a, a], [r0, r], color=colour, linewidth=1.0, zorder=3)
                    ax.plot(a, r, "o" if colour == UNADJ else "^", color=colour,
                            markersize=9 if colour == UNADJ else 10.5, zorder=4)
                    key = ("U" if colour == UNADJ else "A") + ("+" if sign > 0 else "-")
                    tips[key][0].append(a)
                    tips[key][1].append(r)
            ax.text(a, r0 + 1.12 * rmax, f"{row.noiseBinStart:.2f}-{row.noiseBinEnd:.2f}", fontsize=12,
                    rotation=np.rad2deg(a) - 90 if np.cos(a) >= 0 else np.rad2deg(a) + 90, ha="center", va="center")
        shade(ax, tips["U+"][0], tips["U+"][1], tips["A+"][1], r0, +1, ALPHA_OUT)
        shade(ax, tips["U-"][0], tips["U-"][1], tips["A-"][1], r0, -1, ALPHA_IN)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--summaryDir", required=True)
    ap.add_argument("--out", default="Figure5.pdf")
    a = ap.parse_args()
    df = load(a.summaryDir)
    fig, axes = plt.subplots(2, 2, subplot_kw={"projection": "polar"}, figsize=(18, 18))
    for r, compartments in enumerate(ROWS):
        for c, (generator, rmax, rin) in enumerate(COLUMNS):
            draw(axes[r, c], df, compartments, generator, rmax, rin)
    fig.subplots_adjust(left=0.03, right=0.97, top=0.97, bottom=0.03, wspace=0.12, hspace=0.12)
    fig.savefig(a.out, dpi=150)
    print("saved", a.out)


if __name__ == "__main__":
    main()
