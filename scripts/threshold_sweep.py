#!/usr/bin/env python3
"""
threshold_sweep.py — choose RaPDTool's relative-abundance cutoff from evidence.

RaPDTool applies a 1 % cutoff in its user-facing confidence table. This sweeps the
cutoff over a FOCUS profile and, against a known truth set, counts true positives,
false positives, recall and precision at each level, for both RaPDTool modes:

    full   = FOCUS on the assembly   (genome-recovery mode)
    screen = FOCUS on the reads      (matched input vs Kraken2/MetaPhlAn)

Output: a TSV of the sweep (cite exact numbers) and a precision-recall trade-off
figure (SVG + PNG 300 dpi + PDF), one panel per mode. The elbow of each curve is the
cutoff that removes false positives without yet costing true positives.

METHOD NOTE — define on one dataset, validate on another. Optimising and reporting the
cutoff on the same data overfits. Define on the uneven depth series, then confirm the
elbow holds on the equal-coverage control and on ZymoBIOMICS (real data). A cutoff that
moves between them is not robust.

Depth is encoded as one blue hue at increasing intensity (it is an ordered quantity,
not a category); the equal-coverage control is a separate accent colour.

Truth is read from a CAMI gold-standard profile (same taxonomy resolution as OPAL), not
from mock_genomes.list — the four depth datasets share the same 20-taxon set, so any
one dataset's gold_standard.profile works.

Usage:
  ./threshold_sweep.py --truth bench_ln_30M/gold_standard.profile \\
      --datasets ln_3M,ln_10M,ln_30M,equalcov_30M --benchroot . -o figures/threshold_sweep
"""
import argparse
import collections
import os
import sys

THRESHOLDS = [0, 0.05, 0.1, 0.25, 0.5, 0.75, 1.0, 1.5, 2.0]
# screen first: it is the matched-input profiling mode and the one the cutoff is for.
MODES = [("rapdtool_screen", "screen"), ("rapdtool", "full")]

# Depth series: sequential blue, light -> dark. Control: distinct accent.
DEPTH_BLUE = {"ln_3M": "#9ecae1", "ln_10M": "#4292c6", "ln_30M": "#08519c"}
ACCENT = "#eb6834"
INK, MUTED, GRID = "#1a1a1a", "#6b6b6b", "#e4e4e4"


def cami_species(path):
    """{species_taxid: abundance_percent} from a CAMI/BIOBOXES profile.

    Read from the CAMI profile, NOT the raw FOCUS Species_tabular, so the taxonomy
    resolution matches OPAL exactly (both come through profile2cami.py). An earlier
    version resolved FOCUS names with a plain name2taxid, which disagreed with OPAL by
    ~2 taxa and made the recall ceiling look like 18/20 where OPAL sees 20/20."""
    obs = collections.defaultdict(float)
    for line in open(path):
        if line.startswith(("@", "#")) or not line.strip():
            continue
        f = line.rstrip("\n").split("\t")
        if len(f) >= 5 and f[1] == "species":
            try:
                obs[f[0]] += float(f[4])
            except ValueError:
                pass
    return obs


def truth_taxids(path):
    """Species taxids from a CAMI gold-standard profile (the same resolution as OPAL)."""
    return set(cami_species(path))


def sweep(obs, truth):
    rows = []
    for thr in THRESHOLDS:
        kept = {t for t, a in obs.items() if a >= thr}
        tp = len(kept & truth)
        fp = len(kept - truth)
        rec = tp / len(truth) if truth else 0.0
        prec = tp / (tp + fp) if kept else 0.0
        f1 = 2 * prec * rec / (prec + rec) if (prec + rec) else 0.0
        rows.append((thr, tp, fp, rec, prec, f1))
    return rows


def find_profile(benchroot, dataset, mode):
    p = os.path.join(benchroot, "bench_%s" % dataset, "profiles",
                     "rapdtool_%s.profile" % ("screen" if mode == "screen" else "full"))
    return p if os.path.exists(p) else None


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--truth", required=True,
                    help="CAMI gold-standard profile (e.g. bench_ln_30M/gold_standard.profile)")
    ap.add_argument("--datasets", default="ln_3M,ln_10M,ln_30M,equalcov_30M")
    ap.add_argument("--benchroot", default="results")
    ap.add_argument("-o", "--output", default="figures/threshold_sweep")
    args = ap.parse_args()

    truth = truth_taxids(args.truth)
    if not truth:
        sys.exit("ERROR: no truth taxids parsed from %s" % args.truth)
    datasets = [d.strip() for d in args.datasets.split(",") if d.strip()]
    print("[sweep] %d truth taxa; datasets: %s" % (len(truth), ", ".join(datasets)))

    # data[mode][dataset] = sweep rows
    data = {m: {} for _, m in MODES}
    tsv_rows = []
    for ds in datasets:
        for subdir, mode in MODES:
            prof = find_profile(args.benchroot, ds, mode)
            if not prof:
                print("[sweep] %-14s %-6s : (no CAMI profile, skipped)" % (ds, mode))
                continue
            obs = cami_species(prof)
            n_raw = n_res = len(obs)
            rows = sweep(obs, truth)
            data[mode][ds] = rows
            det = next(r[1] for r in rows if r[0] == 0)
            print("[sweep] %-14s %-6s : %d species raw, %d/%d truth detected"
                  % (ds, mode, n_raw, det, len(truth)))
            for thr, tp, fp, rec, prec, f1 in rows:
                tsv_rows.append((ds, mode, thr, tp, fp, rec, prec, f1))

    os.makedirs(os.path.dirname(args.output) or ".", exist_ok=True)
    with open(args.output + ".tsv", "w") as fh:
        fh.write("dataset\tmode\tthreshold\ttp\tfp\trecall\tprecision\tf1\n")
        for r in tsv_rows:
            fh.write("%s\t%s\t%.2f\t%d\t%d\t%.4f\t%.4f\t%.4f\n" % r)
    print("[sweep] wrote %s.tsv" % args.output)

    render(data, datasets, truth, args.output)


def color_for(ds):
    return DEPTH_BLUE.get(ds, ACCENT)


def render(data, datasets, truth, out):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(1, 2, figsize=(11, 5.2), sharex=True, sharey=True)
    for ax, (_, mode) in zip(axes, MODES):
        for ds in datasets:
            rows = data[mode].get(ds)
            if not rows:
                continue
            rec = [r[3] for r in rows]
            prec = [r[4] for r in rows]
            c = color_for(ds)
            ax.plot(rec, prec, "-", color=c, lw=1.8, zorder=2)
            ax.plot(rec, prec, "o", color=c, ms=4, zorder=3)
            # annotate the current default (1 %) and the apparent elbow (0.5 %)
            for thr_mark, dx, dy in ((1.0, 4, -10), (0.5, 5, 6)):
                for r in rows:
                    if abs(r[0] - thr_mark) < 1e-9 and r[4] > 0:
                        ax.annotate("%.2g%%" % thr_mark, (r[3], r[4]),
                                    textcoords="offset points", xytext=(dx, dy),
                                    fontsize=6.5, color=MUTED)
        ax.set_title("FOCUS on %s" % ("reads (screen mode)" if mode == "screen"
                                      else "assembly (full mode)"),
                     fontsize=10.5, color=INK, loc="left", pad=8)
        ax.set_xlabel("Recall  (true species detected / %d)" % len(truth),
                      fontsize=9.5, color=MUTED)
        ax.grid(True, color=GRID, lw=0.7)
        ax.set_axisbelow(True)
        for s in ("top", "right"):
            ax.spines[s].set_visible(False)
        for s in ("left", "bottom"):
            ax.spines[s].set_color(GRID)
        ax.tick_params(colors=MUTED, labelsize=8.5, length=0)
        ax.set_xlim(0.4, 1.02)
        ax.set_ylim(0, 1.05)
    axes[0].set_ylabel("Precision  (true / all reported)", fontsize=9.5, color=MUTED)

    # Legend: depth series (sequential) + control, plus the read on the curve.
    from matplotlib.lines import Line2D
    handles = []
    for ds in datasets:
        if data["screen"].get(ds) or data["full"].get(ds):
            lab = ds.replace("ln_", "").replace("_", " ")
            handles.append(Line2D([0], [0], color=color_for(ds), lw=2.4, marker="o",
                                  ms=4, label=lab))
    axes[1].legend(handles=handles, title="dataset", fontsize=8, title_fontsize=8.5,
                   frameon=False, loc="lower left")

    fig.suptitle("Abundance-cutoff trade-off: precision bought at the cost of recall",
                 fontsize=12.5, color=INK, x=0.045, ha="left", y=0.99, fontweight="bold")
    fig.text(0.045, 0.005,
             "Each point is one cutoff (0–2 %). The elbow — where the curve turns down — "
             "is the cutoff that clears false positives without dropping true species. "
             "Points annotated at 0.5 % and 1 %.",
             fontsize=8, color=MUTED)
    fig.tight_layout(rect=[0, 0.03, 1, 0.96])

    for ext in ("svg", "png", "pdf"):
        fig.savefig("%s.%s" % (out, ext), dpi=300, bbox_inches="tight", facecolor="white")
    print("[sweep] wrote %s.{svg,png,pdf}" % out)


if __name__ == "__main__":
    main()
