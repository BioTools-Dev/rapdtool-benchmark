#!/usr/bin/env python3
"""
plot_f1_threshold.py — F1 vs a UNIFORM abundance threshold, on the FOCUS/competitor
ABUNDANCE profiles (species rank).

CONTEXT: this figure is about the FOCUS *abundance* profile, not RaPDTool's detection.
RaPDTool's species detection is the mash confidence table, which already gives F1 = 1.0
with no threshold (see mash_detection.py). This figure shows that even the FOCUS
abundance profile, which carries a low-abundance tail, overtakes MetaPhlAn once a
relative-abundance threshold is applied uniformly to every tool's output — a caveat for
reading FOCUS composition, not the headline detection result. Title reflects that.

Reads the per-dataset CAMI profiles (OPAL-consistent taxonomy) directly.
"""
import collections
import os

INK, MUTED, GRID = "#1a1a1a", "#6b6b6b", "#e4e4e4"
COLORS = {"RaPDTool_screen": "#2a78d6", "RaPDTool_full": "#1baf7a",
          "Kraken2_full": "#e87ba4", "MetaPhlAn4": "#4a3aa7"}
LABEL = {"RaPDTool_screen": "RaPDTool screen", "RaPDTool_full": "RaPDTool full",
         "Kraken2_full": "Kraken2 full", "MetaPhlAn4": "MetaPhlAn4"}
FILES = {"RaPDTool_screen": "rapdtool_screen", "RaPDTool_full": "rapdtool_full",
         "Kraken2_full": "kraken2_full", "MetaPhlAn4": "metaphlan"}
THRESHOLDS = [0, 0.05, 0.1, 0.25, 0.5, 0.75, 1.0, 1.5, 2.0]
DATASET = "ln_30M"


def cami_species(path):
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


def main():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    B = "results/bench_%s" % DATASET
    truth = set(cami_species("%s/gold_standard.profile" % B))
    fig, ax = plt.subplots(figsize=(8.4, 5.2))

    for t in ["RaPDTool_screen", "RaPDTool_full", "MetaPhlAn4", "Kraken2_full"]:
        obs = cami_species("%s/profiles/%s.profile" % (B, FILES[t]))
        f1s = []
        for thr in THRESHOLDS:
            kept = {x for x, a in obs.items() if a >= thr}
            tp = len(kept & truth); fp = len(kept - truth); fn = len(truth - kept)
            f1s.append(2 * tp / (2 * tp + fp + fn) if (2 * tp + fp + fn) else 0)
        ax.plot(THRESHOLDS, f1s, "-o", color=COLORS[t], lw=2, ms=5,
                label=LABEL[t], zorder=3)

    ax.axvline(0.5, color="#9db8c6", lw=1, ls="--", zorder=1)
    ax.text(0.5, 1.03, "0.5 % cutoff", fontsize=8, color=MUTED, ha="center")
    ax.set_xlabel("uniform relative-abundance threshold (%) applied to every tool",
                  fontsize=9.5, color=MUTED)
    ax.set_ylabel("F1 score (species)", fontsize=9.5, color=MUTED)
    ax.set_title("Even the FOCUS abundance profile overtakes MetaPhlAn on F1 under a "
                 "shared cutoff\n(detection F1 is already 1.0 from the mash table)",
                 fontsize=11, color=INK, loc="left", pad=10, fontweight="bold")
    ax.set_ylim(-0.03, 1.08)
    ax.set_xlim(-0.05, 2.05)
    ax.grid(True, color=GRID, lw=0.7); ax.set_axisbelow(True)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    for s in ("left", "bottom"):
        ax.spines[s].set_color(GRID)
    ax.tick_params(colors=MUTED, labelsize=8.5, length=0)
    ax.legend(fontsize=8.5, frameon=False, loc="lower right")
    fig.text(0.012, 0.008,
             "mock_ln_30M, species rank. MetaPhlAn wins only at 0 % (its output is "
             "already marker-filtered); at any shared threshold ≥0.1 % RaPDTool leads, "
             "reaching F1 = 1.0 at 0.5 %. Above 1 % recall falls: two genomes are <1 % "
             "abundance by design.",
             fontsize=7.6, color=MUTED)
    fig.tight_layout(rect=[0, 0.05, 1, 1])

    os.makedirs("figures", exist_ok=True)
    for ext in ("svg", "png", "pdf"):
        fig.savefig("figures/f1_threshold.%s" % ext, dpi=300, bbox_inches="tight",
                    facecolor="white")
    print("[plot] wrote figures/f1_threshold.{svg,png,pdf}")


if __name__ == "__main__":
    main()
