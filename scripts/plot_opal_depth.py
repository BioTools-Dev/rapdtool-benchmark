#!/usr/bin/env python3
"""
plot_opal_depth.py — profiling accuracy vs sequencing depth (OPAL, species rank).

Two panels over the uneven depth series (3/10/30 M reads, same fixed composition):
  A. Bray-Curtis distance to the gold standard (lower = better)
  B. F1 score (higher = better)
Depth is a sequential quantity, so tools are the categories; the equal-coverage control
(30 M, even) is drawn as a separate marker at the 30 M position for reference.

Reads each dataset's OPAL results.tsv directly (bench_<dataset>/opal/results.tsv), so it
is self-contained — no intermediate file to produce first.
"""
import csv
import os

INK, MUTED, GRID = "#1a1a1a", "#6b6b6b", "#e4e4e4"
# Categorical palette (validated): blue, green, magenta, yellow, aqua, orange.
COLORS = {
    "RaPDTool_screen": "#2a78d6", "RaPDTool_full": "#1baf7a",
    "Kraken2_full": "#e87ba4", "Kraken2_16GB": "#eda100",
    "Kraken2_8GB": "#eb6834", "MetaPhlAn4": "#4a3aa7",
}
LABEL = {"RaPDTool_screen": "RaPDTool screen", "RaPDTool_full": "RaPDTool full",
         "Kraken2_full": "Kraken2 full", "Kraken2_16GB": "Kraken2 16 GB",
         "Kraken2_8GB": "Kraken2 8 GB", "MetaPhlAn4": "MetaPhlAn4"}
DEPTHS = [("ln_3M", 3), ("ln_10M", 10), ("ln_30M", 30)]
ORDER = ["RaPDTool_screen", "RaPDTool_full", "MetaPhlAn4",
         "Kraken2_full", "Kraken2_16GB", "Kraken2_8GB"]


def main():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D

    # data[dataset][tool][metric] straight from each OPAL results.tsv (species rank).
    def load(dataset):
        d = {}
        for r in csv.DictReader(open("results/bench_%s/opal/results.tsv" % dataset), delimiter="\t"):
            if r["rank"] != "species":
                continue
            try:
                d.setdefault(r["tool"], {})[r["metric"]] = float(r["value"])
            except (ValueError, KeyError):
                pass
        return d
    data = {dk: load(dk) for dk, _ in DEPTHS}
    data["equalcov_30M"] = load("equalcov_30M")

    fig, (axA, axB) = plt.subplots(1, 2, figsize=(11, 4.9))

    panels = [(axA, "Bray-Curtis distance", "Bray-Curtis distance", "lower is better", False),
              (axB, "F1 score", "F1 score (species)", "higher is better", True)]
    x = [d for _, d in DEPTHS]
    for ax, metric, ylab, note, higher in panels:
        for t in ORDER:
            y = [data[dk][t][metric] for dk, _ in DEPTHS]
            ax.plot(x, y, "-o", color=COLORS[t], lw=1.8, ms=5, zorder=3)
            # equal-coverage control as a hollow marker at x=30
            yc = data["equalcov_30M"][t][metric]
            ax.plot([30], [yc], marker="D", color=COLORS[t], ms=6, zorder=4,
                    markerfacecolor="white", markeredgewidth=1.4)
        import matplotlib.ticker as mticker
        ax.set_xscale("log")
        ax.set_xticks(x); ax.set_xticklabels(["3 M", "10 M", "30 M"])
        ax.xaxis.set_minor_locator(mticker.NullLocator())   # no 4×10⁰ clutter
        ax.set_xlim(2.6, 34)
        ax.set_xlabel("sequencing depth (read pairs ×2)", fontsize=9.5, color=MUTED)
        ax.set_ylabel(ylab, fontsize=9.5, color=MUTED)
        ax.set_title("%s  ·  %s" % (metric, note), fontsize=10.5, color=INK,
                     loc="left", pad=8)
        ax.grid(True, color=GRID, lw=0.7); ax.set_axisbelow(True)
        for s in ("top", "right"):
            ax.spines[s].set_visible(False)
        for s in ("left", "bottom"):
            ax.spines[s].set_color(GRID)
        ax.tick_params(colors=MUTED, labelsize=8.5, length=0)

    handles = [Line2D([0], [0], color=COLORS[t], lw=2.4, marker="o", ms=5,
                      label=LABEL[t]) for t in ORDER]
    handles.append(Line2D([0], [0], color=MUTED, lw=0, marker="D", ms=6,
                          markerfacecolor="white", markeredgewidth=1.4,
                          label="equal-coverage control (30 M)"))
    axB.legend(handles=handles, fontsize=7.6, frameon=False, loc="center right")

    fig.suptitle("Profiling accuracy vs sequencing depth (OPAL, species rank)",
                 fontsize=12.5, color=INK, x=0.045, ha="left", y=0.99, fontweight="bold")
    fig.text(0.045, 0.005,
             "Uneven community, fixed composition; depth is the only variable. RaPDTool "
             "screen tracks the gold standard closely (Bray-Curtis ~0.05); MetaPhlAn wins "
             "F1 through precision while recovering half the species. Unfiltered profiles.",
             fontsize=7.8, color=MUTED)
    fig.tight_layout(rect=[0, 0.04, 1, 0.96])

    os.makedirs("figures", exist_ok=True)
    for ext in ("svg", "png", "pdf"):
        fig.savefig("figures/opal_depth.%s" % ext, dpi=300, bbox_inches="tight",
                    facecolor="white")
    print("[plot] wrote figures/opal_depth.{svg,png,pdf}")


if __name__ == "__main__":
    main()
