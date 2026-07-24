#!/usr/bin/env python3
"""
plot_micomplete.py — recovered-bin quality: completeness vs contamination.

Two panels (30 M uneven, 30 M equal-coverage). Each point is one recovered genome bin;
the shaded corner is the high-quality region (>=90 % complete, <5 % contamination). The
two reproducible chimeric bins (high redundancy) are labelled — miComplete's redundancy
metric catches them, which is the point.

Reads bench_<dataset>/rapdtool.rep1/workfmbm/outmicomplete/miCompleteOut_*.tab directly.
"""
import glob
import os

INK, MUTED, GRID = "#1a1a1a", "#6b6b6b", "#e4e4e4"
GOOD, BAD = "#1baf7a", "#c94a34"
DATASETS = [("ln_30M", "30 M reads, uneven"), ("equalcov_30M", "30 M reads, equal coverage")]


def load(dataset):
    tab = glob.glob("results/bench_%s/rapdtool.rep1/workfmbm/outmicomplete/miCompleteOut_*.tab"
                    % dataset)
    if not tab:
        return []
    # Bin→species names from the species_bins/ directory, or from the shipped
    # species_bins.list (the bin FASTAs themselves are large and not distributed).
    d = "results/bench_%s/rapdtool.rep1/species_bins" % dataset
    lst = "results/bench_%s/species_bins.list" % dataset
    names = os.listdir(d) if os.path.isdir(d) else (
        [x.strip() for x in open(lst)] if os.path.exists(lst) else [])
    binmap = {}
    for b in names:
        if not b:
            continue
        num = b.split("final_contigs_")[-1].replace(".fna", "")
        binmap["final_contigs_%s" % num] = b.split("__")[0]
    out = []
    for l in open(tab[0]):
        if l.startswith(("#", "Name")):
            continue
        p = l.rstrip("\n").split("\t")
        if len(p) >= 6:
            try:
                out.append((binmap.get(p[0], p[0]), float(p[4]) * 100,
                            (float(p[5]) - 1) * 100))
            except ValueError:
                pass
    return out


def main():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.patches import Rectangle

    fig, axes = plt.subplots(1, 2, figsize=(11, 5.0), sharey=True)
    for ax, (ds, title) in zip(axes, DATASETS):
        rows = load(ds)
        ax.add_patch(Rectangle((90, -2), 15, 7, color=GOOD, alpha=0.06, zorder=0))
        for name, comp, cont in rows:
            chim = cont >= 10
            ax.scatter(comp, cont, s=70, color=BAD if chim else GOOD, zorder=3,
                       edgecolor="white", linewidth=1.1)
            if chim:
                ax.annotate("%s (%.0f%%)" % (name.replace("_", " ")[:22], cont),
                            (comp, cont), fontsize=6.6, color=INK,
                            xytext=(-6, 6), textcoords="offset points", ha="right")
        n = len(rows)
        hq = sum(1 for _, c, ct in rows if c >= 90 and ct < 5)
        ax.set_title("%s\n%d bins · %d high-quality" % (title, n, hq),
                     fontsize=10, color=INK, loc="left", pad=8)
        ax.set_xlabel("completeness (%)", fontsize=9.5, color=MUTED)
        ax.set_xlim(40, 103)
        ax.set_ylim(-4, 120)
        ax.grid(True, color=GRID, lw=0.7)
        ax.set_axisbelow(True)
        for s in ("top", "right"):
            ax.spines[s].set_visible(False)
        for s in ("left", "bottom"):
            ax.spines[s].set_color(GRID)
        ax.tick_params(colors=MUTED, labelsize=8.5, length=0)
        ax.axhline(5, color="#9db8c6", lw=1, ls="--", zorder=1)
    axes[0].set_ylabel("contamination = redundancy − 1 (%)", fontsize=9.5, color=MUTED)

    fig.suptitle("Recovered-bin quality (miComplete) — most bins high-quality, "
                 "chimeras flagged by redundancy",
                 fontsize=12, color=INK, x=0.045, ha="left", y=0.99, fontweight="bold")
    fig.text(0.045, 0.005,
             "Green = high quality (shaded: ≥90 % complete, <5 % contamination). Red = "
             "chimeric bins (≥10 % contamination), the same two organisms in both "
             "datasets — a reproducible binning limitation, correctly flagged by the "
             "redundancy metric. miComplete measures quality, not taxonomic correctness.",
             fontsize=7.6, color=MUTED)
    fig.tight_layout(rect=[0, 0.04, 1, 0.95])

    os.makedirs("figures", exist_ok=True)
    for ext in ("svg", "png", "pdf"):
        fig.savefig("figures/micomplete.%s" % ext, dpi=300, bbox_inches="tight",
                    facecolor="white")
    print("[plot] wrote figures/micomplete.{svg,png,pdf}")


if __name__ == "__main__":
    main()
