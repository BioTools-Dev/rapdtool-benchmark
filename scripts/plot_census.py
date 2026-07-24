#!/usr/bin/env python3
"""
plot_census.py — figure for the type-material database census.

Renders census_full.tsv (produced by check_representation.py) as a two-panel figure:

  A. Overall: what fraction of RaPDTool's 30,209 type-material genomes is absent
     from each competitor database.
  B. Per phylum: fraction absent from BOTH Kraken2 standard and MetaPhlAn4, which
     is the population the 10 conflictive mock genomes were sampled from.

Design notes: this encodes a single measure (percent absent), so it uses one hue at
varying emphasis rather than a categorical palette -- no color-coded identity to
decode, and nothing that breaks under colorblind or greyscale print. Phyla that
contributed a genome to the mock are marked with a filled dot next to the label, so
that information is never carried by color alone.

Usage:
  ./plot_census.py -i census_full.tsv -l mock_genomes.list -o figures/census

Outputs PNG (300 dpi), SVG and PDF for every figure.
"""
import argparse
import collections
import csv
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

INK = "#1a1a1a"
MUTED = "#6b6b6b"
BAR = "#2a6f8f"          # single hue: magnitude, not identity
BAR_SOFT = "#a8c8d8"     # same hue, recessive step
GRID = "#e4e4e4"


def yes(row, key):
    return row[key].strip().upper().startswith("Y")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("-i", "--input", default="data/census_full.tsv")
    ap.add_argument("-l", "--list", default="data/mock_genomes.list")
    ap.add_argument("-o", "--output", default="figures/census")
    ap.add_argument("--min-genomes", type=int, default=50,
                    help="minimum genomes for a phylum to appear in panel B")
    args = ap.parse_args()

    rows = list(csv.DictReader(open(args.input), delimiter="\t"))
    n = len(rows)

    per = collections.defaultdict(lambda: [0, 0])
    for r in rows:
        p = r["phylum"] or "unclassified"
        per[p][1] += 1
        if not yes(r, "kraken_full") and not yes(r, "metaphlan"):
            per[p][0] += 1
    known = sorted(per, key=len, reverse=True)   # longest first, so submatches lose

    # Phyla represented among the mock's conflictive half, for the annotation. Match
    # against the census phylum names rather than parsing position: some list entries
    # carry a trailing parenthetical (e.g. "Archaea (Methanobacteriota)").
    sampled = set()
    conf_section = False
    for raw in open(args.list):
        if raw.lstrip().startswith("#"):
            if raw.lstrip().startswith("# ---") and "CONFLICTIVE" in raw.upper():
                conf_section = True
            continue
        if conf_section and raw.strip() and "#" in raw:
            comment = raw.split("#", 1)[1]
            for p in known:
                if p in comment:
                    sampled.add(p)
                    break

    both = sum(1 for r in rows if not yes(r, "kraken_full") and not yes(r, "metaphlan"))
    overall = [
        ("Kraken2 standard\n(103.7 GB)", 100 * sum(1 for r in rows if not yes(r, "kraken_full")) / n),
        ("Kraken2 capped-8\n(8.1 GB)", 100 * sum(1 for r in rows if not yes(r, "kraken8")) / n),
        ("MetaPhlAn4\n(59.8 GB)", 100 * sum(1 for r in rows if not yes(r, "metaphlan")) / n),
        ("Both Kraken2 std\n+ MetaPhlAn4", 100 * both / n),
    ]

    phyla = sorted(((c / t * 100, p, c, t) for p, (c, t) in per.items()
                    if t >= args.min_genomes), reverse=False)

    fig = plt.figure(figsize=(11, max(5.5, 0.34 * len(phyla) + 2.2)))
    # Panel A holds 4 bars against B's ~19; giving it the full column height would
    # stretch its bars absurdly. Confine A to the top of the left column and use the
    # space beneath it for the sampling note.
    gs = fig.add_gridspec(2, 2, width_ratios=[1, 1.35], height_ratios=[1, 1.05],
                          wspace=0.42, hspace=0.18)

    # --- Panel A ---------------------------------------------------------------
    ax = fig.add_subplot(gs[0, 0])
    labels = [l for l, _ in overall]
    vals = [v for _, v in overall]
    colors = [BAR_SOFT, BAR_SOFT, BAR_SOFT, BAR]  # emphasise the headline bar
    bars = ax.barh(range(len(vals)), vals, color=colors, height=0.6)
    ax.set_yticks(range(len(vals)))
    ax.set_yticklabels(labels, fontsize=9, color=INK)
    ax.invert_yaxis()
    for b, v in zip(bars, vals):
        ax.text(v + 1.2, b.get_y() + b.get_height() / 2, f"{v:.1f}%",
                va="center", fontsize=9.5, color=INK, fontweight="normal")
    ax.set_xlim(0, 72)
    ax.set_xlabel("Type-material genomes absent (%)", fontsize=9.5, color=MUTED)
    ax.set_title(f"A   Absent from competitor databases\n"
                 f"      n = {n:,} type-material genomes",
                 fontsize=10.5, color=INK, loc="left", pad=12)

    # --- Sampling note, under panel A ------------------------------------------
    axn = fig.add_subplot(gs[1, 0])
    axn.axis("off")
    axn.text(
        0, 1,
        f"The {both:,} genomes absent from both databases\n"
        f"({100 * both / n:.1f} % of all type material) are the\n"
        f"population from which the mock community's\n"
        f"10 conflictive genomes were drawn, stratified\n"
        f"by phylum (one per phylum) so that no single\n"
        f"lineage dominates the sample.\n\n"
        f"The 10 reference genomes were drawn from the\n"
        f"{sum(1 for r in rows if yes(r, 'kraken_full') and yes(r, 'metaphlan')):,} "
        f"genomes present in both, as a control.\n\n"
        f"Capping Kraken2 from 103.7 GB to 8.1 GB leaves\n"
        f"the absence rate essentially unchanged: capping\n"
        f"removes minimizers, not taxa.\n\n"
        f"Panel B shows phyla with ≥ {args.min_genomes} genomes; "
        f"{len(sampled - {p for _, p, _, _ in phyla})} sampled\n"
        f"phyla fall below that threshold and are not shown.",
        va="top", ha="left", fontsize=8.4, color=MUTED, linespacing=1.55,
        transform=axn.transAxes)

    # --- Panel B ---------------------------------------------------------------
    ax2 = fig.add_subplot(gs[:, 1])
    pv = [x[0] for x in phyla]
    pl = [x[1] for x in phyla]
    pc = [(x[2], x[3]) for x in phyla]
    b2 = ax2.barh(range(len(pv)), pv, color=BAR, height=0.62)
    ax2.set_yticks(range(len(pv)))
    ax2.set_yticklabels(
        [f"{'● ' if p in sampled else '   '}{p}" for p in pl],
        fontsize=8.5, color=INK)
    for b, v, (c, t) in zip(b2, pv, pc):
        ax2.text(v + 0.35, b.get_y() + b.get_height() / 2, f"{v:.1f}%  ({c:,}/{t:,})",
                 va="center", fontsize=7.8, color=MUTED)
    ax2.set_xlim(0, max(pv) * 1.42)
    ax2.set_xlabel("Absent from both databases (%)", fontsize=9.5, color=MUTED)
    ax2.set_title("B   The gap is taxonomically pervasive\n"
                  "      ● = phylum sampled for the mock community",
                  fontsize=10.5, color=INK, loc="left", pad=12)

    for a in (ax, ax2):
        a.xaxis.grid(True, color=GRID, linewidth=0.7)
        a.set_axisbelow(True)
        for s in ("top", "right", "left"):
            a.spines[s].set_visible(False)
        a.spines["bottom"].set_color(GRID)
        a.tick_params(axis="x", colors=MUTED, labelsize=8.5, length=0)
        a.tick_params(axis="y", length=0)

    fig.suptitle("Standard metagenomic databases omit most characterised type material",
                 fontsize=12.5, color=INK, x=0.055, ha="left", y=0.985,
                 fontweight="bold")

    os.makedirs(os.path.dirname(args.output) or ".", exist_ok=True)
    # PNG for drafts/preview, SVG + PDF as the vector masters journals ask for.
    for ext in ("png", "svg", "pdf"):
        fig.savefig(f"{args.output}.{ext}", dpi=300, bbox_inches="tight",
                    facecolor="white")
    print("[plot_census] wrote %s.{png,svg,pdf}" % args.output)
    print(f"[plot_census] panel B: {len(phyla)} phyla with >= {args.min_genomes} genomes")


if __name__ == "__main__":
    main()
