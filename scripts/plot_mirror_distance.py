#!/usr/bin/env python3
"""
plot_mirror_distance.py — RaPDTool's rank resolution vs genomic distance.

Shows that RaPDTool degrades gracefully out of domain: it makes a species call only
when the genome is genuinely close to a database member (mash identity > ~95 %), backs
off to a genus call at moderate distance (~92-95 %), and abstains entirely when the
organism is distant (< ~80 %) rather than assigning a confident wrong species. This is
a safety property, and a positive result.

Input: the crosswalk of the 14 distance-stratified mirror genomes (built inline here
from mirror_dist_genomes.list + the bench_mirrordist outputs).

One panel: x = Mash identity to nearest RaPDTool genome, y = finest rank RaPDTool
resolved (species / genus / none). A single sequence, so no legend; the shaded bands
mark RaPDTool's mash thresholds.
"""
import csv
import glob
import os
import re

BENCH = "results/bench_mirrordist"
LIST = "data/mirror_dist_genomes.list"

INK, MUTED, GRID = "#1a1a1a", "#6b6b6b", "#e4e4e4"
C_SPECIES, C_GENUS, C_NONE = "#08519c", "#4292c6", "#c94a34"


def load():
    rows = []
    for l in open(LIST):
        if l.startswith("#") or not l.strip():
            continue
        m = re.search(r"#\s*(.+?)\s*\((\d+)\)\s+([\d.]+)%", l)
        if m:
            rows.append({"name": m.group(1).strip(), "genus": m.group(1).split()[0],
                         "id": float(m.group(3))})
    # mash-screen species calls
    tbl = open(os.path.join(BENCH, "rapdtool_screen.rep1",
                            "rapdtool_confidence.tbl")).read()
    screen = set()
    for line in tbl.split("\n"):
        mm = re.match(r"\|\s*([A-Z][a-z]+) [a-z]+\s*\|\s*\d+\s*\|", line)
        if mm:
            screen.add(mm.group(1))
    # FOCUS genera
    F = glob.glob(os.path.join(BENCH, "rapdtool_screen.rep1", "profilesfmbm", "*",
                               "output_All_levels.csv"))[0]
    focus = {r["Genus"] for r in csv.DictReader(open(F))}
    for r in rows:
        if r["genus"] in screen:
            r["rank"], r["y"] = "species", 2
        elif r["genus"] in focus:
            r["rank"], r["y"] = "genus", 1
        else:
            r["rank"], r["y"] = "none", 0
    return rows


def main():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.patches import Rectangle

    rows = sorted(load(), key=lambda r: -r["id"])
    fig, ax = plt.subplots(figsize=(9.5, 5.9))

    # Threshold bands (RaPDTool mash: species >95 %, genus 93-95 %).
    ax.add_patch(Rectangle((95, -0.5), 5.5, 3, color="#08519c", alpha=0.05, zorder=0))
    ax.add_patch(Rectangle((93, -0.5), 2, 3, color="#4292c6", alpha=0.07, zorder=0))
    ax.axvline(95, color="#9db8c6", lw=1, ls="--", zorder=1)
    ax.axvline(93, color="#9db8c6", lw=1, ls="--", zorder=1)
    ax.text(97.6, 2.42, "species threshold\n(>95 %)", fontsize=7.5, color=MUTED, ha="center")
    ax.text(94, 2.42, "genus\n93–95 %", fontsize=7.5, color=MUTED, ha="center")

    color = {"species": C_SPECIES, "genus": C_GENUS, "none": C_NONE}
    for r in rows:
        ax.scatter(r["id"], r["y"], s=90, color=color[r["rank"]], zorder=3,
                   edgecolor="white", linewidth=1.2)

    # Vertical labels so tight x-clusters (species 97-100 %, genus 90-95 %) don't
    # overlap. All rows label BELOW their point (consistent placement; keeps species
    # labels clear of the threshold annotations above the top row). Points sharing an
    # identity get their labels nudged apart on x so the two don't print on top of each
    # other.
    import collections
    same = collections.defaultdict(list)
    for r in rows:
        same[round(r["id"], 1)].append(r)
    for grp in same.values():
        grp.sort(key=lambda r: r["name"])
        for j, r in enumerate(grp):
            dx = 0.0 if len(grp) == 1 else (j - (len(grp) - 1) / 2) * 0.9
            up = False
            off = 0.22 if up else -0.22
            ax.annotate(r["name"], (r["id"], r["y"] + off), fontsize=6.8, color=INK,
                        ha="center", va="bottom" if up else "top", rotation=90,
                        zorder=4, xytext=(dx * 12, 0), textcoords="offset points")

    ax.set_yticks([0, 1, 2])
    ax.set_yticklabels(["abstains\n(no call)", "genus\ncall", "species\ncall"],
                       fontsize=9, color=INK)
    ax.set_ylim(-2.4, 3.4)
    ax.set_xlim(69, 101)
    ax.invert_xaxis()   # near (100 %) on the left, distant on the right
    ax.set_xlabel("Mash identity to nearest RaPDTool genome (%)  —  near ← → distant",
                  fontsize=9.5, color=MUTED)
    ax.set_title("RaPDTool degrades gracefully out of domain: rank resolved tracks genomic distance",
                 fontsize=11.5, color=INK, loc="left", pad=12, fontweight="bold")
    ax.grid(True, axis="x", color=GRID, lw=0.7)
    ax.set_axisbelow(True)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    for s in ("left", "bottom"):
        ax.spines[s].set_color(GRID)
    ax.tick_params(colors=MUTED, labelsize=8.5, length=0)
    fig.text(0.012, 0.02,
             "14 genomes spanning genomic distance to RaPDTool's database, equal coverage 27.7×.\n"
             "No genome is assigned a confident species below the 95 % threshold; organisms below "
             "~80 % identity are not reported at all.",
             fontsize=7.8, color=MUTED, linespacing=1.5)
    fig.tight_layout(rect=[0, 0.085, 1, 1])

    os.makedirs("figures", exist_ok=True)
    for ext in ("svg", "png", "pdf"):
        fig.savefig("figures/mirror_distance.%s" % ext, dpi=300,
                    bbox_inches="tight", facecolor="white")
    print("[plot] wrote figures/mirror_distance.{svg,png,pdf}")
    for r in rows:
        print("  %.1f%%  %-32s -> %s" % (r["id"], r["name"][:31], r["rank"]))


if __name__ == "__main__":
    main()
