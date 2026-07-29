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

One panel: x = Mash identity to nearest RaPDTool genome (measured independently), y =
finest rank RaPDTool resolved (species / genus / none), read from FULL mode so that the
classification uses the cutoffs rapdtool_results.pl actually applies (distance < 0.05
species, 0.05-0.08 genus). The shaded bands are those two cutoffs; a hollow mark is a
genus that only the FOCUS profile placed.
"""
import csv
import glob
import os
import re

BENCH = "results/bench_mirrordist"
LIST = "data/mirror_dist_genomes.list"

INK, MUTED, GRID = "#1a1a1a", "#6b6b6b", "#e4e4e4"
C_SPECIES, C_GENUS, C_NONE = "#08519c", "#4292c6", "#c94a34"


def cells(block):
    """First column of every data row of a Text::SimpleTable block."""
    out = []
    for line in block.split("\n"):
        if not line.startswith("|"):
            continue
        c = line.split("|")[1].strip()
        if not c or "closest-hit" in c or c in ("Species", "Genus"):
            continue                      # header row
        out.append(c)
    return out


def load():
    rows = []
    for l in open(LIST):
        if l.startswith("#") or not l.strip():
            continue
        m = re.search(r"#\s*(.+?)\s*\((\d+)\)\s+([\d.]+)%", l)
        if m:
            rows.append({"name": m.group(1).strip(), "genus": m.group(1).split()[0],
                         "id": float(m.group(3))})
    # The rank each genome is resolved to is read from FULL mode, because that is the
    # path whose cutoffs are documented and fixed: rapdtool_results.pl sorts a bin by
    # its Mash DISTANCE to the nearest reference — < 0.05 into "Species with high
    # confidence" (> 95 % identity) and 0.05-0.08 into "Genus with high confidence"
    # (92-95 %). Screen mode applies a single --screen-identity cutoff (0.95) with no
    # genus tier, so classifying this experiment from screen reports every genus call as
    # FOCUS-derived and leaves the tool's own genus threshold untested.
    # The x axis stays the independently measured minimum Mash distance to the database
    # (data/mirror_distance.tsv); only the classification comes from this table.
    tbl = open(os.path.join(BENCH, "rapdtool.rep1",
                            "rapdtool_confidence.tbl")).read()
    gn_block = tbl.split("Genus with high confidence")[1].split("Species with high")[0]
    sp_block = tbl.split("Species with high confidence")[1].split("FOCUS profile")[0]
    fo_block = tbl.split("FOCUS profile")[1]
    mash_species = {c.split()[0] for c in cells(sp_block)}
    mash_genus = {c.split()[0] for c in cells(gn_block)}
    focus_genus = {c.split("_")[0] for c in cells(fo_block)}

    for r in rows:
        g = r["genus"]
        if g in mash_species:
            r["rank"], r["y"], r["by"] = "species", 2, "mash"
        elif g in mash_genus:
            r["rank"], r["y"], r["by"] = "genus", 1, "mash"
        elif g in focus_genus:
            r["rank"], r["y"], r["by"] = "genus", 1, "focus"
        else:
            r["rank"], r["y"], r["by"] = "none", 0, "none"
    return rows


def main():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.patches import Rectangle

    rows = sorted(load(), key=lambda r: -r["id"])
    fig, ax = plt.subplots(figsize=(11, 6.2))

    # The two bands are the cutoffs rapdtool_results.pl actually applies to a bin's
    # Mash distance: < 0.05 -> species, 0.05-0.08 -> genus. Expressed as identity that
    # is > 95 % and 92-95 %, so the lower edge is 92, not 93.
    ax.add_patch(Rectangle((95, -0.5), 7.6, 3, color="#08519c", alpha=0.05, zorder=0))
    ax.add_patch(Rectangle((92, -0.5), 3, 3, color="#4292c6", alpha=0.07, zorder=0))
    ax.axvline(95, color="#9db8c6", lw=1, ls="--", zorder=1)
    ax.axvline(92, color="#9db8c6", lw=1, ls="--", zorder=1)
    ax.text(99.4, 2.42, "species: Mash distance < 0.05\n(identity > 95 %)",
            fontsize=7.5, color=MUTED, ha="center")
    ax.text(93.5, 2.42, "genus: 0.05–0.08\n(92–95 %)", fontsize=7.5, color=MUTED,
            ha="center")

    # Genomes sharing an identity would otherwise be drawn exactly on top of one
    # another — two pairs do (100.0 % species, 73.7 % abstain), so the figure showed
    # 12 marks for 14 genomes and could not be counted against the rank table. A tie
    # is therefore dodged on the RANK axis, which is categorical; the measured
    # identity, carried by x, is never moved.
    import collections
    same = collections.defaultdict(list)
    for r in rows:
        same[round(r["id"], 1)].append(r)
    for grp in same.values():
        grp.sort(key=lambda r: r["name"])
        n = len(grp)
        for j, r in enumerate(grp):
            r["dy"] = 0.0 if n == 1 else (j - (n - 1) / 2) * 0.34
            r["tied"] = n > 1
        # Labels in this column must clear the LOWEST mark of the column, not each
        # label's own mark, or the dodged partner is drawn on top of the text.
        floor = min(x["dy"] for x in grp)
        for r in grp:
            r["dyfloor"] = floor

    # Hue carries the RANK; fill carries WHICH ENGINE resolved it. A hollow mark is a
    # genus that only the FOCUS profile placed, because the bin sat beyond mash's 0.08
    # genus cutoff — the distinction the shaded genus band exists to make visible.
    color = {"species": C_SPECIES, "genus": C_GENUS, "none": C_NONE}
    for r in rows:
        c = color[r["rank"]]
        solid = r["by"] != "focus"
        ax.scatter(r["id"], r["y"] + r["dy"], s=90, zorder=3, linewidth=1.6,
                   color=c if solid else "white",
                   edgecolor="white" if solid else c)

    # Names are set vertically, so each one occupies a narrow column at its own
    # identity and neighbours closer than MIN_GAP would print into each other — the
    # species cluster (100, 100, 99.5, 98.3, 97.0) and both tied pairs do. Rather
    # than move a label away from its mark on x, which detaches it from the point it
    # names, a crowded label is pushed DOWN to the first free depth and joined to its
    # mark by a leader line. Depth is free space here; horizontal room is not.
    MIN_GAP = 0.45          # identity units needed between two labels at one depth
    DEPTH = (0.20, 2.25)    # spacing exceeds the height of a rotated name (~1.8 units)
    LABEL_DX = 0.95         # sideways room given to the deeper label of a tied pair
    for rank in ("species", "genus", "none"):
        grp = sorted([r for r in rows if r["rank"] == rank], key=lambda r: -r["id"])
        ids = [r["id"] for r in grp]
        last = {}
        for r in grp:
            lvl = 0
            while lvl in last and abs(last[lvl] - r["id"]) < MIN_GAP:
                lvl += 1
            last[lvl] = r["id"]
            r["lvl"] = min(lvl, len(DEPTH) - 1)
            # Depth alone leaves the two names of a tied pair stacked in one column,
            # reading as a single block of text. The deeper one is therefore also moved
            # sideways, AWAY from the rest of the row so it lands in empty axis space —
            # toward 100 % for a tie at the near end, toward the far end for one at the
            # other. A tie shares both identity and rank, so displacing its label
            # asserts nothing about the measurement.
            out = 1.0 if (r["id"] - min(ids)) > (max(ids) - r["id"]) else -1.0
            r["ldx"] = out * LABEL_DX if (r["tied"] and r["lvl"]) else 0.0

    for r in rows:
        lx = r["id"] + r["ldx"]
        ytop = r["y"] + r["dyfloor"] - DEPTH[r["lvl"]]
        # A leader line is drawn only where the pairing actually carries information —
        # i.e. for a crowded but DISTINCT identity. Both members of a tie share an
        # identity and a rank, so which name belongs to which of the two marks says
        # nothing, and a line there would only cross the other label's text.
        if r["lvl"] and not r["tied"]:
            ax.plot([r["id"], r["id"]], [r["y"] + r["dy"] - 0.10, ytop + 0.03],
                    color="#c2c2c2", lw=0.6, zorder=2, solid_capstyle="round")
        ax.annotate(r["name"], (lx, ytop), fontsize=6.8, color=INK,
                    ha="center", va="top", rotation=90, zorder=4)

    n_rank = collections.Counter(r["rank"] for r in rows)
    ax.set_yticks([0, 1, 2])
    ax.set_yticklabels(["abstains\n(no call)  n=%d" % n_rank["none"],
                        "genus\ncall  n=%d" % n_rank["genus"],
                        "species\ncall  n=%d" % n_rank["species"]],
                       fontsize=9, color=INK)
    ax.set_ylim(-4.1, 3.4)   # headroom for the staggered label depths below row 0
    ax.set_xlim(69, 102.6)   # headroom at the near end so the 100 % labels clear the y axis
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
             "~80 % identity are not reported at all.\n"
             "Two pairs share an identity (100.0 %, 73.7 %); each pair's marks are offset vertically "
             "and its names set side by side so both are legible — neither offset carries meaning.\n"
             "Rank is read from full mode, whose cutoffs are fixed: Mash distance < 0.05 -> species, "
             "0.05–0.08 -> genus. Mash resolves the genus of the three\ngenomes inside that window; "
             "beyond it three more are placed to genus by the FOCUS profile alone (hollow marks).",
             fontsize=7.8, color=MUTED, linespacing=1.5)
    from matplotlib.lines import Line2D
    ax.legend(handles=[
        Line2D([], [], marker="o", ls="", ms=7, mfc=C_GENUS, mec="white", mew=1.4,
               label="rank resolved by Mash distance (species or genus tier)"),
        Line2D([], [], marker="o", ls="", ms=7, mfc="white", mec=C_GENUS, mew=1.4,
               label="genus resolved by the FOCUS profile only (beyond the 0.08 cutoff)"),
    ], loc="lower left", bbox_to_anchor=(0.44, 0.015), frameon=False, fontsize=7.6,
        labelcolor=MUTED, handletextpad=0.6, borderpad=0.2)

    fig.tight_layout(rect=[0, 0.155, 1, 1])   # room for the 5-line footnote

    os.makedirs("figures", exist_ok=True)
    for ext in ("svg", "png", "pdf"):
        fig.savefig("figures/mirror_distance.%s" % ext, dpi=300,
                    bbox_inches="tight", facecolor="white")
    print("[plot] wrote figures/mirror_distance.{svg,png,pdf}")
    for r in rows:
        print("  %.1f%%  %-32s -> %s" % (r["id"], r["name"][:31], r["rank"]))


if __name__ == "__main__":
    main()
