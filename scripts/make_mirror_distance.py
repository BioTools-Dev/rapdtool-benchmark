#!/usr/bin/env python3
"""
make_mirror_distance.py — build a DISTANCE-STRATIFIED mirror set.

The point of the mirror experiment is not "abstain vs misassign" but RANK RESOLUTION
as a function of genomic distance. RaPDTool's mash step calls a species above ~95 %
identity (distance < 0.05) and a genus above ~93 % (distance ~0.05-0.07); beyond that
it should back off to a higher rank rather than assign a confident wrong species. The
right test is therefore: sample organisms spanning a range of genomic distances to
RaPDTool's nearest genome, and check that the rank RaPDTool resolves each to tracks the
distance correctly — species only when truly close, genus at moderate distance, higher
or nothing when far. Graceful degradation is a POSITIVE result.

Selection frame (taxonomy is a proxy for distance; the real x-axis is the measured Mash
distance):
  tier A  genus present in RaPDTool   (novel species)        -> expect species/genus call
  tier B  genus absent, family present                       -> expect genus/family call
  tier C  family absent                                      -> expect family+/abstain
plus a couple of tier-0 organisms (a species RaPDTool HAS) as positive controls, to show
a species-level call is made when the organism is genuinely present.

For each candidate: download the genome (NCBI datasets), sketch it, and take the minimum
Mash distance to RaPDTool's mash DB (0726_30213genomes.msh) — the nearest-neighbour
distance. The output table, sorted by distance, is the sampling frame for the final
spread; pick ~2-3 per distance band for the mock.

Usage:
  ./make_mirror_distance.py --per-tier 8 --seed 42 -o mirror_distance.tsv
"""
import argparse
import collections
import os
import random
import subprocess
import sys

ACCMAP = os.environ.get("FOCUS_DB", "/path/to/focus_build") + "/acc_taxid_strain.tsv"
MASH_DB = os.environ.get("MASH_DB", "/path/to/0726_30213genomes.msh")
MASH = os.environ.get("MASH_BIN", "mash")
POOL = "data/mirror_pool.tsv"
FNA_DIR = os.environ.get("MIRROR_FNA_DIR", "/path/to/mirror_db/fna_dist")
TAXONKIT_DB = os.environ.get("TAXONKIT_DB", "/path/to/taxonkit_dump")


def reformat(taxids, fmt):
    out = subprocess.run(
        ["taxonkit", "reformat", "-I", "1", "-f", fmt, "--data-dir", TAXONKIT_DB],
        input="\n".join(taxids), capture_output=True, text=True)
    m = {}
    for l in out.stdout.strip().split("\n"):
        f = l.split("\t")
        if len(f) >= 2 and f[1].strip():
            m[f[0]] = f[1].strip()
    return m


def rapd_taxa(rank_fmt):
    strains = list({l.split("\t")[1].strip() for l in open(ACCMAP)
                    if len(l.split("\t")) >= 2 and l.split("\t")[1].strip().isdigit()})
    return set(reformat(strains, rank_fmt).values())


def tier_of(t, pg, pf, rapd_g, rapd_f):
    g, f = pg.get(t, ""), pf.get(t, "")
    if g and g in rapd_g:
        return "A"
    if f and f in rapd_f:
        return "B"
    if f:
        return "C"
    return "?"


def download(taxid):
    os.makedirs(FNA_DIR, exist_ok=True)
    out = os.path.join(FNA_DIR, "d_%s.fna" % taxid)
    if os.path.exists(out) and os.path.getsize(out) > 0:
        return out
    for mode in (["--reference"],
                 ["--assembly-source", "RefSeq", "--assembly-level", "complete,chromosome"]):
        subprocess.run(["rm", "-f", "tmp_d.zip"], cwd=FNA_DIR)
        r = subprocess.run(["datasets", "download", "genome", "taxon", taxid,
                            *mode, "--include", "genome", "--filename", "tmp_d.zip"],
                           cwd=FNA_DIR, capture_output=True)
        z = os.path.join(FNA_DIR, "tmp_d.zip")
        if r.returncode == 0 and os.path.exists(z) and os.path.getsize(z) > 0:
            subprocess.run(["unzip", "-o", "-q", "tmp_d.zip", "-d", "x_%s" % taxid],
                           cwd=FNA_DIR)
            xd = os.path.join(FNA_DIR, "x_%s" % taxid)
            fnas = []
            for root, _, files in os.walk(xd):
                fnas += [os.path.join(root, f) for f in files if f.endswith(".fna")]
            if fnas:
                os.rename(fnas[0], out)
            subprocess.run(["rm", "-rf", "x_%s" % taxid, "tmp_d.zip"], cwd=FNA_DIR)
            if os.path.exists(out):
                return out
    return None


def mash_nn(fna):
    """Minimum Mash distance from this genome to RaPDTool's DB (nearest neighbour)."""
    r = subprocess.run([MASH, "dist", MASH_DB, fna], capture_output=True, text=True)
    best = None
    for line in r.stdout.strip().split("\n"):
        f = line.split("\t")
        if len(f) >= 3:
            try:
                d = float(f[2])
            except ValueError:
                continue
            if best is None or d < best[0]:
                ref = f[0]
                import re
                m = re.search(r"GC[AF]_\d+\.\d+", ref)
                best = (d, m.group(0) if m else ref)
    return best


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--per-tier", type=int, default=8,
                    help="candidates to sample & download per tier (default 8)")
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("-o", "--output", default="data/mirror_distance.tsv")
    args = ap.parse_args()

    rows = [l.rstrip("\n").split("\t") for l in open(POOL)][1:]
    pool = [r[0] for r in rows if r and r[0].isdigit()]
    name = {r[0]: r[1] for r in rows if len(r) >= 2}
    phyl = {r[0]: r[2] for r in rows if len(r) >= 3}

    print("[dist] classifying %d pool species by tier ..." % len(pool), file=sys.stderr)
    pg, pf = reformat(pool, "{g}"), reformat(pool, "{f}")
    rapd_g, rapd_f = rapd_taxa("{g}"), rapd_taxa("{f}")
    by_tier = collections.defaultdict(list)
    for t in pool:
        by_tier[tier_of(t, pg, pf, rapd_g, rapd_f)].append(t)
    for k in "ABC":
        print("[dist] tier %s: %d species" % (k, len(by_tier[k])), file=sys.stderr)

    rng = random.Random(args.seed)
    picks = []
    for k in "ABC":
        cand = by_tier[k][:]
        rng.shuffle(cand)
        picks += [(k, t) for t in cand[:args.per_tier]]

    results = []
    for tier, t in picks:
        fna = download(t)
        if not fna:
            print("[dist] %-8s download failed" % t, file=sys.stderr)
            continue
        nn = mash_nn(fna)
        if not nn:
            continue
        d, ref = nn
        results.append((t, name.get(t, t), phyl.get(t, "?"), tier, d, ref))
        print("[dist] %-8s tier %s  dist=%.4f (%.1f%%)  %s" %
              (t, tier, d, (1 - d) * 100, name.get(t, t)[:32]), file=sys.stderr)

    results.sort(key=lambda r: r[4])
    with open(args.output, "w") as fh:
        fh.write("taxid\tspecies\tphylum\ttier\tmash_dist\tidentity_pct\tnearest_ref\n")
        for t, nm, ph, tier, d, ref in results:
            fh.write("%s\t%s\t%s\t%s\t%.4f\t%.1f\t%s\n"
                     % (t, nm, ph, tier, d, (1 - d) * 100, ref))
    print("\n[dist] wrote %s (%d genomes, distance %.3f-%.3f)"
          % (args.output, len(results),
             results[0][4] if results else 0, results[-1][4] if results else 0))
    print("[dist] now pick ~2-3 per distance band for the mock; keep the FNAs in %s"
          % FNA_DIR)


if __name__ == "__main__":
    main()
