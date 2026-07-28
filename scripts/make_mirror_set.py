#!/usr/bin/env python3
"""
make_mirror_set.py — select the MIRROR genome set: species Kraken2 has and RaPDTool
does not.

Why this exists
---------------
The conflictive mock asks "what happens to taxa only RaPDTool's database contains?".
This asks the mirror question, which is a SAFETY property rather than a marketing one:

    when RaPDTool is given an organism absent from its database, does it stay silent,
    or does it confidently assign it to the nearest type strain?

Declaring a domain limitation in the manuscript ("not for uncultured organisms") is a
scope statement, not evidence. Users will feed environmental data to the tool whatever
the paper says, so the failure MODE has to be characterised:

    (a) reports nothing / unclassified        -> graceful failure, a stated limitation
    (b) reports a wrong species with a
        plausible-looking Mash distance       -> silent misassignment, which would have
                                                 to be disclosed and guarded against

(a) is what the distance series measured: across fourteen genomes the resolved rank
tracks Mash distance, and no genome below ~80 % identity received any assignment
(benchmark_rationale.md section 4b).

Selection
---------
Population: species present in the Kraken2 standard database but absent from
RaPDTool's type-material set. Sampled stratified by phylum (one per phylum) with a
fixed seed, exactly as the conflictive set was, so the two experiments are symmetric
and equally defensible.

Usage
-----
  ./make_mirror_set.py --sample 10 --seed 42 -o mirror_genomes.tsv
  # then download the genomes it lists (it prints the datasets command)
"""
import argparse
import collections
import os
import random
import subprocess
import sys

FOCUS_LINEAGE = os.environ.get("FOCUS_DB", "/path/to/focus_build") + "/taxid_lineage.tsv"
KRAKEN_INSPECT = os.environ.get("KRAKEN_STD_INSPECT", "/path/to/standard/inspect.txt")
TAXONKIT_DB = os.environ.get("TAXONKIT_DB", "/path/to/taxonkit_dump")

# Patterns marking an organism as uncultured / not formally described. Not a selection
# criterion -- reported so the manuscript can state what fraction of the mirror set is
# genuinely undescribed rather than merely missing from the type-material set.
#
# These MUST be word-boundary anchored. A plain substring test for "bacterium " matches
# inside ordinary genus names (Crypto|bacterium curtum, Exiguo|bacterium, Acido|bacterium)
# and misclassifies validly described species as undescribed: on this pool it reported
# 85.2 % undescribed where the anchored patterns report 84.0 %.
import re as _re
UNCULTURED_RE = _re.compile(
    r"(^candidatus\b"          # Candidatus Xxx
    r"|\buncultured\b"
    r"|\bmetagenome\b"
    r"|^(bacterium|archaeon)\b"      # bare 'bacterium XYZ123'
    r"|\s(bacterium|archaeon)\b"     # '... bacterium XYZ' as its own word
    r"|\ssp\."                        # Genus sp. STRAIN
    r"|\b(endo)?symbiont\b)", _re.I)


def looks_undescribed(name):
    return bool(UNCULTURED_RE.search(name or ""))


ACCMAP = os.environ.get("FOCUS_DB", "/path/to/focus_build") + "/acc_taxid_strain.tsv"


def species_of(taxids):
    """Roll each taxid up to its species-rank taxid via taxonkit (batch)."""
    taxids = [t for t in {str(t) for t in taxids} if t]
    out = subprocess.run(
        ["taxonkit", "reformat", "-I", "1", "-f", "{s}", "-t", "--data-dir", TAXONKIT_DB],
        input="\n".join(taxids), capture_output=True, text=True)
    if out.returncode != 0:
        sys.exit("ERROR: taxonkit reformat failed: %s" % out.stderr[:300])
    m = {}
    for line in out.stdout.rstrip("\n").split("\n"):
        f = line.split("\t")
        m[f[0]] = f[2].strip().split(";")[-1] if len(f) >= 3 and f[2].strip() else f[0]
    return m


def focus_species():
    """RaPDTool's genomes rolled up to SPECIES taxids.

    RaPDTool indexes by strain taxid; comparing those directly against Kraken's species
    taxids reports a species as RaPDTool-only whenever it is present under a strain taxid.
    Rolling up to species is required for the comparison to mean anything — without it,
    4 of 10 mirror picks were species RaPDTool actually contains (verified: same strain,
    mash identity 1.0)."""
    strains = set()
    for line in open(ACCMAP):
        f = line.split("\t")
        if len(f) >= 2 and f[1].strip().isdigit():
            strains.add(f[1].strip())
    sp = species_of(strains)
    return {sp.get(t, t) for t in strains}


def kraken_species():
    """Species-rank taxids in the Kraken2 database (col4 = rank code, col5 = taxid)."""
    s = set()
    for line in open(KRAKEN_INSPECT):
        f = line.rstrip("\n").split("\t")
        if len(f) >= 5 and f[3].strip() == "S":
            s.add(f[4].strip())
    return s


def taxonkit_lineage(taxids):
    """taxid -> (name, phylum, superkingdom). One subprocess call for the whole list."""
    p = subprocess.run(
        ["taxonkit", "lineage", "-R", "--data-dir", TAXONKIT_DB],
        input="\n".join(taxids), capture_output=True, text=True)
    if p.returncode != 0:
        sys.exit("ERROR: taxonkit failed: %s" % p.stderr.strip()[:300])
    out = {}
    for line in p.stdout.rstrip("\n").split("\n"):
        f = line.split("\t")
        if len(f) < 3:
            continue
        tid, lineage, ranks = f[0], f[1].split(";"), f[2].split(";")
        pick = lambda want: next((n for n, r in zip(lineage, ranks) if r == want), "")
        sk = pick("superkingdom") or pick("domain") or "unclassified"
        out[tid] = (lineage[-1] if lineage else tid,
                    pick("phylum") or "unclassified", sk)
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--sample", type=int, default=10, help="genomes to select (default 10)")
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("-o", "--output", default="data/mirror_genomes.tsv")
    ap.add_argument("--pool-out", default="", help="also write the full candidate pool here")
    ap.add_argument("--include-nonprokaryotic", action="store_true",
                    help="do NOT restrict to Bacteria/Archaea. Off by default: Kraken2 "
                         "standard also indexes viruses and human, which RaPDTool's "
                         "prokaryotic type-material database never claimed to cover, so "
                         "including them would test nothing about its failure mode.")
    args = ap.parse_args()

    focus, kraken = focus_species(), kraken_species()
    cand = sorted(kraken - focus, key=int)

    print("[mirror] RaPDTool species (strain taxids rolled up): %d" % len(focus))
    print("[mirror] Kraken2 standard species      : %d" % len(kraken))
    print("[mirror] present in BOTH               : %d" % len(kraken & focus))
    print("[mirror] Kraken2 only (mirror pool)    : %d  (%.1f %% of Kraken2's species)"
          % (len(cand), 100.0 * len(cand) / len(kraken)))

    info = taxonkit_lineage(cand)

    # Break the pool down by superkingdom before selecting: much of Kraken2's species
    # surplus is viral, which is outside what RaPDTool's database ever claimed.
    sk_counts = collections.Counter(info.get(t, ("", "", "unclassified"))[2] for t in cand)
    print("[mirror] mirror pool by superkingdom:")
    for sk, n in sk_counts.most_common():
        print("           %-16s %6d  (%.1f %%)" % (sk or "unclassified", n,
                                                   100.0 * n / len(cand)))

    if not args.include_nonprokaryotic:
        cand = [t for t in cand
                if info.get(t, ("", "", ""))[2] in ("Bacteria", "Archaea")]
        print("[mirror] restricted to Bacteria/Archaea : %d species" % len(cand))
        if not cand:
            sys.exit("ERROR: no prokaryotic candidates left.")

    by_phylum = collections.defaultdict(list)
    for t in cand:
        name, phylum, _ = info.get(t, (t, "unclassified", ""))
        by_phylum[phylum].append((t, name))

    unc = sum(1 for t in cand if looks_undescribed(info.get(t, ("", "", ""))[0]))
    print("[mirror] phyla represented             : %d" % len(by_phylum))
    print("[mirror] names marked uncultured/undescribed: %d (%.1f %%)"
          % (unc, 100.0 * unc / max(len(cand), 1)))

    # Stratified: one per phylum, largest phyla first, then fill.
    rng = random.Random(args.seed)
    order = sorted(by_phylum, key=lambda p: -len(by_phylum[p]))
    picked = []
    for ph in order:
        if len(picked) >= args.sample:
            break
        if ph == "unclassified":
            continue
        picked.append((ph,) + rng.choice(by_phylum[ph]))
    i = 0
    while len(picked) < args.sample and i < len(order):
        ph = order[i]; i += 1
        pool = [x for x in by_phylum[ph] if x[0] not in {p[1] for p in picked}]
        if pool:
            picked.append((ph,) + rng.choice(pool))

    with open(args.output, "w") as fh:
        fh.write("# MIRROR set: species in Kraken2 standard, ABSENT from RaPDTool's "
                 "type-material database.\n")
        fh.write("# Population: %d species (%.1f %% of Kraken2's species); "
                 "stratified by phylum, seed %d.\n"
                 % (len(cand), 100.0 * len(cand) / len(kraken), args.seed))
        fh.write("# Purpose: characterise RaPDTool's OUT-OF-DOMAIN failure mode "
                 "(silence vs. confident misassignment).\n")
        fh.write("taxid\tspecies\tphylum\tlooks_uncultured\n")
        for ph, t, name in picked:
            flag = "yes" if looks_undescribed(name) else "no"
            fh.write("%s\t%s\t%s\t%s\n" % (t, name, ph, flag))

    if args.pool_out:
        with open(args.pool_out, "w") as fh:
            fh.write("taxid\tspecies\tphylum\tsuperkingdom\n")
            for t in cand:
                n, ph, sk = info.get(t, (t, "unclassified", ""))
                fh.write("%s\t%s\t%s\t%s\n" % (t, n, ph, sk))
        print("[mirror] full pool -> %s" % args.pool_out)

    print("\n%-10s %-46s %-24s %s" % ("taxid", "species", "phylum", "uncultured?"))
    for ph, t, name in picked:
        flag = "yes" if looks_undescribed(name) else "no"
        print("%-10s %-46s %-24s %s" % (t, name[:45], ph[:23], flag))

    print("\n[mirror] wrote %s" % args.output)
    print("[mirror] download the genomes with:")
    print("  cut -f1 %s | tail -n +2 | while read t; do" % args.output)
    print("    datasets download genome taxon \"$t\" --reference \\")
    print("      --include genome --filename \"mirror_$t.zip\"; done")


if __name__ == "__main__":
    main()
