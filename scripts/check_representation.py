#!/usr/bin/env python3
"""
check_representation.py — design the benchmark mock by measuring, for candidate
type-material genomes (from the focus_build DB), whether each one's *species* is
present in the competitor databases:

  * Kraken2 standard (full), standard16, standard8  -> parsed from each DB's
    inspect.txt (col 5 = taxid). The three share almost the same taxonomy, so
    "absent from Kraken" is judged against the full DB (robust); the full-vs-8
    difference is a run-time *sensitivity* effect, not a taxid-presence one.
  * MetaPhlAn4 (vJan26 CHOCOPhlAnSGB) -> the .pkl 'taxonomy' dict, whose values
    carry the NCBI taxid path ('2|...|<species_taxid>|'); we collect the last
    taxid of each path (+ the s__ species name as a fallback).

Each focus genome's GCA accession -> taxid (acc_taxid_strain.tsv) -> species name
and phylum (taxid_lineage.tsv). RaPDTool indexes by STRAIN taxid, so every taxid
(RaPDTool's and the competitors') is rolled up to its SPECIES taxid via taxonkit
before any membership test — see the CRITICAL note below. Requires $TAXONKIT_DB.

Verdict per genome:
  reference        present in Kraken(full+16+8) AND MetaPhlAn  -> "fair fight" set
  conflictive-both absent from Kraken(full) AND MetaPhlAn      -> strongest for RaPDTool
  conflictive-mpa  in Kraken but NOT MetaPhlAn
  conflictive-kraken NOT in Kraken but in MetaPhlAn

Usage:
  ./check_representation.py --sample 1500 [--seed 42] [-o report.tsv]
  ./check_representation.py GCA_000006945.2 GCA_000007745.1 ...
  ./check_representation.py -l candidates.txt        # one GCA per line
"""
import argparse, bz2, os, pickle, random, re, subprocess, sys

# Paths come from the environment — see config.sh (source it first). Fallbacks are
# obvious placeholders so a missing config fails loudly rather than silently.
FOCUS   = os.environ.get("FOCUS_DB", "/path/to/focus_build")
ACCMAP  = f"{FOCUS}/acc_taxid_strain.tsv"
LINEAGE = f"{FOCUS}/taxid_lineage.tsv"
GENDIR  = f"{FOCUS}/db"
TAXONKIT_DB = os.environ.get("TAXONKIT_DB", "/path/to/taxonkit_dump")

# CRITICAL: RaPDTool indexes genomes by STRAIN taxid (e.g. 469378 "C. curtum DSM 15641"),
# while Kraken2/MetaPhlAn index by SPECIES taxid (84163 "C. curtum"). Comparing the two
# directly reports a species as absent whenever it is present only under a strain taxid,
# which systematically over-counts absence. Every membership test here is therefore done
# at SPECIES level: each genome's strain taxid is rolled up to its species taxid, and the
# competitor taxid sets are rolled up too, before intersecting.
KRAKEN  = {
    "full": os.environ.get("KRAKEN_STD_INSPECT", "/path/to/standard/inspect.txt"),
    "16":   os.environ.get("KRAKEN_16_INSPECT",  "/path/to/standard16/inspect.txt"),
    "8":    os.environ.get("KRAKEN_8_INSPECT",   "/path/to/standard8/inspect.txt"),
}
MPA = os.environ.get("MPA_PKL", "/path/to/metaphlan_db/mpa_*.pkl")

# genera to prefer when picking the 'reference' set (well-known, in every DB)
PREFERRED_REF_GENERA = {
    "Escherichia", "Staphylococcus", "Pseudomonas", "Klebsiella", "Bacillus",
    "Salmonella", "Enterococcus", "Streptococcus", "Listeria", "Acinetobacter",
    "Mycobacterium", "Clostridium", "Campylobacter", "Helicobacter", "Vibrio",
}
GCA_RE = re.compile(r"(GC[AF]_\d+\.\d+)")


def norm(s):
    return re.sub(r"[^a-z0-9]", "", s.lower())


def species_taxid(taxids):
    """Map each taxid to its species-rank taxid via taxonkit (batch). A taxid already
    at or above species maps to itself where taxonkit yields no species; those are kept
    as-is so nothing is silently dropped."""
    taxids = [t for t in {str(t) for t in taxids} if t]
    if not taxids:
        return {}
    out = subprocess.run(
        ["taxonkit", "reformat", "-I", "1", "-f", "{s}", "-t",
         "--data-dir", TAXONKIT_DB],
        input="\n".join(taxids), capture_output=True, text=True)
    if out.returncode != 0:
        sys.exit("ERROR: taxonkit reformat failed: %s" % out.stderr[:300])
    m = {}
    for line in out.stdout.rstrip("\n").split("\n"):
        f = line.split("\t")
        if len(f) >= 3 and f[2].strip():
            m[f[0]] = f[2].strip().split(";")[-1]
        elif f:
            m[f[0]] = f[0]        # no species rank -> keep the taxid itself
    return m


def roll(taxid_set, mapping):
    """Species-level projection of a taxid set."""
    return {mapping.get(t, t) for t in taxid_set}


def load_lineage():
    lin = {}
    with open(LINEAGE) as fh:
        for l in fh:
            f = l.rstrip("\n").split("\t")
            if len(f) >= 8:
                lin[f[0]] = {"phylum": f[2], "genus": f[6], "species": f[7]}
    return lin


def load_accmap():
    a = {}
    with open(ACCMAP) as fh:
        for l in fh:
            f = l.rstrip("\n").split("\t")
            if len(f) >= 2 and f[1].strip().isdigit():
                a[f[0].strip()] = f[1].strip()
    return a


def load_kraken(path):
    s = set()
    with open(path) as fh:
        for l in fh:
            if l.startswith("#"):
                continue
            f = l.rstrip("\n").split("\t")
            if len(f) >= 5:
                s.add(f[4].strip())
    return s


def load_mpa(path):
    try:
        db = pickle.load(open(path, "rb"))
    except Exception:
        db = pickle.load(bz2.open(path, "rb"))
    taxids, names = set(), set()
    for key, val in db.get("taxonomy", {}).items():
        try:
            tp = [x for x in str(val[0]).split("|") if x]
            if tp:
                taxids.add(tp[-1])
        except Exception:
            pass
        m = re.search(r"s__([^|]+)", key)
        if m:
            names.add(norm(m.group(1)))
    # merged taxid remaps (old->new); index both directions for lookups
    merged = db.get("merged_taxon", {}) or {}
    return taxids, names, merged


def classify(gca, accmap, lineage, krak, mpa_ids, mpa_names, merged, sp_of):
    tid = accmap.get(gca)
    if not tid:
        return None
    info = lineage.get(tid)
    if not info:
        return None
    # Roll this genome's strain taxid up to its species taxid; all competitor sets in
    # `krak`/`mpa_ids` are already species-level (see main()), so we compare like with like.
    sp = sp_of.get(tid, tid)
    mpa_tids = {sp}
    if tid in merged:
        v = merged[tid]
        mpa_tids.add(str(v[0]) if isinstance(v, (list, tuple)) else str(v))
    in_mpa = bool(mpa_tids & mpa_ids) or (norm(info["species"]) in mpa_names)
    return {
        "gca": gca, "taxid": tid, "species_taxid": sp, "phylum": info["phylum"],
        "genus": info["genus"], "species": info["species"],
        "k_full": sp in krak["full"], "k16": sp in krak["16"], "k8": sp in krak["8"],
        "mpa": in_mpa,
    }


def verdict(r):
    if r["k_full"] and r["k16"] and r["k8"] and r["mpa"]:
        return "reference"
    if not r["k_full"] and not r["mpa"]:
        return "conflictive-both"
    if r["k_full"] and not r["mpa"]:
        return "conflictive-mpa"
    if not r["k_full"] and r["mpa"]:
        return "conflictive-kraken"
    return "partial"


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("gca", nargs="*", help="GCA/GCF accessions or genome paths")
    ap.add_argument("-l", "--list", help="file with one GCA/path per line")
    ap.add_argument("--sample", type=int, help="randomly sample N genomes from focus_build/db")
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("-o", "--output", help="write full TSV report here")
    ap.add_argument("--n-ref", type=int, default=10, help="reference genomes to suggest")
    ap.add_argument("--n-conf", type=int, default=10, help="conflictive genomes to suggest")
    args = ap.parse_args()

    print("[load] lineage, acc-map, Kraken inspect.txt (x3), MetaPhlAn pkl ...",
          file=sys.stderr)
    lineage = load_lineage()
    accmap = load_accmap()
    krak = {k: load_kraken(v) for k, v in KRAKEN.items()}
    mpa_ids, mpa_names, merged = load_mpa(MPA)

    # Roll every taxid set up to species so all comparisons are species-vs-species.
    print("[load] rolling taxids up to species (taxonkit) ...", file=sys.stderr)
    all_tids = set(accmap.values()) | mpa_ids
    for s in krak.values():
        all_tids |= s
    sp_of = species_taxid(all_tids)
    krak = {k: roll(v, sp_of) for k, v in krak.items()}
    mpa_ids = roll(mpa_ids, sp_of)

    print("[load] kraken full/16/8 species=%d/%d/%d  mpa species=%d names=%d"
          % (len(krak["full"]), len(krak["16"]), len(krak["8"]),
             len(mpa_ids), len(mpa_names)), file=sys.stderr)

    # collect candidate accessions
    cands = []
    for x in args.gca:
        m = GCA_RE.search(os.path.basename(x))
        cands.append(m.group(1) if m else x)
    if args.list:
        with open(args.list) as fh:
            for l in fh:
                l = l.strip()
                if not l or l.startswith("#"):
                    continue
                m = GCA_RE.search(os.path.basename(l))
                cands.append(m.group(1) if m else l)
    if args.sample:
        files = [f for f in os.listdir(GENDIR) if f.endswith(".fna")]
        random.seed(args.seed)
        random.shuffle(files)
        for f in files:
            m = GCA_RE.search(f)
            if m:
                cands.append(m.group(1))
            if len(cands) >= args.sample:
                break
    if not cands:
        ap.error("provide accessions, -l, or --sample")

    rows = []
    for gca in cands:
        r = classify(gca, accmap, lineage, krak, mpa_ids, mpa_names, merged, sp_of)
        if r:
            r["verdict"] = verdict(r)
            rows.append(r)

    # tallies
    from collections import Counter
    tally = Counter(r["verdict"] for r in rows)
    print("\n[classified %d/%d genomes]" % (len(rows), len(cands)), file=sys.stderr)
    for k in ("reference", "conflictive-both", "conflictive-mpa",
              "conflictive-kraken", "partial"):
        print("  %-18s %d" % (k, tally.get(k, 0)), file=sys.stderr)

    # optional full TSV
    hdr = ["gca", "taxid", "species_taxid", "phylum", "genus", "species",
           "kraken_full", "kraken16", "kraken8", "metaphlan", "verdict"]
    if args.output:
        with open(args.output, "w") as out:
            out.write("\t".join(hdr) + "\n")
            for r in sorted(rows, key=lambda r: (r["verdict"], r["phylum"], r["species"])):
                out.write("\t".join([
                    r["gca"], r["taxid"], r["species_taxid"], r["phylum"], r["genus"],
                    r["species"], "Y" if r["k_full"] else "N", "Y" if r["k16"] else "N",
                    "Y" if r["k8"] else "N", "Y" if r["mpa"] else "N", r["verdict"],
                ]) + "\n")
        print("[wrote] %s" % args.output, file=sys.stderr)

    # ---- suggest a balanced pick -------------------------------------------
    def pick_diverse(pool, n, prefer_genera=None):
        chosen, seen_phyla, seen_sp = [], set(), set()
        pool = list(pool)
        if prefer_genera:
            pool.sort(key=lambda r: (r["genus"] not in prefer_genera, r["phylum"]))
        # first pass: one per phylum, unique species
        for r in pool:
            if len(chosen) >= n:
                break
            if r["species"] in seen_sp:
                continue
            if r["phylum"] not in seen_phyla:
                chosen.append(r); seen_phyla.add(r["phylum"]); seen_sp.add(r["species"])
        # second pass: fill remaining, still unique species
        for r in pool:
            if len(chosen) >= n:
                break
            if r["species"] in seen_sp:
                continue
            chosen.append(r); seen_sp.add(r["species"])
        return chosen

    ref_pool = [r for r in rows if r["verdict"] == "reference"]
    conf_both = [r for r in rows if r["verdict"] == "conflictive-both"]
    conf_other = [r for r in rows if r["verdict"] in ("conflictive-mpa", "conflictive-kraken")]

    ref = pick_diverse(ref_pool, args.n_ref, PREFERRED_REF_GENERA)
    conf = pick_diverse(conf_both, args.n_conf)
    if len(conf) < args.n_conf:  # top up with single-tool-conflictive
        conf += pick_diverse([r for r in conf_other if r not in conf],
                             args.n_conf - len(conf))

    def show(title, items):
        print("\n" + title)
        print("  %-18s %-12s %-11s %-9s %s" %
              ("species", "taxid", "kraken(f/16/8)", "metaphlan", "phylum"))
        for r in items:
            print("  %-18.18s %-12s %d/%d/%d        %-9s %s" % (
                r["species"], r["taxid"], r["k_full"], r["k16"], r["k8"],
                "Y" if r["mpa"] else "N", r["phylum"]))
        print("  GCA list:", " ".join(GENDIR + "/" + r["gca"] + ".fna" for r in items))

    show("=== REFERENCE set (present in all DBs) ===", ref)
    show("=== CONFLICTIVE set (absent from competitors; RaPDTool-only) ===", conf)

    print("\nAll %d GCA for make_mock.sh:" % (len(ref) + len(conf)))
    print(" ".join(GENDIR + "/" + r["gca"] + ".fna" for r in ref + conf))


if __name__ == "__main__":
    main()
