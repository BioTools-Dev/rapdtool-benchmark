#!/usr/bin/env python3
"""
mash_detection.py — RaPDTool's species DETECTION from the mash confidence table.

RaPDTool has two species-level outputs and they must not be conflated:
  * the mash-screen confidence table (rapdtool_confidence.tbl) — its confident species
    DETECTION (high-containment reference genomes), analogous to MetaPhlAn's
    marker-filtered species list;
  * the FOCUS profile (output_All_levels.csv) — its relative-ABUNDANCE estimate, which
    RaPDTool itself labels "be cautious at species taxonomic level".

Detection metrics (recall / precision / F1, and the reference/conflictive split) must
come from the mash table; abundance metrics (Bray-Curtis, L1) from FOCUS. Using FOCUS
for detection understates RaPDTool badly (on the mocks: FOCUS ~204 species, ~184 false
positives; mash: exactly the true species, 0 false positives).

This tool parses the mash table (screen mode), rolls taxids to species via taxonkit,
and scores detection against a CAMI gold profile (and, for the mocks, splits the result
into the reference and conflictive halves from mock_genomes.list).

Usage:
  ./mash_detection.py --bench bench_ln_30M --gold bench_ln_30M/gold_standard.profile \\
      --split mock_genomes.list
  ./mash_detection.py --bench bench_zymo_even --gold bench_zymo_even/gold_dna.profile \\
      --extra-true 96241:1423   # credit B. spizizenii as the Zymo B. subtilis member
"""
import argparse
import os
import re
import subprocess
import sys

TAXONKIT_DB = os.environ.get("TAXONKIT_DB", "/path/to/taxonkit_dump")


def mash_taxids(conf_tbl):
    """taxids in the 'Reference genomes detected (mash screen)' block."""
    if not os.path.exists(conf_tbl):
        return set()
    txt = open(conf_tbl).read()
    if "Reference genomes detected" not in txt:
        return set()          # full mode has no mash-screen table
    block = txt.split("Reference genomes detected")[1].split("FOCUS")[0]
    ids = set()
    for line in block.split("\n"):
        m = re.match(r"\|\s*[A-Za-z].*?\|\s*(\d+)\s*\|", line)
        if m:
            ids.add(m.group(1))
    return ids


def gold_species(profile):
    s = set()
    for line in open(profile):
        if line.startswith(("@", "#")) or not line.strip():
            continue
        f = line.rstrip("\n").split("\t")
        if len(f) >= 2 and f[1] == "species":
            s.add(f[0])
    return s


def rollup(taxids):
    taxids = [t for t in {str(x) for x in taxids} if t]
    if not taxids:
        return {}
    out = subprocess.run(
        ["taxonkit", "reformat", "-I", "1", "-f", "{s}", "-t", "--data-dir", TAXONKIT_DB],
        input="\n".join(taxids), capture_output=True, text=True)
    m = {}
    for line in out.stdout.rstrip("\n").split("\n"):
        f = line.split("\t")
        m[f[0]] = f[2].split(";")[-1] if len(f) >= 3 and f[2].strip() else f[0]
    return m


def split_halves(list_path):
    ref, conf, is_conf = [], [], False
    for line in open(list_path):
        if line.lstrip().startswith("#"):
            if "---" in line and "CONFLICTIVE" in line.upper():
                is_conf = True
            continue
        m = re.search(r"\((\d+)\)", line)
        if m:
            (conf if is_conf else ref).append(m.group(1))
    return ref, conf


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--bench", required=True, help="bench_<dataset> directory")
    ap.add_argument("--gold", required=True, help="CAMI gold profile")
    ap.add_argument("--split", help="mock_genomes.list, for the reference/conflictive split")
    ap.add_argument("--extra-true", default="",
                    help="comma-sep pairs detected:goldmember to credit a reclassified "
                         "taxon, e.g. 96241:1423 (B. spizizenii counts as the B. subtilis "
                         "member)")
    args = ap.parse_args()

    conf_tbl = os.path.join(args.bench, "rapdtool_screen.rep1", "rapdtool_confidence.tbl")
    detected = mash_taxids(conf_tbl)
    gold = gold_species(args.gold)

    equiv = {}
    for pair in filter(None, args.extra_true.split(",")):
        a, b = pair.split(":")
        equiv[a] = b

    ref, conf = ([], [])
    if args.split:
        ref, conf = split_halves(args.split)

    # Roll EVERY taxid we will compare (detected, gold, split halves, equiv) so vintage
    # differences (e.g. Gloeothece 497965 vs 2546359) collapse to one species taxid.
    everything = detected | gold | set(equiv) | set(equiv.values()) | set(ref) | set(conf)
    sp = rollup(everything)

    # Credit a reclassified detection as its gold member: REPLACE the detected taxid
    # with the gold-side equivalent (e.g. B. spizizenii 96241 -> B. subtilis 1423), so
    # it is counted once as a true positive, not once as TP and once as FP.
    det_raw = set(detected)
    for a, b in equiv.items():
        if a in det_raw:
            det_raw.discard(a)
            det_raw.add(b)
    det = {sp.get(t, t) for t in det_raw}
    golds = {sp.get(t, t) for t in gold}

    tp = len(det & golds)
    fp = len(det - golds)
    fn = len(golds - det)
    rec = tp / len(golds) if golds else 0
    prec = tp / (tp + fp) if det else 0
    f1 = 2 * prec * rec / (prec + rec) if (prec + rec) else 0

    print("[mash-detection] %s" % args.bench)
    print("  detected=%d  gold=%d  TP=%d  FP=%d  FN=%d" % (len(det), len(golds), tp, fp, fn))
    print("  recall=%.3f  precision=%.3f  F1=%.3f" % (rec, prec, f1))

    if args.split:
        rs = {sp.get(t, t) for t in ref}
        cs = {sp.get(t, t) for t in conf}
        print("  reference detected: %d/%d" % (len(det & rs), len(rs)))
        print("  conflictive detected: %d/%d" % (len(det & cs), len(cs)))


if __name__ == "__main__":
    main()
