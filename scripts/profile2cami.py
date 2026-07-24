#!/usr/bin/env python3
"""
profile2cami.py — convert taxonomic profiles to the CAMI/BIOBOXES profiling
format consumed by OPAL (https://github.com/CAMI-challenge/OPAL).

Supported inputs (--format):
  focus      FOCUS 'output_All_levels.csv' (RaPDTool profiling output)
  bracken    Bracken '*.bracken' table
  kraken     Kraken2 '--report' file
  metaphlan  MetaPhlAn 3/4 profile ('-o' output)
  truth      gold-standard mock composition: '<taxid|name><TAB|,><abundance>'
             (build the CAMI gold standard for OPAL -g from your own mock)
  auto       guess from the header (default; never resolves to 'truth')

Method (identical for every input, so the outputs are comparable):
  1. Each input is reduced to a list of (leaf, abundance) at its deepest rank
     (species for all four formats here).
  2. The leaf is resolved to an NCBI taxid — directly (bracken/kraken/metaphlan
     already carry taxids) or by name lookup (FOCUS gives names only).
  3. The *full* standard-rank lineage is derived from NCBI for that taxid, so the
     TAXPATH is authoritative and consistent across tools (it does not trust each
     tool's own rank strings — important because FOCUS uses updated phylum names
     such as 'Pseudomonadota').
  4. Abundance is propagated to every ancestor rank; each rank is normalised to
     sum to 100 %.

Taxonomy backend (--backend, default 'auto'):
  taxonkit   uses the `taxonkit` binary + a local NCBI dump (dir via --taxonkit-db
             or $TAXONKIT_DB, else taxonkit's ~/.taxonkit). No pip install, no
             download — recommended when the dump is already on disk.
  ete3       uses ete3.NCBITaxa (pip install ete3); first use downloads the NCBI
             taxonomy to ~/.etetoolkit unless --taxdb points at a pinned sqlite.
  auto       taxonkit if its binary + a usable dump are found, else ete3.
Both derive the SAME standard-rank lineage, so pick whichever is available; use a
dump recent enough that renamed phyla resolve.

The converter never fabricates abundances; unmapped entries are dropped and the
mapped fraction is reported to stderr so coverage is transparent.
"""

import argparse
import csv
import os
import shutil
import subprocess
import sys

STANDARD_RANKS = ["superkingdom", "phylum", "class", "order",
                  "family", "genus", "species"]


# --------------------------------------------------------------------------- #
# Input parsers: each returns a list of (leaf, abundance) where `leaf` is either
# an int taxid or a str name, and abundance is any non-negative number (scale is
# irrelevant — it is renormalised later).
# --------------------------------------------------------------------------- #
def parse_focus(path):
    """FOCUS output_All_levels.csv: Kingdom..Strain,<sample> ; names, no taxids."""
    leaves = []
    with open(path, newline="") as fh:
        reader = csv.reader(fh)
        header = next(reader)
        # abundance is the last column; species is column index 6
        for row in reader:
            if len(row) < 9 or not row[6]:
                continue
            species = row[6].replace("_", " ").strip()
            genus = row[5].replace("_", " ").strip()
            try:
                ab = float(row[-1])
            except ValueError:
                continue
            if ab <= 0:
                continue
            # keep genus as fallback if the species name cannot be resolved
            leaves.append((species, ab, genus))
    return leaves


def parse_bracken(path):
    """Bracken table: name, taxonomy_id, taxonomy_lvl, ..., fraction_total_reads."""
    leaves = []
    with open(path) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        idx_tid = header.index("taxonomy_id")
        # prefer the fraction column; fall back to new_est_reads
        col = "fraction_total_reads" if "fraction_total_reads" in header else "new_est_reads"
        idx_ab = header.index(col)
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) <= max(idx_tid, idx_ab):
                continue
            try:
                tid = int(f[idx_tid]); ab = float(f[idx_ab])
            except ValueError:
                continue
            if ab > 0:
                leaves.append((tid, ab, None))
    return leaves


def parse_kraken(path):
    """Kraken2 report: pct, reads_clade, reads_taxon, rank_code, taxid, name.
    Uses clade-level reads of species ('S') rows."""
    leaves = []
    with open(path) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 6:
                continue
            rank_code = f[3].strip()
            if rank_code != "S":          # species only (skip S1/S2 strain rows)
                continue
            try:
                reads_clade = float(f[1]); tid = int(f[4])
            except ValueError:
                continue
            if reads_clade > 0:
                leaves.append((tid, reads_clade, None))
    return leaves


def parse_metaphlan(path):
    """MetaPhlAn 3/4: clade_name, NCBI_tax_id, relative_abundance, ...
    Uses species rows (clade ending in s__..., no t__ SGB suffix)."""
    leaves = []
    with open(path) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 3:
                continue
            clade = f[0]
            if "|s__" not in clade and not clade.startswith("s__"):
                continue
            if "|t__" in clade:            # skip strain/SGB leaf rows
                continue
            taxpath_ids = f[1].split("|")
            try:
                tid = int(taxpath_ids[-1]); ab = float(f[2])
            except ValueError:
                continue
            if ab > 0:
                leaves.append((tid, ab, None))
    return leaves


def parse_truth(path):
    """Gold-standard mock composition: two columns, tab- or comma-separated,
    <taxid|name><delim><abundance>. Column 1 is treated as a taxid when it is a
    plain integer, otherwise as an organism name (underscores -> spaces). The
    abundance scale is arbitrary — it is renormalised per rank to 100 %. Lines
    starting with '#', blank lines, and a leading header row are skipped."""
    leaves = []
    with open(path) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 2:
                f = line.rstrip("\n").split(",")
            if len(f) < 2:
                continue
            key = f[0].strip()
            try:
                ab = float(f[1].strip())
            except ValueError:
                continue          # header row or non-numeric abundance
            if ab <= 0:
                continue
            if key.isdigit():
                leaves.append((int(key), ab, None))
            else:
                name = key.replace("_", " ").strip()
                leaves.append((name, ab, None))
    return leaves


PARSERS = {"focus": parse_focus, "bracken": parse_bracken,
           "kraken": parse_kraken, "metaphlan": parse_metaphlan,
           "truth": parse_truth}


def guess_format(path):
    with open(path) as fh:
        head = fh.readline()
    low = head.lower()
    if head.startswith("Kingdom,") or "output_all_levels" in os.path.basename(path).lower():
        return "focus"
    if "fraction_total_reads" in low or "\ttaxonomy_lvl" in low:
        return "bracken"
    if low.startswith("#mpa") or "clade_name" in low:
        return "metaphlan"
    # Kraken report has no header: 6 tab fields, 4th is a rank code letter
    f = head.rstrip("\n").split("\t")
    if len(f) >= 6 and f[3].strip()[:1].isalpha():
        return "kraken"
    raise SystemExit("ERROR: could not auto-detect format; pass --format explicitly")


# --------------------------------------------------------------------------- #
# Lineage accumulation + CAMI writer
# --------------------------------------------------------------------------- #
def resolve_leaves(leaves, ncbi):
    """Map (leaf, ab, fallback) -> leaf taxid; return (taxid, ab) and coverage."""
    resolved, mapped, total = [], 0.0, 0.0
    # Batch-resolve names for speed
    name_cache = {}
    names = [l[0] for l in leaves if isinstance(l[0], str)]
    fbs = [l[2] for l in leaves if l[2]]
    if names or fbs:
        name_cache = ncbi.get_name_translator(list(set(names + fbs)))
    for leaf, ab, fb in leaves:
        total += ab
        tid = None
        if isinstance(leaf, int):
            tid = leaf
        else:
            hit = name_cache.get(leaf)
            if hit:
                tid = hit[0]
            elif fb and name_cache.get(fb):
                tid = name_cache[fb][0]
        if tid:
            resolved.append((tid, ab)); mapped += ab
        else:
            print("  unresolved: %-45s (abundance %.4g)" % (str(leaf)[:45], ab),
                  file=sys.stderr)
    return resolved, mapped, total


def accumulate(resolved, ncbi):
    nodes = {r: {} for r in STANDARD_RANKS}
    for tid, ab in resolved:
        try:
            lineage = ncbi.get_lineage(tid)
        except Exception:
            print("  no lineage for taxid %s (abundance %.4g)" % (tid, ab),
                  file=sys.stderr)
            continue
        if not lineage:
            continue
        ranks = ncbi.get_rank(lineage)
        names = ncbi.get_taxid_translator(lineage)
        path_ids, path_sn = [], []
        for t in lineage:
            r = ranks.get(t)
            if r not in STANDARD_RANKS:
                continue
            path_ids.append(str(t))
            path_sn.append(names.get(t, str(t)))
            node = nodes[r].setdefault(
                t, {"taxpath": "|".join(path_ids),
                    "taxpathsn": "|".join(path_sn), "pct": 0.0})
            node["pct"] += ab
    return nodes


def write_cami(nodes, out, sample_id):
    with open(out, "w") as fh:
        fh.write("# Taxonomic Profiling Output\n")
        fh.write("@SampleID:%s\n" % sample_id)
        fh.write("@Version:0.10.0\n")
        fh.write("@Ranks:%s\n" % "|".join(STANDARD_RANKS))
        fh.write("@TaxonomyID:ncbi\n")
        fh.write("@@TAXID\tRANK\tTAXPATH\tTAXPATHSN\tPERCENTAGE\n")
        for rank in STANDARD_RANKS:
            entries = nodes[rank]
            total = sum(n["pct"] for n in entries.values())
            if total <= 0:
                continue
            for tid, n in sorted(entries.items(),
                                 key=lambda kv: kv[1]["pct"], reverse=True):
                pct = n["pct"] / total * 100.0
                fh.write("%s\t%s\t%s\t%s\t%.6f\n" %
                         (tid, rank, n["taxpath"], n["taxpathsn"], pct))


# --------------------------------------------------------------------------- #
# Taxonomy backends. Both expose the subset of the ete3.NCBITaxa API that
# resolve_leaves()/accumulate() use: get_name_translator, get_lineage, get_rank,
# get_taxid_translator. The taxonkit backend restricts get_lineage to the standard
# ranks (all accumulate() keeps anyway), sourced from a local NCBI dump.
# --------------------------------------------------------------------------- #
class TaxonkitBackend:
    """NCBITaxa-compatible taxonomy backed by the `taxonkit` binary + a local dump."""

    def __init__(self, data_dir=None):
        if shutil.which("taxonkit") is None:
            raise RuntimeError("taxonkit not found on PATH")
        self.data_dir = data_dir
        self._chain = {}      # leaf taxid -> [(rank, taxid, name), ...] (standard ranks)
        self._rank = {}       # taxid -> rank
        self._name = {}       # taxid -> scientific name

    def _tk(self, args, inp):
        cmd = ["taxonkit"] + args
        if self.data_dir:
            cmd += ["--data-dir", self.data_dir]
        return subprocess.run(cmd, input=inp, capture_output=True, text=True).stdout

    def get_name_translator(self, names):
        names = [n for n in names]
        if not names:
            return {}
        out = self._tk(["name2taxid"], "\n".join(names) + "\n")
        res = {}
        for line in out.splitlines():
            f = line.split("\t")
            if len(f) >= 2 and f[1].strip():
                res.setdefault(f[0], []).append(int(f[1]))
        return res

    def prime(self, taxids):
        """Batch-resolve the standard-rank chain for many taxids in two calls."""
        ids = sorted({int(t) for t in taxids if t not in self._chain})
        if not ids:
            return
        lineage_out = self._tk(["lineage"], "\n".join(map(str, ids)) + "\n")
        out = self._tk(["reformat", "-t", "-f", "{d};{p};{c};{o};{f};{g};{s}"],
                       lineage_out)
        for line in out.splitlines():
            f = line.split("\t")
            if len(f) < 4:
                continue
            leaf = int(f[0])
            chain = []
            for rank, nm, ti in zip(STANDARD_RANKS, f[2].split(";"), f[3].split(";")):
                if ti.strip() and nm.strip():
                    t2 = int(ti)
                    chain.append((rank, t2, nm))
                    self._rank[t2] = rank
                    self._name[t2] = nm
            self._chain[leaf] = chain

    def get_lineage(self, tid):
        if tid not in self._chain:
            self.prime([tid])
        return [c[1] for c in self._chain.get(tid, [])] or None

    def get_rank(self, lineage):
        return {t: self._rank.get(t, "no rank") for t in lineage}

    def get_taxid_translator(self, lineage):
        return {t: self._name.get(t, str(t)) for t in lineage}


def _taxonkit_dump_ok(data_dir):
    """True if taxonkit has a usable dump: explicit dir, $TAXONKIT_DB, or ~/.taxonkit."""
    for d in (data_dir, os.environ.get("TAXONKIT_DB"),
              os.path.expanduser("~/.taxonkit")):
        if d and os.path.exists(os.path.join(d, "names.dmp")):
            return True
    return False


def make_backend(args):
    """Build the taxonomy backend per --backend (auto prefers taxonkit if usable)."""
    want = args.backend
    if want == "auto":
        want = ("taxonkit" if shutil.which("taxonkit")
                and _taxonkit_dump_ok(args.taxonkit_db) else "ete3")
    if want == "taxonkit":
        print("[profile2cami] taxonomy backend: taxonkit", file=sys.stderr)
        return TaxonkitBackend(data_dir=args.taxonkit_db)
    try:
        from ete3 import NCBITaxa
    except ImportError:
        raise SystemExit("ERROR: ete3 backend requested but ete3 is not installed "
                         "(pip install ete3), and no taxonkit dump was found. "
                         "Either install ete3 or provide --taxonkit-db / $TAXONKIT_DB.")
    print("[profile2cami] taxonomy backend: ete3", file=sys.stderr)
    return NCBITaxa(dbfile=args.taxdb) if args.taxdb else NCBITaxa()


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("input", help="profile file to convert")
    ap.add_argument("-o", "--output", required=True, help="CAMI profile to write")
    ap.add_argument("-f", "--format", choices=list(PARSERS) + ["auto"],
                    default="auto", help="input format (default: auto)")
    ap.add_argument("-s", "--sample-id", default=None,
                    help="SampleID for the CAMI header (default: output basename)")
    ap.add_argument("--backend", choices=["auto", "taxonkit", "ete3"], default="auto",
                    help="taxonomy backend (default: auto — taxonkit if a local dump "
                         "is available, else ete3)")
    ap.add_argument("--taxonkit-db", default=os.environ.get("TAXONKIT_DB"),
                    help="taxonkit dump dir (default: $TAXONKIT_DB, else ~/.taxonkit)")
    ap.add_argument("--taxdb", default=os.environ.get("PROFILE2CAMI_TAXDB"),
                    help="path to a pinned ete3 NCBITaxa sqlite (default: "
                         "$PROFILE2CAMI_TAXDB, else ete3's ~/.etetoolkit). Build a "
                         "pinned one with setup_taxdb.py for reproducibility.")
    args = ap.parse_args()

    ncbi = make_backend(args)

    fmt = args.format if args.format != "auto" else guess_format(args.input)
    print("[profile2cami] format=%s  input=%s" % (fmt, args.input), file=sys.stderr)

    leaves = PARSERS[fmt](args.input)
    if not leaves:
        raise SystemExit("ERROR: no usable rows parsed from %s" % args.input)

    resolved, mapped, total = resolve_leaves(leaves, ncbi)
    if hasattr(ncbi, "prime"):                      # batch taxonkit lookups
        ncbi.prime([t for t, _ in resolved])
    nodes = accumulate(resolved, ncbi)
    sample = args.sample_id or os.path.splitext(os.path.basename(args.output))[0]
    write_cami(nodes, args.output, sample)

    cov = 100.0 * mapped / total if total else 0.0
    print("[profile2cami] leaves=%d  resolved=%d  abundance mapped=%.2f%%"
          % (len(leaves), len(resolved), cov), file=sys.stderr)
    print("[profile2cami] wrote %s" % args.output, file=sys.stderr)
    if cov < 90:
        print("[profile2cami] WARNING: <90%% of abundance mapped — check unresolved "
              "names above and consider a newer NCBI taxonomy dump.", file=sys.stderr)


if __name__ == "__main__":
    main()
