#!/usr/bin/env python3
"""
make_abundance.py — build a FIXED per-contig abundance vector for InSilicoSeq.

Why this exists: 'iss generate -a lognormal' draws a NEW random composition on every
run, so two runs at different depths would differ in composition as well as in depth,
confounding the two variables. Writing the vector once, to a versioned file, makes the
depth series a clean one-variable experiment and makes the mock reproducible without
depending on iss's RNG.

Design of the vector:
  * Log-spaced genome abundances spanning --range (default 20x). This keeps the
    community realistically uneven -- real metagenomes are never uniform -- while
    avoiding the pathological draw we got from plain lognormal (one taxon at 61 %,
    1500x dynamic range), where the rare members stayed below 0.25x coverage even at
    30 M reads and could not assemble at any affordable depth.
  * Abundance ranks are INTERLEAVED between the 'reference' and 'conflictive' halves
    of the genome list, so both groups span the same abundance range. Without this,
    abundance would correlate with database representation and the headline result
    (conflictive taxa detected only by RaPDTool) could be dismissed as an abundance
    artefact rather than a database-coverage effect.
  * Within a genome, abundance is split across contigs proportionally to contig
    length, so coverage is even along the genome (iss computes coverage per record).

Output: '<record_id>\\tabundance' summing to 1.0 -- the format 'iss --abundance_file'
expects, identical to the 'reads_abundance.txt' iss writes itself.

Usage:
  ./make_abundance.py -l mock_genomes.list -o mock_abundance.txt --depths 3 10 30
"""
import argparse
import os
import sys

RL = 126  # iss hiseq model read length, for the coverage preview


def parse_list(path):
    """Return [(fasta_path, label, group)] preserving file order.

    Group is taken from the '--- CONFLICTIVE' banner in mock_genomes.list, so the
    split stays in sync with the list instead of being duplicated here.
    """
    entries, group = [], "reference"
    for raw in open(path):
        if raw.lstrip().startswith("#"):
            # Only the '# --- CONFLICTIVE ...' section banner flips the group; a bare
            # mention of the word (the file header says "10 reference + 10 conflictive")
            # must not.
            if raw.lstrip().startswith("# ---") and "CONFLICTIVE" in raw.upper():
                group = "conflictive"
            continue
        code, _, comment = raw.partition("#")
        code = code.strip()
        if not code:
            continue
        label = comment.split("(")[0].strip() or code.split("/")[-1]
        # The shipped lists write genome paths as ${FOCUS_DB}/db/GCA_x.fna so they are
        # portable; expand them here the way make_mock.sh does, and fail loudly if the
        # variable is unset rather than trying to open a literal '${FOCUS_DB}/...'.
        path_ = os.path.expandvars(os.path.expanduser(code))
        if "$" in path_:
            sys.exit("ERROR: unexpanded variable in %s\n       %s\n"
                     "       Did you 'source config.sh'? (see README, Configuration)"
                     % (path, path_))
        entries.append((path_, label, group))
    return entries


def contigs(path):
    """Return [(record_id, length)] for a fasta."""
    out, rid, n = [], None, 0
    for line in open(path):
        if line.startswith(">"):
            if rid is not None:
                out.append((rid, n))
            rid, n = line[1:].split()[0], 0
        else:
            n += len(line.strip())
    if rid is not None:
        out.append((rid, n))
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("-l", "--list", required=True, help="mock_genomes.list")
    ap.add_argument("-o", "--output", required=True, help="abundance file for iss -b")
    ap.add_argument("--range", type=float, default=20.0,
                    help="dynamic range, most:least abundant genome (default 20)")
    ap.add_argument("--equal-coverage", action="store_true",
                    help="give every genome the SAME COVERAGE (abundance proportional "
                         "to genome length), for the positive control. Note that equal "
                         "abundance is NOT equal coverage: a large genome given the "
                         "same read share is sequenced more shallowly. Overrides --range.")
    ap.add_argument("--depths", type=float, nargs="*", default=[3, 10, 30],
                    help="read totals in millions, for the coverage preview")
    args = ap.parse_args()

    entries = parse_list(args.list)
    if not entries:
        sys.exit("ERROR: no genomes parsed from %s" % args.list)

    sizes = {p: contigs(p) for p, _, _ in entries}
    total_bp = {p: sum(n for _, n in c) for p, c in sizes.items()}

    n = len(entries)
    if args.equal_coverage:
        # Coverage_i = N * frac_i * RL / size_i, so frac_i proportional to size_i
        # makes coverage identical for every genome.
        gb = sum(total_bp.values())
        assign = {p: total_bp[p] / gb for p, _, _ in entries}
        write_vector(entries, sizes, total_bp, assign, args)
        report(entries, total_bp, assign, args, equal=True)
        return

    # Log-spaced weights: rank 0 = most abundant, rank n-1 = least.
    weights = [args.range ** (-i / (n - 1)) for i in range(n)]
    s = sum(weights)
    fracs = [w / s for w in weights]

    # Interleave the two groups across the abundance ranks so neither group is
    # systematically more abundant than the other.
    ref = [e for e in entries if e[2] == "reference"]
    conf = [e for e in entries if e[2] == "conflictive"]
    # Alternate which group takes the higher rank of each consecutive pair, so neither
    # group is systematically one rank more abundant than the other.
    ordered = []
    for i in range(max(len(ref), len(conf))):
        first, second = (conf, ref) if i % 2 == 0 else (ref, conf)
        for src in (first, second):
            if i < len(src):
                ordered.append(src[i])

    assign = {p: fracs[i] for i, (p, _, _) in enumerate(ordered)}
    write_vector(entries, sizes, total_bp, assign, args)

    print("[make_abundance] dynamic range %.0fx  (max %.2f %%, min %.2f %%)"
          % (args.range, max(fracs) * 100, min(fracs) * 100))
    print("[make_abundance] reference/conflictive ranks interleaved "
          "(mean abundance %.2f %% vs %.2f %%)"
          % (100 * sum(assign[p] for p, _, g in entries if g == "reference") / max(len(ref), 1),
             100 * sum(assign[p] for p, _, g in entries if g == "conflictive") / max(len(conf), 1)))

    report(entries, total_bp, assign, args)


def write_vector(entries, sizes, total_bp, assign, args):
    """Split each genome's fraction across its contigs, proportional to length.

    This is the whole point of generating the vector ourselves: iss's own --abundance
    distributions operate PER RECORD, so for draft multi-contig genomes the realised
    per-genome abundance ends up tracking contig count rather than the requested
    distribution. (Observed: 'iss --abundance uniform' on this genome set gave one
    1,491-contig genome 61 % of the reads and a single-contig genome 0.04 %.)
    """
    lines, check = [], 0.0
    for path, _, _ in entries:
        gfrac, gbp = assign[path], total_bp[path]
        for rid, ln in sizes[path]:
            f = gfrac * ln / gbp
            lines.append("%s\t%.12f" % (rid, f))
            check += f
    # Absorb float drift into the last record so iss sees a vector summing to 1.
    if lines and abs(check - 1.0) > 1e-12:
        rid, val = lines[-1].split("\t")
        lines[-1] = "%s\t%.12f" % (rid, float(val) + (1.0 - check))

    with open(args.output, "w") as fh:
        fh.write("\n".join(lines) + "\n")
    print("[make_abundance] %d genomes, %d contigs -> %s" %
          (len(entries), len(lines), args.output))


def report(entries, total_bp, assign, args, equal=False):
    if equal:
        gb = sum(total_bp.values())
        print("[make_abundance] EQUAL-COVERAGE mode: abundance proportional to genome "
              "length,\n                 so every genome receives the same coverage "
              "(positive control).")
        for d in args.depths:
            print("                 %gM reads -> %.1fx for all %d genomes"
                  % (d, d * 1e6 * RL / gb, len(entries)))
    if not args.depths:
        return
    hdr = "".join("%12s" % ("%gM" % d) for d in args.depths)
    print("\n%-34s%8s%9s%s" % ("genome", "Mb", "abund%", hdr))
    for path, label, group in sorted(entries, key=lambda e: -assign[e[0]]):
        cov = "".join("%11.1fx" % (d * 1e6 * assign[path] * RL / total_bp[path])
                      for d in args.depths)
        tag = "*" if group == "conflictive" else " "
        print("%s%-33s%8.2f%9.2f%s" %
              (tag, label[:33], total_bp[path] / 1e6, assign[path] * 100, cov))
    print("\n* = conflictive (absent from Kraken2 and MetaPhlAn databases)")
    gb = sum(total_bp.values())
    for d in args.depths:
        print("  %gM reads -> %.0f Mbp, mean coverage %.1fx, ~%.1f GB fastq"
              % (d, d * 1e6 * RL / 1e6, d * 1e6 * RL / gb, d * 1e6 * 268 / 1e9))


if __name__ == "__main__":
    main()
