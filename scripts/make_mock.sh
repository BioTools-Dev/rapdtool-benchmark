#!/usr/bin/env bash
#
# make_mock.sh — build a benchmark mock community from reference genomes:
#   simulate reads (InSilicoSeq) -> assemble contigs (MegaHit) -> write the CAMI
#   gold-standard input 'mock_composition.tsv' (taxid <TAB> abundance), by mapping
#   each genome's GCA/GCF accession to its taxid and summing InSilicoSeq's
#   per-contig abundances per genome.
#
# Output (under -o, default ./mock):
#   reads_R1.fastq, reads_R2.fastq   -> for Kraken2 / MetaPhlAn (and iss truth)
#   asm/final.contigs.fa(sta)        -> for RaPDTool (full mode); .fasta copy for FOCUS
#   mock_composition.tsv             -> gold standard: feed to
#                                       'profile2cami.py -f truth'
#   make_mock.log
#
# Requires: iss (insilicoseq), megahit, python3. All in the rapdtool_bench env.
#
# Usage:
#   ./make_mock.sh -o mock GCA_x.fna GCA_y.fna GCA_z.fna ...
#   ./make_mock.sh -o mock -l genomes.list            # one genome path per line
# Options:
#   -o DIR         output dir (default ./mock)
#   -l FILE        file with one genome fasta path per line (instead of positionals)
#   -n N           total read pairs to simulate (default 2000000)
#   -m MODEL       iss error model: hiseq|miseq|novaseq (default hiseq)
#   -a DIST        iss abundance distribution: lognormal|uniform|exponential|
#                  halfnormal|zero_inflated_lognormal (default lognormal)
#   -b FILE        FIXED abundance vector (iss --abundance_file), overrides -a.
#                  Build it with make_abundance.py. Use this for a depth series:
#                  iss redraws a random composition on every run, so without a fixed
#                  vector two depths would differ in composition as well as depth.
#   --seed N       seed iss's RNGs (reproducible reads for a given -n/-b)
#   -t N           threads (default: all cores)
#   --acc-map F    GCA/GCF -> taxid table (default:
#                  $FOCUS_DB/acc_taxid_strain.tsv)
#   --no-assembly  skip MegaHit (only reads + gold standard)

set -euo pipefail

OUT="./mock"; LIST=""; NREADS=2000000; MODEL="hiseq"; ABUND="lognormal"
ABUNDFILE=""; SEED=""
THREADS="$(nproc 2>/dev/null || echo 4)"
ACCMAP="${FOCUS_DB:-/path/to/focus_build}/acc_taxid_strain.tsv"
DO_ASM=1
GENOMES=()

while [[ $# -gt 0 ]]; do
  case "$1" in
    -o) OUT="$2"; shift 2;;
    -l) LIST="$2"; shift 2;;
    -n) NREADS="$2"; shift 2;;
    -m) MODEL="$2"; shift 2;;
    -a) ABUND="$2"; shift 2;;
    -b|--abundance-file) ABUNDFILE="$2"; shift 2;;
    --seed) SEED="$2"; shift 2;;
    -t) THREADS="$2"; shift 2;;
    --acc-map) ACCMAP="$2"; shift 2;;
    --no-assembly) DO_ASM=0; shift;;
    -h|--help) sed -n '2,40p' "$0"; exit 0;;
    -*) echo "unknown option: $1" >&2; exit 2;;
    *) GENOMES+=("$1"); shift;;
  esac
done

# Collect genomes from -l if given
if [[ -n "$LIST" ]]; then
  while IFS= read -r line; do
    line="${line%%#*}"                              # strip inline comments
    line="${line#"${line%%[![:space:]]*}"}"         # ltrim whitespace
    line="${line%"${line##*[![:space:]]}"}"         # rtrim whitespace
    [[ -z "$line" ]] && continue
    eval "line=\"$line\""                           # expand ${FOCUS_DB} etc. in the path
    GENOMES+=("$line")
  done < "$LIST"
fi

if [[ ${#GENOMES[@]} -lt 2 ]]; then
  echo "ERROR: need at least 2 genome fasta files (positional or via -l)." >&2
  exit 2
fi
for g in "${GENOMES[@]}"; do
  [[ -f "$g" ]] || { echo "ERROR: genome not found: $g" >&2; exit 2; }
done
[[ -f "$ACCMAP" ]] || { echo "ERROR: acc-map not found: $ACCMAP" >&2; exit 2; }
[[ -n "$ABUNDFILE" && ! -f "$ABUNDFILE" ]] && { echo "ERROR: abundance file not found: $ABUNDFILE" >&2; exit 2; }
command -v iss  >/dev/null || { echo "ERROR: iss (insilicoseq) not on PATH." >&2; exit 2; }
[[ "$DO_ASM" == 1 ]] && { command -v megahit >/dev/null || { echo "ERROR: megahit not on PATH." >&2; exit 2; }; }

mkdir -p "$OUT"
LOG="$OUT/make_mock.log"
if [[ -n "$ABUNDFILE" ]]; then
  echo "[make_mock] $(date) genomes=${#GENOMES[@]} n_reads=$NREADS model=$MODEL abundance_file=$ABUNDFILE seed=${SEED:-none}" | tee "$LOG"
else
  echo "[make_mock] $(date) genomes=${#GENOMES[@]} n_reads=$NREADS model=$MODEL abundance=$ABUND seed=${SEED:-none}" | tee "$LOG"
fi

# --- 1. Simulate reads -------------------------------------------------------
echo "[make_mock] InSilicoSeq: simulating reads..." | tee -a "$LOG"
# -b/--abundance-file and -a are mutually exclusive in iss; a fixed vector wins.
iss_args=(--genomes "${GENOMES[@]}" --n_reads "$NREADS" --model "$MODEL"
          --cpus "$THREADS" --output "$OUT/reads")
if [[ -n "$ABUNDFILE" ]]; then
  iss_args+=(--abundance_file "$ABUNDFILE")
else
  iss_args+=(--abundance "$ABUND")
fi
[[ -n "$SEED" ]] && iss_args+=(--seed "$SEED")
iss generate "${iss_args[@]}" 2>&1 | tee -a "$LOG"
# iss writes: reads_R1.fastq, reads_R2.fastq, reads_abundance.txt

ABUND_FILE="$OUT/reads_abundance.txt"
# With --abundance_file, iss may not re-emit reads_abundance.txt (the composition was
# ours to begin with). The gold standard is built from that file, so supply it.
if [[ ! -f "$ABUND_FILE" && -n "$ABUNDFILE" ]]; then
  cp "$ABUNDFILE" "$ABUND_FILE"
  echo "[make_mock] iss did not emit reads_abundance.txt; used $ABUNDFILE" | tee -a "$LOG"
fi
[[ -f "$ABUND_FILE" ]] || { echo "ERROR: iss abundance file not found: $ABUND_FILE" >&2; exit 1; }
rm -f "$OUT"/reads.iss.tmp.*.vcf   # drop InSilicoSeq per-genome temp VCFs

# --- 2. Gold standard: per-genome abundance -> taxid -------------------------
echo "[make_mock] building gold standard mock_composition.tsv..." | tee -a "$LOG"
python3 - "$ABUND_FILE" "$ACCMAP" "$OUT/mock_composition.tsv" "${GENOMES[@]}" <<'PY' 2>&1 | tee -a "$LOG"
import re, sys
abund_file, accmap, out_tsv = sys.argv[1], sys.argv[2], sys.argv[3]
genomes = sys.argv[4:]
ACC = re.compile(r'(GC[AF]_\d+\.\d+)')

# record header -> GCA accession (from the genome file the header lives in)
hdr2acc = {}
acc_of = {}
for g in genomes:
    m = ACC.search(g.split('/')[-1])
    acc = m.group(1) if m else g.split('/')[-1]
    acc_of[g] = acc
    with open(g) as fh:
        for line in fh:
            if line.startswith('>'):
                rid = line[1:].split()[0]
                hdr2acc[rid] = acc

# sum iss per-record abundances per accession
per_acc = {}
with open(abund_file) as fh:
    for line in fh:
        f = line.rstrip('\n').split('\t')
        if len(f) < 2:
            continue
        try:
            ab = float(f[1])
        except ValueError:
            continue                       # header/non-numeric
        acc = hdr2acc.get(f[0])
        if acc is None:                    # try loose contig-id match
            acc = hdr2acc.get(f[0].split()[0])
        if acc:
            per_acc[acc] = per_acc.get(acc, 0.0) + ab

# accession -> taxid (col0=GCA, col1=taxid)
acc2taxid = {}
with open(accmap) as fh:
    for line in fh:
        f = line.rstrip('\n').split('\t')
        if len(f) >= 2 and f[1].strip().isdigit():
            acc2taxid[f[0].strip()] = f[1].strip()

# aggregate to taxid (two chosen strains may share a species taxid -> sum)
per_taxid = {}
unmapped = []
for acc, ab in per_acc.items():
    tid = acc2taxid.get(acc)
    if tid:
        per_taxid[tid] = per_taxid.get(tid, 0.0) + ab
    else:
        unmapped.append(acc)

total = sum(per_taxid.values())
with open(out_tsv, 'w') as out:
    out.write("# taxid\tabundance  (gold standard; from InSilicoSeq)\n")
    for tid, ab in sorted(per_taxid.items(), key=lambda kv: kv[1], reverse=True):
        out.write("%s\t%.6f\n" % (tid, 100.0 * ab / total if total else 0.0))

print("  genomes=%d  mapped_taxids=%d  unmapped_acc=%d"
      % (len(genomes), len(per_taxid), len(unmapped)))
if unmapped:
    print("  WARNING unmapped accessions (not in acc-map): %s" % ", ".join(unmapped[:8]))
PY

# --- 3. Assemble for RaPDTool ------------------------------------------------
if [[ "$DO_ASM" == 1 ]]; then
  echo "[make_mock] MegaHit: assembling contigs..." | tee -a "$LOG"
  rm -rf "$OUT/asm"
  megahit -1 "$OUT/reads_R1.fastq" -2 "$OUT/reads_R2.fastq" \
          -t "$THREADS" -o "$OUT/asm" 2>&1 | tee -a "$LOG"
  # RaPDTool full runs FOCUS, whose bundled v1.8 only accepts .fasta/.fna/.fastq
  # (NOT MegaHit's .fa). Provide a .fasta copy so full mode doesn't abort at FOCUS.
  cp "$OUT/asm/final.contigs.fa" "$OUT/asm/final.contigs.fasta"
fi

# --- 4. Report / next steps --------------------------------------------------
cat <<EOF | tee -a "$LOG"

[make_mock] done. Set these in run_benchmark.sh CONFIG:
  ASSEMBLY=$OUT/asm/final.contigs.fasta   (.fasta, not .fa — FOCUS needs it)
  READS_R1=$OUT/reads_R1.fastq
  READS_R2=$OUT/reads_R2.fastq
Gold standard -> CAMI:
  python3 scripts/profile2cami.py $OUT/mock_composition.tsv -f truth -o results/bench_out/gold_standard.profile -s gold
EOF
