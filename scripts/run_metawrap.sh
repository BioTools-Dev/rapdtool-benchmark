#!/usr/bin/env bash
# run_metawrap.sh <dataset> <mock_dir>
#
# Fair MAG-recovery comparison. Runs MetaWRAP binning + bin_refinement on the SAME
# assembly RaPDTool consumed (only the binning/refinement stage differs), then RE-SCORES
# both bin sets with the SAME evaluator (miComplete/Bact105, inside RaPDTool's SIF) so the
# comparison does not reward MetaWRAP's own CheckM objective. Containerised end to end —
# no dependency install. Records MAGs, wall time and peak RSS to results/metawrap/.
#
# Portable via env (see config.sh); defaults suit the study machine:
#   METAWRAP_SIF   apptainer pull docker://quay.io/biocontainers/metawrap-mg:1.3.0--hdfd78af_1
#   CHECKM_DB      https://data.ace.uq.edu.au/public/CheckM_databases/checkm_data_2015_01_16.tar.gz
#   RAPDTOOL_SIF   RaPDTool Apptainer image (provides miComplete for the common re-scoring)
#   THREADS        (default 16)   METAWRAP_OUT (default /tmp/metawrap_<dataset>)
#
# Outputs (results/metawrap/): <ds>.checkm.stats, <ds>.micomplete.tab, summary.csv
set -uo pipefail
DS="${1:?usage: run_metawrap.sh <dataset> <mock_dir>}"
M="${2:?mock dir with asm/final.contigs.fasta + reads_R1/R2.fastq}"

SIF="${METAWRAP_SIF:-/path/to/metawrap.sif}"
CHECKM_DB="${CHECKM_DB:-/path/to/checkm_db}"
RAPDTOOL_SIF="${RAPDTOOL_SIF:-/path/to/rapdtool_v2.3.0.sif}"
THREADS="${THREADS:-16}"
OUT="${METAWRAP_OUT:-/tmp/metawrap_$DS}"
TIME_BIN=/usr/bin/time

# results dir derived from this script's location (no hard-coded absolute path)
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"      # the benchmark-kit root
RES="${BENCH_ROOT:-$HERE}/results/metawrap"
SUMMARY="$RES/summary.csv"
mkdir -p "$OUT" "$RES"
[[ -f "$SUMMARY" ]] || echo "tool,dataset,step,wall_seconds,peak_rss_gb,exit_code" > "$SUMMARY"

# --- preflight ----------------------------------------------------------------------
for f in "$M/asm/final.contigs.fasta" "$M/reads_R1.fastq" "$M/reads_R2.fastq"; do
  [[ -f "$f" ]] || { echo "[metawrap] missing input: $f" >&2; exit 2; }
done
[[ -f "$SIF" ]]          || { echo "[metawrap] MetaWRAP SIF not found: $SIF" >&2; exit 2; }
[[ -f "$RAPDTOOL_SIF" ]] || { echo "[metawrap] RaPDTool SIF (miComplete) not found: $RAPDTOOL_SIF" >&2; exit 2; }
[[ -d "$CHECKM_DB/genome_tree" ]] || { echo "[metawrap] CheckM DB not ready at $CHECKM_DB" >&2; exit 1; }

# Bind mounts cover every path the containers must see (inputs, output, DBs) — portable,
# not hard-wired to a fixed mount.
MWBINDS="-B $M -B $OUT -B $CHECKM_DB"
APP="apptainer exec $MWBINDS --env LC_ALL=C --env CHECKM_DATA_PATH=$CHECKM_DB $SIF"

# MetaWRAP wants reads named *_1.fastq / *_2.fastq
ln -sf "$M/reads_R1.fastq" "$OUT/reads_1.fastq"
ln -sf "$M/reads_R2.fastq" "$OUT/reads_2.fastq"

run() {  # run <step> -- <cmd...> ; times it, appends a summary.csv row
  local step="$1"; shift 2
  local tlog="$OUT/${step}.time.txt"
  echo "[metawrap] $(date) $DS :: $step"
  $TIME_BIN -v -o "$tlog" "$@" > "$OUT/${step}.log" 2>&1
  local ec=$?
  local rss wall wsec rgb
  rss=$(awk -F': ' '/Maximum resident set size/{print $2}' "$tlog")
  wall=$(awk -F': ' '/Elapsed \(wall clock\)/{print $2}' "$tlog")
  wsec=$(awk -v t="$wall" 'BEGIN{n=split(t,a,":");s=(n==3)?a[1]*3600+a[2]*60+a[3]:(n==2)?a[1]*60+a[2]:a[1];printf "%.1f",s}')
  rgb=$(awk -v k="${rss:-0}" 'BEGIN{printf "%.2f",k/1048576}')
  echo "metawrap,$DS,$step,$wsec,$rgb,$ec" >> "$SUMMARY"
  return $ec
}

# --- 1. binning (3-binner ensemble) + refinement ------------------------------------
echo "[metawrap] binning (metabat2 + maxbin2 + concoct) on RaPDTool's assembly"
run binning -- $APP metawrap binning -o "$OUT/binning" -t "$THREADS" \
    -a "$M/asm/final.contigs.fasta" --metabat2 --maxbin2 --concoct \
    "$OUT/reads_1.fastq" "$OUT/reads_2.fastq" || { echo "[metawrap] binning failed"; exit 1; }

echo "[metawrap] bin_refinement (CheckM; keep >=50% complete, <10% contam)"
run refine -- $APP metawrap bin_refinement -o "$OUT/refined" -t "$THREADS" \
    -A "$OUT/binning/metabat2_bins" -B "$OUT/binning/maxbin2_bins" -C "$OUT/binning/concoct_bins" \
    -c 50 -x 10 || { echo "[metawrap] refinement failed"; exit 1; }

BINS="$OUT/refined/metawrap_50_10_bins"
cp "$OUT/refined/metawrap_50_10_bins.stats" "$RES/${DS}.checkm.stats"

# --- 2. common-evaluator re-scoring: miComplete/Bact105 inside RaPDTool's SIF ---------
# MetaWRAP's refinement optimises for CheckM; scoring both bin sets with the SAME tool
# RaPDTool uses (miComplete) is what makes the completeness/contamination numbers
# comparable. miComplete wants a <path>\t<ext> list and .fna extensions.
echo "[metawrap] re-scoring MetaWRAP bins with miComplete/Bact105 (RaPDTool's evaluator)"
MC="$OUT/micomplete"; mkdir -p "$MC"; : > "$MC/bins.tab"
for f in "$BINS"/*.fa; do
  ln -sf "$f" "$MC/$(basename "${f%.fa}").fna"
  printf '%s\tfna\n' "$MC/$(basename "${f%.fa}").fna" >> "$MC/bins.tab"
done
apptainer exec -B "$OUT" "$RAPDTOOL_SIF" bash -lc \
  "cd $MC && miComplete $MC/bins.tab --hmms Bact105 --threads $THREADS > $MC/micomplete.tab 2> $MC/micomplete.log"
cp "$MC/micomplete.tab" "$RES/${DS}.micomplete.tab"

# --- 3. report ----------------------------------------------------------------------
echo "[metawrap] $(date) done — $DS"
awk -F'\t' '$1!="Name"&&$1!~/^#/&&NF>=6{comp=$5*100;red=($6-1)*100;n++;c[n]=comp;ct[n]=red;
    if(comp>=90&&red<5)hq++}
  END{asort(c);asort(ct);mc=(n%2)?c[(n+1)/2]:(c[n/2]+c[n/2+1])/2;
      mt=(n%2)?ct[(n+1)/2]:(ct[n/2]+ct[n/2+1])/2;
      printf "  MAGs=%d  median completeness=%.1f%%  median contamination=%.1f%%  HQ=%d/%d\n",n,mc,mt,hq,n}' \
  "$RES/${DS}.micomplete.tab"
echo "  resources -> $SUMMARY"; grep ",$DS," "$SUMMARY"
