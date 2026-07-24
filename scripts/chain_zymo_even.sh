#!/usr/bin/env bash
# ZymoBIOMICS even community (D6300, ERR2984773) end-to-end:
# NOTE: example orchestration used for this study. Set the dataset paths via config.sh
# (MOCK_ROOT, ZYMO_DIR) and `source config.sh` first (KRAKEN_DB_FULL, MPA_DB come from it).
#   decompress -> assemble (for RaPDTool full / MAGs on real data) -> benchmark.
# Kraken full + MetaPhlAn + RaPDTool full+screen (capped Kraken skipped: the DB-size
# curve is already on the mocks; Zymo is the real-data competence test).
set -uo pipefail
cd "$(dirname "$0")/.."   # repo root
[ -f config.sh ] && source config.sh
ENVBIN=$HOME/miniconda3/envs/rapdtool_bench/bin
export PATH="$ENVBIN:$PATH"

Z=${ZYMO_DIR:?source config.sh}
M=${MOCK_ROOT:?source config.sh}/zymo_even
mkdir -p "$M"

echo "[chain] $(date) decompressing reads"
[ -s "$M/reads_R1.fastq" ] || zcat "$Z/ERR2984773_1.fastq.gz" > "$M/reads_R1.fastq"
[ -s "$M/reads_R2.fastq" ] || zcat "$Z/ERR2984773_2.fastq.gz" > "$M/reads_R2.fastq"
echo "[chain] pairs: $(( $(wc -l < "$M/reads_R1.fastq") / 4 ))"

echo "[chain] $(date) assembling (MEGAHIT)"
if [ ! -s "$M/asm/final.contigs.fasta" ]; then
  rm -rf "$M/asm"
  "$ENVBIN/megahit" -1 "$M/reads_R1.fastq" -2 "$M/reads_R2.fastq" -t 16 -o "$M/asm" \
    > "$M/megahit.log" 2>&1
  cp "$M/asm/final.contigs.fa" "$M/asm/final.contigs.fasta"   # FOCUS needs .fasta
fi
echo "[chain] assembly: $(grep -v '^>' "$M/asm/final.contigs.fasta" | tr -d '\n' | wc -c) bp, $(grep -c '^>' "$M/asm/final.contigs.fasta") contigs"

echo "[chain] $(date) running benchmark"
BENCH_ENV_BIN=$ENVBIN \
RAPDTOOL=${RAPDTOOL:?source config.sh} \
THREADS=16 REPEATS=1 OUTDIR=$PWD/results/bench_zymo_even \
RUN_KRAKEN_CAP16=0 RUN_KRAKEN_CAP8=0 \
ASSEMBLY=$M/asm/final.contigs.fasta \
READS_R1=$M/reads_R1.fastq READS_R2=$M/reads_R2.fastq \
KRAKEN_DB_FULL=${KRAKEN_DB_FULL:?source config.sh} \
MPA_DB=${MPA_DB:?source config.sh} \
scripts/run_benchmark.sh

echo "[chain] $(date) done"
cat bench_zymo_even/summary.csv 2>/dev/null
