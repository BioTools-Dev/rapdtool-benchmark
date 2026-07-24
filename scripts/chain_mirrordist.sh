#!/usr/bin/env bash
# Wait for the distance-mirror mock to finish, verify it, then run the benchmark.
# NOTE: example orchestration used for this study. Set the dataset paths via config.sh
# (MOCK_ROOT) and `source config.sh` first (KRAKEN_DB_FULL, MPA_DB come from it too).
# Runs unattended; every step is logged to chain_mirrordist.log.
set -uo pipefail
cd "$(dirname "$0")/.."   # repo root
[ -f config.sh ] && source config.sh

MOCKLOG=$MOCK_ROOT/mock_mirrordist_15M.log
ASM=$MOCK_ROOT/mock_mirrordist_15M/asm/final.contigs.fasta
COMP=$MOCK_ROOT/mock_mirrordist_15M/mock_composition.tsv
M=$MOCK_ROOT/mock_mirrordist_15M

echo "[chain] $(date) waiting for mock to finish (timeout ~45 min)"
ok=0
for i in $(seq 1 270); do          # 270 * 10s = 45 min
  if grep -q 'make_mock] done' "$MOCKLOG" 2>/dev/null && [ -s "$ASM" ]; then ok=1; break; fi
  sleep 10
done

if [ "$ok" != 1 ] || [ ! -s "$ASM" ]; then
  echo "[chain] MOCK DID NOT COMPLETE — aborting, not running benchmark"; exit 1
fi

ntax=$(grep -vc '^#' "$COMP" 2>/dev/null || echo 0)
asmbp=$(grep -v '^>' "$ASM" | tr -d '\n' | wc -c)
echo "[chain] mock done: gold-standard taxids=$ntax  assembly=$asmbp bp"
if [ "$ntax" -lt 10 ]; then
  echo "[chain] gold standard has <10 taxids — aborting"; exit 1
fi

echo "[chain] $(date) launching benchmark (RaPDTool full+screen + Kraken control)"
BENCH_ENV_BIN=$HOME/miniconda3/envs/rapdtool_bench/bin \
RAPDTOOL=${RAPDTOOL:?source config.sh} \
THREADS=16 REPEATS=1 OUTDIR=$PWD/results/bench_mirrordist \
RUN_KRAKEN_CAP16=0 RUN_KRAKEN_CAP8=0 RUN_METAPHLAN=0 \
ASSEMBLY=$M/asm/final.contigs.fasta \
READS_R1=$M/reads_R1.fastq READS_R2=$M/reads_R2.fastq \
KRAKEN_DB_FULL=${KRAKEN_DB_FULL:?source config.sh} \
scripts/run_benchmark.sh

rc=$?
echo "[chain] $(date) benchmark finished (exit $rc)"
echo "[chain] summary:"
cat bench_mirrordist/summary.csv 2>/dev/null
exit $rc
