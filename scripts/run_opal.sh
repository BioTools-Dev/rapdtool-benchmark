#!/usr/bin/env bash
# Convert every profile to CAMI and run OPAL, for each depth-series dataset.
# All profiles (including the gold standard) MUST carry the same SampleID (-s mock) or
# OPAL silently skips them. Conversions run in rapdtool_bench (taxonkit); OPAL runs from
# its own env's python directly (avoids conda-activate issues).
set -uo pipefail
cd "$(dirname "$0")/.."   # repo root
[ -f config.sh ] && source config.sh
P2C="python3 scripts/profile2cami.py"
OPALPY="${OPAL_PY:-/path/to/envs/opal/bin/python}"
OPAL="${OPAL_BIN:-/path/to/envs/opal/bin/opal.py}"

for D in "$@"; do
  B="results/bench_$D"; PR="$B/profiles"; mkdir -p "$PR"
  GOLD="${MOCK_ROOT:-/path/to/mock_datasets}/mock_$D/mock_composition.tsv"
  echo "=========================================================="
  echo "[opal] $D"
  # convert() aborts the dataset if a conversion fails or yields an empty profile,
  # instead of letting OPAL run on a stale/missing file.
  convert() {  # convert <out> -- <profile2cami args...>
    local out="$1"; shift 2
    if ! $P2C "$@" -o "$out" 2>"$out.err" || [ ! -s "$out" ]; then
      echo "[opal] CONVERSION FAILED -> $out"; sed -n '1,3p' "$out.err"; return 1
    fi
    rm -f "$out.err"
  }
  convert "$B/gold_standard.profile"    -- "$GOLD" -f truth -s mock          || continue
  convert "$PR/rapdtool_screen.profile" -- "$B"/rapdtool_screen.rep1/profilesfmbm/*/output_All_levels.csv -f focus -s mock || continue
  convert "$PR/rapdtool_full.profile"   -- "$B"/rapdtool.rep1/profilesfmbm/*/output_All_levels.csv        -f focus -s mock || continue
  for k in full cap16 cap8; do
    convert "$PR/kraken2_${k}.profile"  -- "$B/kraken2_${k}.rep1.bracken" -f bracken -s mock || continue 2
  done
  convert "$PR/metaphlan.profile"       -- "$B/metaphlan.rep1.profile" -f metaphlan -s mock || continue

  "$OPALPY" "$OPAL" -g "$B/gold_standard.profile" -o "$B/opal" \
    "$PR/rapdtool_screen.profile" "$PR/rapdtool_full.profile" \
    "$PR/kraken2_full.profile" "$PR/kraken2_cap16.profile" "$PR/kraken2_cap8.profile" \
    "$PR/metaphlan.profile" \
    -l "RaPDTool_screen,RaPDTool_full,Kraken2_full,Kraken2_16GB,Kraken2_8GB,MetaPhlAn4" \
    2>&1 | tail -2
done
echo "[opal] done"
