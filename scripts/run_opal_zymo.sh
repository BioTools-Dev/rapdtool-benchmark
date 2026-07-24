#!/usr/bin/env bash
# OPAL for a ZymoBIOMICS community, respecting the sequence- vs cell-abundance split.
#
# Kraken2/Bracken and FOCUS estimate SEQUENCE abundance     -> gold_dna.profile
# MetaPhlAn estimates CELL abundance (genome-size normalised) -> gold_genomecopy.profile
# Detection (recall/precision/F1 of which species) is basis-independent; only the
# abundance-error metrics (L1, Bray-Curtis) differ, so each tool is scored against its
# matching basis. Run OPAL twice, once per basis, then read the appropriate rows.
#
# Usage: ./run_opal_zymo.sh bench_zymo_even
set -uo pipefail
cd "$(dirname "$0")/.."   # repo root
[ -f config.sh ] && source config.sh
B="${1:?usage: run_opal_zymo.sh <bench_dir>}"
PR="$B/profiles"; mkdir -p "$PR"
P2C="python3 scripts/profile2cami.py"
OPALPY="${OPAL_PY:-/path/to/envs/opal/bin/python}"
OPAL="${OPAL_BIN:-/path/to/envs/opal/bin/opal.py}"

convert() { local out="$1"; shift 2
  if ! $P2C "$@" -o "$out" 2>"$out.err" || [ ! -s "$out" ]; then
    echo "CONVERSION FAILED -> $out"; sed -n '1,3p' "$out.err"; exit 1; fi
  rm -f "$out.err"; }

convert "$PR/rapdtool_screen.profile" -- "$B"/rapdtool_screen.rep1/profilesfmbm/*/output_All_levels.csv -f focus -s mock
convert "$PR/rapdtool_full.profile"   -- "$B"/rapdtool.rep1/profilesfmbm/*/output_All_levels.csv        -f focus -s mock
convert "$PR/kraken2_full.profile"    -- "$B/kraken2_full.rep1.bracken" -f bracken   -s mock
convert "$PR/metaphlan.profile"       -- "$B/metaphlan.rep1.profile"    -f metaphlan -s mock

# Sequence-abundance tools vs the DNA-basis gold
"$OPALPY" "$OPAL" -g "$B/gold_dna.profile" -o "$B/opal_dna" \
  "$PR/rapdtool_screen.profile" "$PR/rapdtool_full.profile" "$PR/kraken2_full.profile" \
  -l "RaPDTool_screen,RaPDTool_full,Kraken2_full" 2>&1 | tail -1
# Cell-abundance tool vs the genome-copy-basis gold
"$OPALPY" "$OPAL" -g "$B/gold_genomecopy.profile" -o "$B/opal_genomecopy" \
  "$PR/metaphlan.profile" -l "MetaPhlAn4" 2>&1 | tail -1
echo "[zymo-opal] done: $B/opal_dna (seq-abundance tools) + $B/opal_genomecopy (MetaPhlAn)"
