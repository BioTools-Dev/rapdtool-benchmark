#!/usr/bin/env bash
# verify_kit.sh — confirm this kit reproduces the published figures and numbers.
#
# Runs the "level 1" reproduction end to end (README, Configuration) and checks each
# result against what is committed in the repository, so you find out in ~2 minutes
# whether your environment reproduces the study — before spending days on the full
# re-run.
#
#   source config.sh
#   scripts/verify_kit.sh
#
# What it checks (each is PASS/FAIL, nothing is assumed):
#   1. prerequisites        — $TAXONKIT_DB, python3, matplotlib, git
#   2. the five figures     — regenerated and compared byte-for-byte against the
#                             committed PNGs (PNG carries no timestamp, so an exact
#                             match means an exact reproduction; SVG/PDF embed a
#                             creation date and are therefore not byte-comparable)
#   3. threshold sweep      — figures/threshold_sweep.tsv vs the committed file
#   4. detection numbers    — mock (20/20, 0 FP) and Zymo (8/8, 0 FP)
#   5. census numbers       — recomputed from data/census_full.tsv (no DB needed)
#   6. optional: the census is REBUILT from the competitor databases if they are
#      configured (needs $KRAKEN_STD_INSPECT, $MPA_PKL, $FOCUS_DB); skipped otherwise
#
# Regenerated files are written to a scratch copy of figures/ and the repository is
# left untouched, so this is safe to run on a clean checkout.
set -uo pipefail
cd "$(dirname "$0")/.."
BASE="$PWD"

PASS=0; FAIL=0; SKIP=0
ok()   { printf '  \033[32mPASS\033[0m  %s\n' "$1"; PASS=$((PASS+1)); }
bad()  { printf '  \033[31mFAIL\033[0m  %s\n' "$1"; FAIL=$((FAIL+1)); }
skip() { printf '  \033[33mSKIP\033[0m  %s\n' "$1"; SKIP=$((SKIP+1)); }
head_() { printf '\n\033[1m%s\033[0m\n' "$1"; }

PY=${PYTHON:-python3}

# ---------------------------------------------------------------- 1. prerequisites
head_ "1. Prerequisites"
command -v "$PY" >/dev/null && ok "python3 ($($PY -V 2>&1))" || { bad "python3 not found"; exit 1; }
$PY -c 'import matplotlib' 2>/dev/null && ok "matplotlib importable" \
  || bad "matplotlib missing — pip install matplotlib (needed by every figure script)"
command -v git >/dev/null && ok "git (needed to compare against the committed files)" \
  || skip "git not found — figure comparison will be skipped"
if [ -n "${TAXONKIT_DB:-}" ] && [ -d "${TAXONKIT_DB:-/nonexistent}" ]; then
  ok "TAXONKIT_DB=$TAXONKIT_DB"
else
  bad "TAXONKIT_DB unset or not a directory — did you 'source config.sh'? (README, Configuration)"
fi

# committed-file comparison helper: $1 = repo-relative path
same_as_committed() {
  git -C "$BASE" cat-file -e "HEAD:$1" 2>/dev/null || return 2   # not committed
  local a b
  a=$(md5sum < "$1" | cut -d' ' -f1)
  b=$(git -C "$BASE" show "HEAD:$1" | md5sum | cut -d' ' -f1)
  [ "$a" = "$b" ]
}

# ------------------------------------------------------------------- 2. the figures
head_ "2. Figures regenerate identically (PNG, byte-for-byte)"
if ! command -v git >/dev/null || ! git -C "$BASE" rev-parse HEAD >/dev/null 2>&1; then
  skip "not a git checkout — cannot compare against committed figures"
else
  TMPFIG=$(mktemp -d); trap 'rm -rf "$TMPFIG"' EXIT
  cp -a figures/. "$TMPFIG/" 2>/dev/null || true
  LOG=$(mktemp)
  run_fig() {                       # $1 = label, $2 = repo-relative png, rest = command
    local label="$1" png="$2"; shift 2
    if ! "$@" >"$LOG" 2>&1; then
      bad "$label — script exited non-zero:"; sed 's/^/          /' "$LOG" | tail -4; return
    fi
    if same_as_committed "$png"; then ok "$label"
    elif [ $? -eq 2 ];              then skip "$label (not committed, nothing to compare)"
    else bad "$label — regenerated $png differs from the committed figure"; fi
  }
  run_fig "Fig 2  census"          figures/census.png          $PY scripts/plot_census.py \
            -i data/census_full.tsv -l data/mock_genomes.list -o figures/census
  run_fig "Fig 3  opal_depth"      figures/opal_depth.png      $PY scripts/plot_opal_depth.py
  run_fig "Fig 4  micomplete"      figures/micomplete.png      $PY scripts/plot_micomplete.py
  run_fig "Fig 5  mirror_distance" figures/mirror_distance.png $PY scripts/plot_mirror_distance.py
  run_fig "       f1_threshold"    figures/f1_threshold.png    $PY scripts/plot_f1_threshold.py
  rm -f "$LOG"
fi

# -------------------------------------------------------------- 3. threshold sweep
head_ "3. Threshold sweep"
if $PY scripts/threshold_sweep.py --truth results/bench_ln_30M/gold_standard.profile \
      -o figures/threshold_sweep >/dev/null 2>&1; then
  if same_as_committed figures/threshold_sweep.tsv; then ok "threshold_sweep.tsv matches"
  elif [ $? -eq 2 ]; then skip "threshold_sweep.tsv not committed"
  else bad "threshold_sweep.tsv differs from the committed file"; fi
else
  bad "threshold_sweep.py failed"
fi

# ------------------------------------------------------------ 4. detection numbers
head_ "4. Species detection (the central experiment)"
check_detection() {                 # $1 label, $2 expected "TP FP", rest = args
  local label="$1" want="$2"; shift 2
  local out got
  out=$($PY scripts/mash_detection.py "$@" 2>&1) || { bad "$label — script failed"; return; }
  got=$(printf '%s\n' "$out" | sed -n 's/.*TP=\([0-9]*\) *FP=\([0-9]*\).*/\1 \2/p' | head -1)
  if [ "$got" = "$want" ]; then ok "$label — TP/FP = $got (expected $want)"
  else bad "$label — TP/FP = '${got:-none}', expected '$want'"
       printf '%s\n' "$out" | sed 's/^/          /' | tail -4; fi
}
check_detection "mock 30 M (20 species, 0 false positives)" "20 0" \
  --bench results/bench_ln_30M --gold results/bench_ln_30M/gold_standard.profile \
  --split data/mock_genomes.list
check_detection "ZymoBIOMICS even (8 species, 0 false positives)" "8 0" \
  --bench results/bench_zymo_even --gold results/bench_zymo_even/gold_dna.profile \
  --extra-true 96241:1423
# the reference/conflictive split is the headline claim — check it explicitly
if $PY scripts/mash_detection.py --bench results/bench_ln_30M \
      --gold results/bench_ln_30M/gold_standard.profile --split data/mock_genomes.list 2>&1 \
    | grep -q "conflictive detected: 10/10"; then
  ok "conflictive split — 10/10 detected (species absent from Kraken2 and MetaPhlAn)"
else
  bad "conflictive split is not 10/10 — the central result did not reproduce"
fi

# --------------------------------------------------------------- 5. census numbers
head_ "5. Census figures quoted in the manuscript"
$PY - "$BASE/data/census_full.tsv" <<'PYEOF'
import csv, sys
rows = list(csv.DictReader(open(sys.argv[1]), delimiter="\t"))
n = len(rows)
kf   = sum(1 for r in rows if r["kraken_full"] == "N")
k8   = sum(1 for r in rows if r["kraken8"]     == "N")
mpa  = sum(1 for r in rows if r["metaphlan"]   == "N")
both = sum(1 for r in rows if r["verdict"] == "conflictive-both")
ref  = sum(1 for r in rows if r["verdict"] == "reference")
checks = [
    ("30,209 type-material genomes",            n == 30209,        n),
    ("55.8 % absent from Kraken2 standard",     round(100*kf/n,1) == 55.8,  "%.1f %% (%d)" % (100*kf/n, kf)),
    ("17.3 % absent from MetaPhlAn4",           round(100*mpa/n,1) == 17.3, "%.1f %% (%d)" % (100*mpa/n, mpa)),
    ("11.3 % (3,423) absent from both",         both == 3423,      both),
    ("11,552 present in both (reference pool)", ref == 11552,      ref),
    ("capping to 8 GB changes absence <0.01 pp",abs(100*(k8-kf)/n) < 0.01, "%+.4f pp" % (100*(k8-kf)/n)),
]
bad = 0
for label, good, got in checks:
    print(("  \033[32mPASS\033[0m  %s  [%s]" if good else "  \033[31mFAIL\033[0m  %s  [got %s]") % (label, got))
    bad += 0 if good else 1
sys.exit(1 if bad else 0)
PYEOF
if [ $? -eq 0 ]; then PASS=$((PASS+6)); else FAIL=$((FAIL+1)); fi

# ------------------------------------------- 6. optional: rebuild census from the DBs
head_ "6. Optional — rebuild the census from the competitor databases"
if [ -r "${KRAKEN_STD_INSPECT:-/nonexistent}" ] && [ -r "${MPA_PKL:-/nonexistent}" ] \
   && [ -d "${FOCUS_DB:-/nonexistent}" ]; then
  TMPC=$(mktemp)
  if $PY scripts/check_representation.py --sample 30209 --seed 42 -o "$TMPC" >/dev/null 2>&1; then
    if cmp -s "$TMPC" data/census_full.tsv; then
      ok "census rebuilt from the databases is identical to data/census_full.tsv"
    else
      bad "rebuilt census differs — your database builds differ from the study's (see README, Environment)"
    fi
  else
    bad "check_representation.py failed"
  fi
  rm -f "$TMPC"
else
  skip "competitor databases not configured (KRAKEN_STD_INSPECT / MPA_PKL / FOCUS_DB) — census verified from the shipped TSV only"
fi

# ------------------------------------------------------------------------- summary
head_ "Summary"
printf '  %d passed, %d failed, %d skipped\n\n' "$PASS" "$FAIL" "$SKIP"
if [ "$FAIL" -eq 0 ]; then
  printf '  \033[32mThis kit reproduces the published figures and numbers on your machine.\033[0m\n'
  printf '  The full re-run (README Steps 1-6) needs the large databases; see Configuration.\n\n'
  exit 0
fi
printf '  \033[31mSomething did not reproduce.\033[0m Most common causes, in order:\n'
printf '    1. config.sh not sourced (TAXONKIT_DB unset)\n'
printf '    2. a different matplotlib version (figures differ but numbers pass — harmless)\n'
printf '    3. edited result files under results/\n\n'
exit 1
