#!/usr/bin/env bash
#
# run_benchmark.sh — reproducible resource + accuracy benchmark for the RaPDTool paper.
#
# Measures, for each tool, on the SAME input and hardware:
#   - on-disk database size (du)
#   - peak resident memory + wall-clock time (/usr/bin/time -v)
# and lays out the folders for the accuracy step (OPAL / AMBER), which you run after.
#
# It does NOT invent any numbers: it only records what the tools actually do.
# Edit the CONFIG block, then run.  Re-run REPEATS times and take the median.
#
# Requirements on PATH: /usr/bin/time (GNU), du, and whichever tools you enable below
# (rapdtool, kraken2, bracken, metaphlan, motus).  OPAL/AMBER are used in the last step.

set -uo pipefail

# ------------------------------- CONFIG --------------------------------------
THREADS="${THREADS:-16}"                 # keep identical across all tools
REPEATS="${REPEATS:-3}"                  # runs per tool; report the median
OUTDIR="${OUTDIR:-$PWD/results/bench_out}"       # results root
SUMMARY="$OUTDIR/summary.csv"

# Inputs (edit these) -------------------------------------------------------
ASSEMBLY="${ASSEMBLY:-/path/to/assembly.fasta}"   # for RaPDTool (contigs)
READS_R1="${READS_R1:-/path/to/reads_R1.fastq.gz}"  # for read/marker classifiers
READS_R2="${READS_R2:-/path/to/reads_R2.fastq.gz}"

# Databases (edit these) ----------------------------------------------------
KRAKEN_DB_FULL="${KRAKEN_DB_FULL:-/path/to/kraken2_standard}"    # accuracy ceiling
KRAKEN_DB_CAP16="${KRAKEN_DB_CAP16:-/path/to/kraken2_16GB}"      # capped, memory-limited
KRAKEN_DB_CAP8="${KRAKEN_DB_CAP8:-/path/to/kraken2_8GB}"         # capped, closest to RaPDTool
MPA_DB="${MPA_DB:-}"                      # MetaPhlAn --db_dir (was --bowtie2db; empty = default)
READ_LEN="${READ_LEN:-150}"              # for Bracken

# Tool launchers (edit if not on PATH) --------------------------------------
RAPDTOOL="${RAPDTOOL:-rapdtool}"         # RaPDTool launcher. RaPDTool ships as an Apptainer
                                         # wrapper installed in its OWN conda env ('rapdtool').
                                         # Preferred: the absolute wrapper path
                                         # (from $RAPDTOOL, set in config.sh) run from this
                                         # env — apptainer is on PATH here (/usr/bin/apptainer) and
                                         # the wrapper exec's it, so /usr/bin/time measures RaPDTool's
                                         # real peak RSS.
RAPDTOOL_ENV="${RAPDTOOL_ENV:-}"         # Optional: run RaPDTool via 'conda run -n <env>' (e.g.
                                         # RAPDTOOL_ENV=rapdtool). Use only if the wrapper won't run
                                         # here. NOTE: conda run inserts a persistent python parent,
                                         # so /usr/bin/time may under-report RaPDTool's peak RSS.

# The other tools default to bare names (resolved on PATH), but PATH resolution has
# already silently broken one full run: if the conda env is not active -- or if PATH
# contains relative entries such as '.', './bin' (as this host's .zshrc does), so that
# resolution depends on the current directory -- every tool exits 126/127 in ~0 s and
# the matrix completes with no data. Pass ABSOLUTE paths to make the run independent of
# shell state; BENCH_ENV_BIN sets them all at once:
#   BENCH_ENV_BIN=$HOME/miniconda3/envs/rapdtool_bench/bin ./run_benchmark.sh
BENCH_ENV_BIN="${BENCH_ENV_BIN:-}"
KRAKEN2="${KRAKEN2:-${BENCH_ENV_BIN:+$BENCH_ENV_BIN/}kraken2}"
BRACKEN="${BRACKEN:-${BENCH_ENV_BIN:+$BENCH_ENV_BIN/}bracken}"
METAPHLAN="${METAPHLAN:-${BENCH_ENV_BIN:+$BENCH_ENV_BIN/}metaphlan}"

# Absolute launchers are not sufficient on their own: several tools shell out to
# SIBLING executables that they resolve on PATH. MetaPhlAn invokes 'bowtie2' and dies
# with "[Errno 13] Permission denied: 'bowtie2'" if it is not resolvable, after having
# already loaded ~7 GB -- i.e. it fails late and looks like a different problem.
# Prepending the env's bin to PATH fixes every such sibling lookup at once, while the
# absolute launchers above still pin WHICH tool we measure.
if [[ -n "$BENCH_ENV_BIN" ]]; then
  export PATH="$BENCH_ENV_BIN:$PATH"
fi
# RaPDTool uses its own cached DBs (mash + FOCUS); size them via rapdtool --where.
# full + cap16 + cap8 give the Kraken accuracy-vs-DB-size curve against RaPDTool's <0.5 GB.

# Toggle tools on/off (all overridable from the environment: RUN_X=0 ./run_benchmark.sh) --
RUN_RAPDTOOL_FULL="${RUN_RAPDTOOL_FULL:-1}"     # full mode (assembly): binning + per-bin classification + MAGs
RUN_RAPDTOOL_SCREEN="${RUN_RAPDTOOL_SCREEN:-1}" # screen mode (reads): mash-screen + FOCUS, read-based profiling
                                                # — apples-to-apples input vs Kraken/MetaPhlAn for OPAL
RUN_KRAKEN_FULL="${RUN_KRAKEN_FULL:-1}"
RUN_KRAKEN_CAP16="${RUN_KRAKEN_CAP16:-1}"
RUN_KRAKEN_CAP8="${RUN_KRAKEN_CAP8:-1}"
RUN_METAPHLAN="${RUN_METAPHLAN:-1}"
RUN_MOTUS="${RUN_MOTUS:-0}"
# -----------------------------------------------------------------------------

mkdir -p "$OUTDIR"
if [[ ! -f "$SUMMARY" ]]; then
  echo "tool,replicate,wall_seconds,peak_rss_kb,peak_rss_gb,exit_code" > "$SUMMARY"
fi

# --- Preflight: fail NOW, not six hours from now ------------------------------
# A run once completed its whole matrix in minutes because kraken2/metaphlan were not
# resolvable, exiting 126 in 0.0 s per tool. Every enabled tool is checked up front.
preflight_fail=0
check_tool() {   # check_tool <enabled> <label> <launcher...>
  [[ "$1" == 1 ]] || return 0
  local label="$2"; shift 2
  if ! command -v "$1" >/dev/null 2>&1 && [[ ! -x "$1" ]]; then
    echo "PREFLIGHT: $label -> '$1' not found or not executable" >&2
    preflight_fail=1
  fi
}
check_input() {  # check_input <enabled> <label> <path>
  [[ "$1" == 1 ]] || return 0
  [[ -e "$3" ]] || { echo "PREFLIGHT: $2 missing: $3" >&2; preflight_fail=1; }
}

check_tool "$RUN_KRAKEN_FULL$RUN_KRAKEN_CAP16$RUN_KRAKEN_CAP8" kraken2 "$KRAKEN2" 2>/dev/null
[[ "$RUN_KRAKEN_FULL" == 1 || "$RUN_KRAKEN_CAP16" == 1 || "$RUN_KRAKEN_CAP8" == 1 ]] && {
  check_tool 1 kraken2 "$KRAKEN2"; check_tool 1 bracken "$BRACKEN"
  check_input "$RUN_KRAKEN_FULL"  "kraken DB (full)"  "$KRAKEN_DB_FULL"
  check_input "$RUN_KRAKEN_CAP16" "kraken DB (cap16)" "$KRAKEN_DB_CAP16"
  check_input "$RUN_KRAKEN_CAP8"  "kraken DB (cap8)"  "$KRAKEN_DB_CAP8"
}
check_tool "$RUN_METAPHLAN" metaphlan "$METAPHLAN"
# MetaPhlAn shells out to bowtie2 by name; check the sibling too, or it fails ~30 s in.
check_tool "$RUN_METAPHLAN" "bowtie2 (MetaPhlAn dependency)" bowtie2
[[ "$RUN_METAPHLAN" == 1 && -n "$MPA_DB" ]] && check_input 1 "MetaPhlAn DB" "$MPA_DB"
check_tool "$RUN_RAPDTOOL_FULL$RUN_RAPDTOOL_SCREEN" rapdtool "$RAPDTOOL" 2>/dev/null
[[ "$RUN_RAPDTOOL_FULL" == 1 || "$RUN_RAPDTOOL_SCREEN" == 1 ]] && {
  [[ -n "$RAPDTOOL_ENV" ]] || check_tool 1 rapdtool "$RAPDTOOL"
}
check_input "$RUN_RAPDTOOL_FULL" "assembly" "$ASSEMBLY"
[[ "$RUN_KRAKEN_FULL$RUN_KRAKEN_CAP16$RUN_KRAKEN_CAP8$RUN_METAPHLAN$RUN_RAPDTOOL_SCREEN" == *1* ]] && {
  check_input 1 "reads R1" "$READS_R1"; check_input 1 "reads R2" "$READS_R2"
}

if [[ "$preflight_fail" == 1 ]]; then
  cat >&2 <<'EOF'

ABORTING before any tool runs.

Most often the conda env is not active, or PATH resolution is unreliable (this host's
.zshrc puts '.', './bin', './scripts' first, so resolution depends on the current
directory). Pass absolute paths instead of relying on PATH:

  BENCH_ENV_BIN=$HOME/miniconda3/envs/rapdtool_bench/bin \
  RAPDTOOL=$RAPDTOOL \
  ... ./run_benchmark.sh
EOF
  exit 1
fi
echo "[preflight] all enabled tools and inputs present" >&2

# Portable GNU time locator
TIME_BIN="$(command -v time)"; [[ -x /usr/bin/time ]] && TIME_BIN=/usr/bin/time
if ! "$TIME_BIN" -v true >/dev/null 2>&1; then
  echo "ERROR: GNU /usr/bin/time (with -v) is required." >&2; exit 1
fi

# run_timed <tool-label> <replicate> -- <command...>
run_timed() {
  local label="$1" rep="$2"; shift 3   # drop label, rep, and the literal '--'
  local tlog="$OUTDIR/${label}.rep${rep}.time.txt"
  echo ">>> [$label] replicate $rep: $*" >&2
  "$TIME_BIN" -v -o "$tlog" "$@" > "$OUTDIR/${label}.rep${rep}.stdout" \
                                 2> "$OUTDIR/${label}.rep${rep}.stderr"
  local ec=$?
  # Parse GNU time fields
  local rss wall
  rss=$(awk -F': ' '/Maximum resident set size/{print $2}' "$tlog")
  wall=$(awk -F': ' '/Elapsed \(wall clock\)/{print $2}' "$tlog")
  # wall may be h:mm:ss or m:ss.ss -> seconds
  local wsec
  wsec=$(awk -v t="$wall" 'BEGIN{
    n=split(t,a,":"); s=0;
    if(n==3){s=a[1]*3600+a[2]*60+a[3]} else if(n==2){s=a[1]*60+a[2]} else {s=a[1]}
    printf "%.1f", s }')
  local rgb; rgb=$(awk -v k="${rss:-0}" 'BEGIN{printf "%.2f", k/1048576}')
  echo "$label,$rep,$wsec,${rss:-NA},$rgb,$ec" >> "$SUMMARY"
}

db_size() {  # db_size <label> <path>  -> GB with 2 decimals (from bytes)
  [[ -e "$2" ]] || return 0
  local b; b=$(du -sb "$2" 2>/dev/null | cut -f1)
  awk -v n="${b:-0}" -v l="$1" \
      'BEGIN{printf "db_size,%s,%.2f GB\n", l, n/1073741824}' >> "$OUTDIR/db_sizes.csv"
}

run_kraken() {  # run_kraken <label> <db> <rep>  -> timed kraken2 + (untimed) Bracken
  local label="$1" db="$2" rep="$3"
  run_timed "$label" "$rep" -- \
    "$KRAKEN2" --db "$db" --threads "$THREADS" --paired \
            --report "$OUTDIR/${label}.rep${rep}.report" \
            --output "$OUTDIR/${label}.rep${rep}.kraken" "$READS_R1" "$READS_R2"
  # Bracken is fast post-processing (low RAM); timed separately from the classifier
  # footprint, which is dominated by Kraken2 loading its database into memory.
  "$BRACKEN" -d "$db" -i "$OUTDIR/${label}.rep${rep}.report" \
          -o "$OUTDIR/${label}.rep${rep}.bracken" -r "$READ_LEN" -l S \
          > "$OUTDIR/${label}.rep${rep}.bracken.log" 2>&1 || true
}

# RaPDTool launcher as an array (optional 'conda run -n <env>' prefix) --------
RAPD=()
[[ -n "$RAPDTOOL_ENV" ]] && RAPD=(conda run --no-capture-output -n "$RAPDTOOL_ENV")
RAPD+=("$RAPDTOOL")

# --- Database sizes (once) ---------------------------------------------------
: > "$OUTDIR/db_sizes.csv"
db_size kraken2_full  "$KRAKEN_DB_FULL"
db_size kraken2_cap16 "$KRAKEN_DB_CAP16"
db_size kraken2_cap8  "$KRAKEN_DB_CAP8"
"${RAPD[@]}" --where > "$OUTDIR/rapdtool_where.txt" 2>&1 || true   # record RaPDTool DB paths → du them

# RaPDTool screen mode profiles reads with a single -i; hand it the SAME data as the read
# classifiers by concatenating R1+R2 once. This is untimed input prep (analogous to giving
# the two read files to Kraken/MetaPhlAn) — only the rapdtool_screen run itself is timed.
SCREEN_READS="$OUTDIR/screen_reads.fastq"
if [[ "$RUN_RAPDTOOL_SCREEN" == 1 && ! -s "$SCREEN_READS" ]]; then
  cat "$READS_R1" "$READS_R2" > "$SCREEN_READS"
fi

# --- Timed runs --------------------------------------------------------------
for rep in $(seq 1 "$REPEATS"); do

  if [[ "$RUN_RAPDTOOL_FULL" == 1 ]]; then
    run_timed rapdtool "$rep" -- \
      "${RAPD[@]}" -i "$ASSEMBLY" -o "$OUTDIR/rapdtool.rep${rep}" -m full -t "$THREADS" --force
  fi

  if [[ "$RUN_RAPDTOOL_SCREEN" == 1 ]]; then
    run_timed rapdtool_screen "$rep" -- \
      "${RAPD[@]}" -i "$SCREEN_READS" -o "$OUTDIR/rapdtool_screen.rep${rep}" -m screen -t "$THREADS" --force
  fi

  [[ "$RUN_KRAKEN_FULL"  == 1 ]] && run_kraken kraken2_full  "$KRAKEN_DB_FULL"  "$rep"
  [[ "$RUN_KRAKEN_CAP16" == 1 ]] && run_kraken kraken2_cap16 "$KRAKEN_DB_CAP16" "$rep"
  [[ "$RUN_KRAKEN_CAP8"  == 1 ]] && run_kraken kraken2_cap8  "$KRAKEN_DB_CAP8"  "$rep"

  if [[ "$RUN_METAPHLAN" == 1 ]]; then
    mpa_args=("$METAPHLAN" "$READS_R1","$READS_R2" --input_type fastq --nproc "$THREADS"
              --mapout "$OUTDIR/metaphlan.rep${rep}.mapout"    # MetaPhlAn 4.1+: --mapout (was --bowtie2out)
              -o "$OUTDIR/metaphlan.rep${rep}.profile")
    [[ -n "$MPA_DB" ]] && mpa_args+=(--db_dir "$MPA_DB")   # MetaPhlAn 4.1+: --db_dir (was --bowtie2db)
    run_timed metaphlan "$rep" -- "${mpa_args[@]}"
  fi

  if [[ "$RUN_MOTUS" == 1 ]]; then
    run_timed motus "$rep" -- \
      motus profile -f "$READS_R1" -r "$READS_R2" -t "$THREADS" \
                    -o "$OUTDIR/motus.rep${rep}.profile"
  fi

done

# --- Median summary ----------------------------------------------------------
echo
echo "=== Median wall-clock (s) and peak RSS (GB) per tool ==="
awk -F, 'NR>1 && $6==0 {
           w[$1]=w[$1]" "$3; g[$1]=g[$1]" "$5 }
         END{
           for(t in w){
             n=split(w[t],A," ");   # bubble-sort the small arrays, then take the middle
             for(i=1;i<=n;i++)for(j=i+1;j<=n;j++){if(A[j]+0<A[i]+0){x=A[i];A[i]=A[j];A[j]=x}}
             m=split(g[t],B," ");
             for(i=1;i<=m;i++)for(j=i+1;j<=m;j++){if(B[j]+0<B[i]+0){x=B[i];B[i]=B[j];B[j]=x}}
             printf "  %-14s  wall=%ss  RSS=%sGB\n", t, A[int((n+1)/2)], B[int((m+1)/2)]
           }
         }' "$SUMMARY"

cat <<EOF

Raw per-replicate metrics : $SUMMARY
Database sizes            : $OUTDIR/db_sizes.csv   (add RaPDTool DBs from rapdtool_where.txt via 'du -sBG')
Per-run logs/outputs      : $OUTDIR/<tool>.repN.*

------------------------------------------------------------------------------
NEXT (accuracy — run manually once profiles exist):

0) Pin the NCBI taxonomy (reproducible name/taxid resolution; ensures the updated
   phylum names from FOCUS resolve). Two options — use the SAME dump the FOCUS DB
   was built from, already on disk:

   (a) taxonkit backend — no pip install, no download (recommended here):
       export TAXONKIT_DB=/path/to/taxonkit_dump   # or: source config.sh
       # profile2cami.py auto-selects taxonkit when this dump is present.

   (b) ete3 backend — build a pinned sqlite once (needs 'pip install ete3'):
       python3 setup_taxdb.py --url /path/to/taxdump.tar.gz \\
                              -o \$OUTDIR/taxdb/taxa.sqlite
       export PROFILE2CAMI_TAXDB=\$OUTDIR/taxdb/taxa.sqlite

   (dump 2026-07-10, sha256 c1b91199…; verified 10/10 updated phyla + 100% of the
   FOCUS output species resolve). On a machine without the dump, pass setup_taxdb.py
   a dated archive URL: .../taxdump_archive/taxdmp_<YYYY-MM-DD>.zip

1) Convert each profile to CAMI/BIOBOXES format with profile2cami.py (auto-detects
   FOCUS / Bracken / Kraken2 report / MetaPhlAn), then OPAL vs the gold standard:

     # RaPDTool full (assembly-based) and screen (read-based, matched input vs the classifiers):
     python3 scripts/profile2cami.py \$OUTDIR/rapdtool.rep1/profilesfmbm/*/output_All_levels.csv \\
             -o \$OUTDIR/rapdtool_full.profile -s rapdtool_full
     python3 scripts/profile2cami.py \$OUTDIR/rapdtool_screen.rep1/profilesfmbm/*/output_All_levels.csv \\
             -o \$OUTDIR/rapdtool_screen.profile -s rapdtool_screen
     python3 scripts/profile2cami.py \$OUTDIR/kraken2_full.rep1.bracken  -o \$OUTDIR/kraken_full.profile  -s kraken2_full
     python3 scripts/profile2cami.py \$OUTDIR/kraken2_cap16.rep1.bracken -o \$OUTDIR/kraken_cap16.profile -s kraken2_cap16
     python3 scripts/profile2cami.py \$OUTDIR/kraken2_cap8.rep1.bracken  -o \$OUTDIR/kraken_cap8.profile  -s kraken2_cap8
     python3 scripts/profile2cami.py \$OUTDIR/metaphlan.rep1.profile     -o \$OUTDIR/mpa.profile          -s metaphlan

   For CAMI datasets the gold standard ships in CAMI format already. For your OWN
   mock, build it from a 2-column table ('<taxid|name><TAB>abundance'):
     python3 scripts/profile2cami.py mock_composition.tsv -f truth -o \$OUTDIR/gold_standard.profile -s gold

     scripts/... opal.py -g \$OUTDIR/gold_standard.profile \\
             \$OUTDIR/rapdtool_screen.profile \$OUTDIR/rapdtool_full.profile \\
             \$OUTDIR/kraken_full.profile \$OUTDIR/kraken_cap16.profile \$OUTDIR/kraken_cap8.profile \\
             \$OUTDIR/mpa.profile \\
             -o \$OUTDIR/opal   # rapdtool_screen = matched-input (reads) head-to-head;
                                # rapdtool_full = assembly-based; Kraken curve full/16/8
   -> genus/species recall, precision, L1 norm, Bray-Curtis, weighted UniFrac.
   (profile2cami.py resolves taxonomy via taxonkit or ete3 per step 0; no other setup.)

2) Bin accuracy (RaPDTool only; others produce no bins) vs CAMI binning gold standard:
     amber.py -g gold_standard_binning.tsv rapdtool_binning.tsv -o $OUTDIR/amber

3) MAG quality already in RaPDTool output (miComplete); cross-check with CheckM if desired.

Fill the manuscript results tables (Tables 1–5) with the medians above + the OPAL/AMBER outputs.
Do not fabricate any value.
EOF
