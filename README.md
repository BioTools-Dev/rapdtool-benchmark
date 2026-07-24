# RaPDTool benchmark kit

Reproducible benchmark of RaPDTool against read/marker-based classifiers
(Kraken2 (+Bracken) and MetaPhlAn4; mOTUs is wired in as an optional comparator
but was not part of the reported runs) for the manuscript. It measures
**footprint** (DB size, peak RAM, wall time), **profiling accuracy** (OPAL vs a
gold standard), and **genome-recovery / type-material** capabilities that read
classifiers do not provide.

> Analysis code, protocol, community-defining inputs, and summary result files for the
> RaPDTool benchmark. Configure paths in `config.sh` (see Configuration), then follow
> Steps 0 → 6.

**Reading order.** For *why* the benchmark is designed this way — the flow diagram,
what each test measures, and an honest account of what is strong and what is missing —
read **`benchmark_rationale.md`** first. This file is the operational how-to: follow
Steps 0 → 6 below and a reviewer can reproduce every number in the manuscript.

## Directory layout

```
rapdtool-benchmark/
  *.md          docs (this README + rationale, methods, mock_design, results tables)
  scripts/      all analysis scripts (.py, .sh)
  data/         inputs: genome lists, abundance vectors, the census TSV, flow diagram
  results/      bench_<dataset>/ tool outputs, profiles, OPAL results
  figures/      generated figures (SVG + PNG 300 dpi + PDF) — manuscript deliverables
```

**Run every command from this base directory** (e.g. `scripts/mash_detection.py …`,
`python3 scripts/plot_census.py`). Scripts resolve `data/`, `results/` and `figures/`
relative to it; the shell scripts `cd` here themselves.

## Configuration

External paths (reference databases, tool locations) are read from environment
variables, not hardcoded. **Once, before running anything that touches a database:**

```bash
cp config.sh.example config.sh     # then edit config.sh with your paths
source config.sh
```

`config.sh` is git-ignored; only the template is committed.

**Two levels of reproduction:**

1. **Verify the figures and reported numbers from the shipped summaries.** The small
   result files under `results/` (CAMI profiles, OPAL `results.tsv`, `summary.csv`, mash
   confidence tables, miComplete tables) and the inputs under `data/` are included, so
   the plotting and detection scripts regenerate every figure and number directly. This
   needs only `$TAXONKIT_DB` (a pinned NCBI taxdump) for the taxid rollups:
   ```bash
   source config.sh
   # the five figure scripts — each regenerates from the shipped result files alone:
   python3 scripts/plot_census.py -i data/census_full.tsv -l data/mock_genomes.list \
       -o figures/census                          # Fig 2  -> figures/census.*
   python3 scripts/plot_opal_depth.py             # Fig 3  -> figures/opal_depth.*
   python3 scripts/plot_micomplete.py             # Fig 4  -> figures/micomplete.*
   python3 scripts/plot_mirror_distance.py        # Fig 5  -> figures/mirror_distance.*
   python3 scripts/plot_f1_threshold.py           # threshold analysis -> figures/f1_threshold.*
   # detection numbers (recall / precision / F1) from the mash confidence table:
   python3 scripts/mash_detection.py --bench results/bench_ln_30M \
       --gold results/bench_ln_30M/gold_standard.profile --split data/mock_genomes.list
   ```
   (`scripts/threshold_sweep.py --truth <gold_standard.profile>` writes the
   precision–recall sweep, `figures/threshold_sweep.*`; the flow diagram is rendered
   from `data/benchmark_flow.dot` with `dot`.)
2. **Full re-run from scratch.** Regenerate the simulated reads
   (`scripts/make_mock.sh … --seed 42`, from the shipped genome lists and abundance
   vectors), run the tools (`scripts/run_benchmark.sh`), and re-derive everything. This
   needs the reference databases named in `config.sh` (Kraken2 ≈104 GB, MetaPhlAn ≈60 GB,
   RaPDTool ≈0.5 GB) and the ZymoBIOMICS reads (ENA **PRJEB29504**). Large regenerable
   data (reads, assemblies, raw Kraken/Bracken output) are **not** shipped — they are
   git-ignored and reproduced from the recipe.

Large read/assembly data live outside the kit at `$MOCK_ROOT/mock_*` and `$ZYMO_DIR`.

## Files

| File | Purpose |
|------|---------|
| **`benchmark_rationale.md`** | **why the benchmark is shaped this way** — flow diagram, what each test measures, honest strengths/weaknesses. Read this first. |
| `check_representation.py` | census a genome set against Kraken2 (full/16/8) and MetaPhlAn4: is each species in each database? Produces both the population census and the mock genome selection |
| `data/census_full.tsv` | the census of all 30,209 type-material genomes (Phase 0 result) |
| `plot_census.py` / `figures/census.*` | census figure (PNG 300 dpi + SVG + PDF) |
| `data/mock_genomes.list` / `mock_design.md` | the 20 selected genomes (10 reference + 10 conflictive) and the sampling rationale |
| `make_mirror_set.py` | select the **mirror** set: species Kraken2 has and RaPDTool lacks, to characterise out-of-domain failure mode |
| `data/mirror_genomes.tsv` / `.list` / `data/mirror_pool.tsv` | the 10 selected mirror genomes, and the full 9,640-species pool they came from |
| `make_abundance.py` | build a **fixed** per-contig abundance vector — required, see the InSilicoSeq caveat below |
| `make_mirror_distance.py` / `data/mirror_distance.tsv` / `data/mirror_dist_genomes.list` | distance-stratified mirror: sample across genomic distance, measure Mash distance to RaPDTool's DB |
| `plot_mirror_distance.py` / `figures/mirror_distance.*` | rank-resolution-vs-distance figure |
| `data/mock_abundance.txt` / `data/mock_abundance_equalcov.txt` / `data/mirror_dist_abundance.txt` | the vectors actually used (uneven 20×; equal-coverage; mirror equal-coverage) |
| `make_mock.sh`       | build a mock community: simulate reads (iss) → assemble (MEGAHIT) → write the `mock_composition.tsv` gold standard |
| `run_benchmark.sh`   | run each tool on the same input/hardware; record RAM, time, DB sizes; report medians |
| `setup_taxdb.py`     | pin the NCBI taxonomy dump (build an ete3 SQLite) for reproducible name/taxid resolution |
| `profile2cami.py`    | convert FOCUS / Bracken / Kraken2 / MetaPhlAn profiles **and** a gold standard to CAMI/BIOBOXES for OPAL (taxonkit or ete3 backend) |
| `mash_detection.py`  | RaPDTool species **detection** from the mash confidence table (recall/precision/F1 + ref/conflictive split); detection = mash, abundance = FOCUS |
| `run_opal.sh` / `run_opal_zymo.sh` | convert every profile (matching SampleID) and run OPAL, per dataset; the Zymo variant handles the sequence- vs cell-abundance split |
| `run_metawrap.sh`    | MAG-recovery comparison vs MetaWRAP on the **same assembly**: binning + bin_refinement, then re-scores both bin sets with the same evaluator (miComplete/Bact105) → `results/metawrap/` (Step 5d, rationale §4d) |
| `threshold_sweep.py` / `figures/threshold_sweep.*` | precision–recall vs abundance cutoff (reads the CAMI profiles; picks the operating point) |
| `plot_f1_threshold.py` / `figures/f1_threshold.*` | F1 vs a uniform abundance threshold — the MetaPhlAn crossover |
| `plot_opal_depth.py` / `figures/opal_depth.*` | profiling accuracy (Bray–Curtis, F1) vs sequencing depth |
| `config.sh.example` | template of external paths (databases, tools); copy to `config.sh` and edit |
| `README.md`          | this file — the operational step-by-step |

## The datasets this benchmark uses

| Dataset | What it is | Answers |
|---|---|---|
| *(census, no reads)* | all 30,209 type-material genomes vs each competitor DB | what each database **contains** |
| `mock_ln_3M` | 20 genomes, uneven (20× range), 3 M reads | depth response, low |
| `mock_ln_10M` | **same vector**, 10 M reads | depth response, mid |
| `mock_ln_30M` | **same vector**, 30 M reads | depth response, high |
| `mock_equalcov_30M` | 20 genomes, equal coverage (38.8×), 30 M reads | capability with coverage removed as a limit |
| `mock_mirrordist_15M` | 14 genomes spanning **genomic distance** to RaPDTool's DB (100 % → ~70 % id), equal coverage (27.7×), 15 M reads | out-of-domain **rank resolution vs distance**: does it degrade gracefully? |
| `ZymoBIOMICS_db` | **real** Illumina data, community defined by a third party | does it work on data we did not simulate? |

The three `ln` sets share one fixed abundance vector so that **depth is the only
variable** between them. Locations: simulated mocks under `$MOCK_ROOT/mock_*`, real data
under `$ZYMO_DIR`, mirror genomes downloaded to `$MIRROR_FNA_DIR`.

Together they answer four different questions, and none substitutes for another:
coverage of the databases (census), performance under realistic unevenness (depth
series), capability when coverage is not limiting (equal-coverage control), safety
outside the intended domain (mirror), and external validity (Zymo). See
`benchmark_rationale.md` for why each exists.

## Prerequisites

- The tools you compare, reachable by **absolute path** (see the warning in Step 1 —
  do not rely on `PATH`): `rapdtool`, `kraken2`, `bracken`, `metaphlan` (and/or
  `motus`); GNU `/usr/bin/time`.
- Accuracy step: OPAL (`cami-opal`) and AMBER (`cami-amber`), installed in a
  **separate Python 3.11 conda env** (see A below — they pin numpy/pandas versions
  that fail to build on Python ≥3.13).
- A taxonomy backend for `profile2cami.py` (auto-selected):
  - **taxonkit** + a local NCBI dump — no pip install, no download. Already set up
    here: `export TAXONKIT_DB=$TAXONKIT_DB` (the dump the FOCUS DB was
    built from, 2026-07-10, sha256 `c1b91199…`; 10/10 updated phyla + 100% of FOCUS
    output species resolve). **Recommended.**
  - or **ete3** (`pip install ete3`) with a pinned sqlite from `setup_taxdb.py`.

---

## Getting started (first time on this machine)

**A. Install the tools (once).** RaPDTool you already have; the rest go in a
dedicated conda env:

```bash
# Tools env (any recent Python is fine here)
conda create -n rapdtool_bench -c conda-forge -c bioconda \
    kraken2 bracken metaphlan insilicoseq megahit taxonkit -y
conda activate rapdtool_bench
# RaPDTool runs via its launcher wrapper; if `rapdtool` is not on PATH in this env,
# pass its absolute path as RAPDTOOL=... in Step 1 (e.g. .../RaPDTool/scripts/rapdtool)
```

**OPAL/AMBER go in a SEPARATE Python 3.11 env.** `cami-opal` pins exact
`numpy==2.0.1` / `pandas==2.2.2`, which have no wheels for Python ≥3.13 and then
fail to compile from source with new GCC. The OPAL wheel itself is `…-py311-…`, so
use 3.11 and pip installs prebuilt wheels (no compilation):

```bash
conda create -n opal -c conda-forge python=3.11 -y
conda activate opal
pip install cami-opal cami-amber
```

`profile2cami.py` needs neither numpy nor pandas (stdlib + taxonkit only), so run it
in `rapdtool_bench`; switch to the `opal` env only for the `opal.py` / `amber.py`
calls. (Add `ete3` to whichever env runs `profile2cami.py` only if you are not
using the taxonkit backend.)

**B. Download the databases.**
- Kraken2 **Standard** (full, ~104 GB on disk), **Standard-16** (~15 GB) and
  **Standard-8** (~7.5 GB) from the prebuilt indexes at
  https://benlangmead.github.io/aws-indexes/k2 (the k2 tarballs already include the
  Bracken `*.kmer_distrib` files). Full + the two capped sizes give the Kraken
  accuracy-vs-DB-size curve vs RaPDTool's 0.5 GB. Needs ~130 GB free disk; loading the
  full DB needs ~104 GB RAM. **This run used the 2026-06-26 build (place it at `$KRAKEN_STD_DB` etc.):**
  ```bash
  mkdir -p /path/to/kraken2_db/{standard,standard16,standard8} && cd /path/to/kraken2_db
  wget https://genome-idx.s3.amazonaws.com/kraken/k2_standard_20260626.tar.gz
  wget https://genome-idx.s3.amazonaws.com/kraken/k2_standard_16_GB_20260626.tar.gz
  wget https://genome-idx.s3.amazonaws.com/kraken/k2_standard_08_GB_20260626.tar.gz
  tar -xzf k2_standard_20260626.tar.gz    -C standard
  tar -xzf k2_standard_16_GB_20260626.tar.gz -C standard16
  tar -xzf k2_standard_08_GB_20260626.tar.gz -C standard8
  ```
- MetaPhlAn DB (**~60 GB** for the current `vJan26_CHOCOPhlAnSGB` SGB build — ~21 GB
  marker `.fna` + ~39 GB prebuilt Bowtie2 index `.bt2l`), once:
  `metaphlan --install --db_dir ~/metaphlan_db`
  (MetaPhlAn 4.1+ renamed `--bowtie2db` → `--db_dir`; pick any path with space, e.g.
  `$MPA_DB`. Then set `MPA_DB=<that dir>` in `run_benchmark.sh` CONFIG.)
- RaPDTool DBs: already cached (`rapdtool --where`).

**C. Build the datasets.** The 20 genomes are already selected in `data/mock_genomes.list`;
`make_mock.sh` simulates reads, assembles them, and writes the gold standard.

> ### ⚠ Read this before generating any mock
>
> **`iss --abundance <dist>` assigns abundance per FASTA *record* (contig), not per
> genome.** With draft multi-contig genomes the realised per-genome abundance therefore
> tracks **contig count**, largely regardless of which distribution you request.
> Measured on this genome set with `--abundance uniform`: a 1,491-contig genome got
> **61.3 %** of the reads and a single-contig genome got **0.04 %** — a 1,000-fold
> coverage range from a request for *uniform*.
>
> **Always pass a fixed vector with `-b`**, built by `make_abundance.py`, which assigns
> per-genome fractions explicitly and splits each across that genome's contigs
> proportionally to length. Never rely on `-a`.

```bash
# run from the repository base directory
conda activate rapdtool_bench

# C.1 — build the two abundance vectors (fast, no reads involved)
scripts/make_abundance.py -l data/mock_genomes.list -o data/mock_abundance.txt --depths 3 10 30
scripts/make_abundance.py -l data/mock_genomes.list -o data/mock_abundance_equalcov.txt \
                    --equal-coverage --depths 30
```

Both commands print the per-genome coverage they will produce at each depth — check
that table before spending an hour on simulation.

```bash
# C.2 — the depth series: ONE vector, three depths (~30 min + assemblies)
for N in 3000000 10000000 30000000; do
  scripts/make_mock.sh -o $MOCK_ROOT/mock_ln_$((N/1000000))M -t 16 -n $N \
    -b data/mock_abundance.txt --seed 42 -l data/mock_genomes.list
done

# C.3 — the equal-coverage positive control (38.8x for all 20 genomes)
scripts/make_mock.sh -o $MOCK_ROOT/mock_equalcov_30M -t 16 -n 30000000 \
  -b data/mock_abundance_equalcov.txt --seed 42 -l data/mock_genomes.list
```

Each run produces, under its `-o` directory:
`reads_R{1,2}.fastq`, `asm/final.contigs.fasta`, `mock_composition.tsv` (the gold
standard), and `make_mock.log`.

Notes that will save you a debugging session:

- **Use `final.contigs.fasta`, not MEGAHIT's `final.contigs.fa`.** The FOCUS version
  bundled in the RaPDTool container accepts only `.fasta`/`.fna`/`.fastq`; given `.fa`
  it aborts with a misleading `NameError: name 'sys' is not defined` that hides the
  real message. `make_mock.sh` writes the `.fasta` copy for you.
- **`iss` does not write `reads_abundance.txt` when `-b` is used** (the composition was
  yours to begin with). `make_mock.sh` copies your vector into place so the gold
  standard still builds.
- `-n` is **total reads**, not pairs: `-n 30000000` gives 15 M pairs.
- `--seed` makes the read simulation reproducible for a given `-n` and `-b`.

Database paths are passed as environment variables at run time (Step 1); the `CONFIG`
block in `run_benchmark.sh` holds only defaults, so you never edit the script.

**How the 20-genome list was produced (and how to regenerate/customise it).**
`data/mock_genomes.list` is not hand-waved — it was derived with `check_representation.py`,
which reads each competitor DB's *own* taxonomy (no downloads) and reports, per
genome species, presence in Kraken2 (full/16/8) and MetaPhlAn4. To reproduce it:

1. **Confirm the DB paths** the script reads (hard-coded constants near the top of
   `check_representation.py`): `KRAKEN` → `$KRAKEN_DB_FULL{,16,8}/inspect.txt`,
   `MPA` → the MetaPhlAn `.pkl`, and `FOCUS` → `$FOCUS_DB`
   (`acc_taxid_strain.tsv`, `taxid_lineage.tsv`, `db/`). Edit them if yours differ.

2. **Sample and classify** a pool of type-material genomes (deterministic via `--seed`):
   ```bash
   conda activate rapdtool_bench
   python3 scripts/check_representation.py --sample 1500 --seed 42 -o data/bench_out_repcheck.tsv
   ```
   It prints per-category tallies (`reference`, `conflictive-both`, `conflictive-mpa`,
   `conflictive-kraken`), then **suggests 10 reference + 10 conflictive** genomes with a
   ready-to-paste GCA path list; `data/bench_out_repcheck.tsv` holds the full per-genome table
   (columns: gca, taxid, phylum, genus, species, kraken_full/16/8, metaphlan, verdict).

3. **Pick your 20.** Either take the suggested list as-is, or curate — this kit used
   **5 clinical + 5 phylum-diversifiers** (all `reference`, i.e. present in every DB) and
   **10 `conflictive-both`** (absent from Kraken *and* MetaPhlAn), one per phylum. To add
   a specific canonical species not caught by the random sample, check it directly:
   ```bash
   python3 scripts/check_representation.py GCA_001281725.1 GCA_002902205.1   # E. coli, S. aureus
   #   -> confirms verdict=reference (1/1/1, Y) before you add it
   ```
   Options: positional GCAs / paths, `-l FILE`, `--n-ref`, `--n-conf`, `--sample`, `--seed`.

4. **Write `data/mock_genomes.list`** — one genome path per line (inline `#` comments and blank
   lines are ignored by both `check_representation.py -l` and `make_mock.sh -l`). The
   verdicts and full rationale for the shipped list are in `mock_design.md`.

> Tip: dry-run the light path first — the census (Step 0b) and the figure/number
> verification need only `$TAXONKIT_DB`, no large database download. Confirm those work
> before committing disk and time to the ~104 GB Kraken2 index for a full re-run.

---

## The complete workflow

### Step 0 — Pin the taxonomy (once)

Harmonises every profile to a single taxonomy so OPAL matches taxa consistently
and the updated phylum names from FOCUS resolve, using the **same dump the FOCUS DB
was built from** (already on disk — no download). Pick one backend:

**(a) taxonkit — recommended, no install/download:**
```bash
export TAXONKIT_DB=$TAXONKIT_DB   # profile2cami.py auto-selects taxonkit
```

**(b) ete3 — build a pinned sqlite once (`pip install ete3`):**
```bash
# run from the repository base directory
python3 scripts/setup_taxdb.py --url $TAXONKIT_DB/../taxdump.tar.gz \
                       -o data/taxdb/taxa.sqlite
export PROFILE2CAMI_TAXDB=$PWD/data/taxdb/taxa.sqlite
```
Writes `taxdb.version` (source + SHA-256 + phylum-resolution check) — keep it with
the analysis. On another machine, pass a dated archive URL instead:
`https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump_archive/taxdmp_<YYYY-MM-DD>.zip`.

Force a backend explicitly with `profile2cami.py --backend {taxonkit,ete3}`.

### Step 0b — The database census (no reads; run this first)

The strongest result in the study and the cheapest to reproduce. It measures what each
database *contains*, independently of any simulated community:

```bash
conda activate rapdtool_bench
scripts/check_representation.py --sample 30209 --seed 42 -o data/census_full.tsv   # ~3 min
$BENCH_ENV_BIN/python scripts/plot_census.py \
    -i data/census_full.tsv -l data/mock_genomes.list -o figures/census
```

Produces `data/census_full.tsv` (one row per genome: gca, taxid, phylum, genus, species,
kraken_full/16/8, metaphlan, verdict) and the figure in PNG + SVG + PDF.

Result on this system: of 30,209 type-material genomes, **55.8 % are absent from the
Kraken2 standard database** (103.7 GB), 17.3 % from MetaPhlAn4, and **11.3 %
(3,423) from both**. Capping Kraken2 to 8.1 GB changes the absence rate by <0.01 pp.

> This step also defines the population the mock genomes are sampled from, which is
> what makes the selection defensible rather than cherry-picked. See
> `benchmark_rationale.md` §2 and `mock_design.md`.

### Step 1 — Run the tools (footprint metrics)

Pass the inputs/DBs as environment variables (they override the `CONFIG` defaults, so
you don't edit the script). First a **`REPEATS=1` smoke test** to validate the
plumbing and that Kraken-full fits in RAM, then re-run with `REPEATS=3` for the paper
medians — same command, only `REPEATS` changes:

> ### ⚠ Pass ABSOLUTE tool paths — do not rely on `PATH`
> This cost a full overnight run. With the conda env inactive — or on any host whose
> `PATH` contains relative entries such as `.`, `./bin`, `./scripts` (as this one's
> `.zshrc` does), making resolution depend on the current directory — every tool exits
> 126/127 in 0.0 s, and the whole matrix "completes" in minutes with no data.
>
> `BENCH_ENV_BIN` points at the env's `bin/` and sets `kraken2`, `bracken` and
> `metaphlan` at once; `RAPDTOOL` points at the wrapper. **With these, `conda activate`
> is not needed.** A preflight check now aborts in seconds if anything is missing —
> you should see `[preflight] all enabled tools and inputs present` before the first
> tool starts.

```bash
M=$MOCK_ROOT/mock_ln_30M                   # one dataset; full matrix below

BENCH_ENV_BIN=$BENCH_ENV_BIN \
RAPDTOOL=$RAPDTOOL \
THREADS=16 REPEATS=1 OUTDIR=$PWD/results/bench_ln_30M \
ASSEMBLY=$M/asm/final.contigs.fasta \
READS_R1=$M/reads_R1.fastq \
READS_R2=$M/reads_R2.fastq \
KRAKEN_DB_FULL=$KRAKEN_DB_FULL \
KRAKEN_DB_CAP16=$KRAKEN_DB_CAP16 \
KRAKEN_DB_CAP8=$KRAKEN_DB_CAP8 \
MPA_DB=$MPA_DB \
scripts/run_benchmark.sh
```

Individual launchers can also be set one by one (`KRAKEN2=`, `BRACKEN=`,
`METAPHLAN=`) if the tools live in different places.

**Give each dataset its own `OUTDIR`** — `run_benchmark.sh` *appends* to
`summary.csv`, so a shared output directory silently mixes datasets into one file.
To run the whole matrix (~6 h; MetaPhlAn scales linearly with read count):

```bash
# run from the repository base directory

for D in ln_3M ln_10M ln_30M equalcov_30M; do
  M=$MOCK_ROOT/mock_$D
  BENCH_ENV_BIN=$BENCH_ENV_BIN \
  RAPDTOOL=$RAPDTOOL \
  THREADS=16 REPEATS=1 OUTDIR=$PWD/results/bench_$D \
  ASSEMBLY=$M/asm/final.contigs.fasta \
  READS_R1=$M/reads_R1.fastq READS_R2=$M/reads_R2.fastq \
  KRAKEN_DB_FULL=$KRAKEN_DB_FULL \
  KRAKEN_DB_CAP16=$KRAKEN_DB_CAP16 \
  KRAKEN_DB_CAP8=$KRAKEN_DB_CAP8 \
  MPA_DB=$MPA_DB \
  scripts/run_benchmark.sh
done
```

**Sanity check while it runs:** the first `summary.csv` row must show a wall time
well above 0.0 s. A row at 0.0 s with `exit_code` 126 or 127 means the tool never
started — stop immediately rather than letting the night pass.

**Re-running after a partial failure:** prefer deleting the affected `bench_*`
directories and starting that dataset again, over editing `summary.csv` by hand. The
outputs of one directory should always come from one command with one configuration —
hand-pruned rows leave orphan `.stderr`/`.time.txt` files that are impossible to
interpret later. The mock datasets under `$MOCK_ROOT/mock_*` are expensive and validated;
never delete those to fix a benchmark problem.

**`REPEATS=1` is correct for accuracy**; all compared tools are deterministic, so
replicating accuracy measures nothing. Use `REPEATS=3` on **one** dataset to get the
resource medians for the paper, and state which dataset (and therefore which depth)
those medians came from.

Individual tools can be switched off — every toggle is environment-overridable:
`RUN_RAPDTOOL_FULL`, `RUN_RAPDTOOL_SCREEN`, `RUN_KRAKEN_FULL`, `RUN_KRAKEN_CAP16`,
`RUN_KRAKEN_CAP8`, `RUN_METAPHLAN`, `RUN_MOTUS`. e.g. re-run only RaPDTool:

```bash
RUN_KRAKEN_FULL=0 RUN_KRAKEN_CAP16=0 RUN_KRAKEN_CAP8=0 RUN_METAPHLAN=0 ... scripts/run_benchmark.sh
```

> **RaPDTool launcher.** RaPDTool ships as an Apptainer wrapper installed in its own
> conda env (`rapdtool`). The benchmark itself runs from `rapdtool_bench` (which has
> `kraken2`/`metaphlan`/… and also `/usr/bin/apptainer`). Passing the **absolute wrapper
> path** as `RAPDTOOL=$RAPDTOOL` (above) runs it correctly
> from `rapdtool_bench` and lets `/usr/bin/time` measure RaPDTool's real peak RSS (the
> wrapper `exec`s apptainer). If that ever fails, fall back to running it in its own env
> with `RAPDTOOL_ENV=rapdtool` (drop the `RAPDTOOL=…` line) — but note `conda run` adds a
> persistent python parent, so peak RSS for RaPDTool may be under-reported.

Check that `summary.csv` in the OUTDIR has `exit_code=0` on every row before trusting
any result. Produces, under the OUTDIR:
- `summary.csv` — per-replicate wall-clock + peak RSS for each tool
- `db_sizes.csv` — on-disk database sizes
- median wall-time / RAM printed to the terminal
- each tool's outputs (`rapdtool.repN/…`, `*.bracken`, `*.report`, `*.profile`)

Also record for RaPDTool the DB sizes from `rapdtool --where` (mash + FOCUS) into
`db_sizes.csv`.

### Step 2 — Convert every profile to CAMI

`profile2cami.py` auto-detects the format; it resolves each species leaf to an NCBI
taxid and derives the **full** standard-rank lineage from the pinned taxonomy, so all
tools land on the same tree.

> ### ⚠ Every profile must carry the SAME `SampleID`
> Pass **`-s mock` to every conversion, including the gold standard.** OPAL matches a
> profile to the gold standard by sample ID; any profile whose ID differs is *silently
> skipped* with a warning, and if all of them differ OPAL reports
> `No profile could be evaluated`. Do not name the samples after the tools — the tool
> name is carried by `-l` at the OPAL step, not by the sample ID.

```bash
export TAXONKIT_DB=$TAXONKIT_DB    # Step 0 backend (see above)
B=results/bench_ln_30M                                  # the OUTDIR from Step 1
mkdir -p $B/profiles

# RaPDTool full (assembly) and screen (reads, matched input vs the classifiers):
python3 scripts/profile2cami.py $B/rapdtool.rep1/profilesfmbm/*/output_All_levels.csv \
        -f focus -s mock -o $B/profiles/rapdtool_full.profile
python3 scripts/profile2cami.py $B/rapdtool_screen.rep1/profilesfmbm/*/output_All_levels.csv \
        -f focus -s mock -o $B/profiles/rapdtool_screen.profile
for k in full cap16 cap8; do
  python3 scripts/profile2cami.py $B/kraken2_${k}.rep1.bracken -f bracken -s mock \
          -o $B/profiles/kraken2_${k}.profile
done
python3 scripts/profile2cami.py $B/metaphlan.rep1.profile -f metaphlan -s mock \
        -o $B/profiles/metaphlan.profile
```

Check the `abundance mapped=…%` line on stderr for each; investigate any
`unresolved:` names if coverage < 90 %. All six should report 100 %.

> ### ⚠ Use the UNFILTERED FOCUS output, not the confidence table
> RaPDTool applies a **1 % abundance cutoff** in `rapdtool_confidence.tbl` / `.txt`.
> Feed OPAL the raw profile — `profilesfmbm/*/output_All_levels.csv` (207 species on
> `mock_ln_30M`) — **not** the confidence table (16 species).
>
> Filtering RaPDTool at 1 % while Bracken stays unfiltered would manufacture the
> precision result: 8,129 species vs 16 compares a filtered output against an
> unfiltered one, which measures nothing.
>
> Report instead a **threshold sweep applied to every tool alike** — none / 0.1 % / 1 %.
> The 1 % row is then RaPDTool's default operating point evaluated fairly, with the
> same filter applied to the competitors.
>
> **Note for the 1 % row:** two of the twenty genomes are below 1 % by design (0.89 %
> and 0.76 %), so maximum recall at that threshold is **18/20 for every tool**.

### Step 3 — Build the gold standard

- **CAMI datasets**: the gold-standard profile ships in CAMI format already — use it directly.
- **Your own mock**: from the 2-column table `make_mock.sh` wrote
  (`<taxid|name><TAB>abundance`). **Same `-s mock` as every other profile:**

```bash
python3 scripts/profile2cami.py $MOCK_ROOT/mock_ln_30M/mock_composition.tsv -f truth \
        -s mock -o results/bench_ln_30M/gold_standard.profile
```

Each dataset has its **own** gold standard — the equal-coverage control and the depth
series are different communities. Never evaluate a profile against another dataset's
gold standard.

`mock_composition.tsv` example (tab- or comma-separated; abundance scale is
renormalised per rank):

```
# taxid or name <TAB> relative_abundance
470	60
Staphylococcus aureus	20
Acinetobacter_johnsonii	20
```

### Step 4 — Profiling accuracy with OPAL

**Shortcut — Steps 2 + 4 together for every dataset:** `run_opal.sh` converts all six
profiles + the gold standard (matching SampleID) and runs OPAL, aborting a dataset if
any conversion yields an empty profile:

```bash
conda activate rapdtool_bench            # profile2cami needs taxonkit
export TAXONKIT_DB=$TAXONKIT_DB
scripts/run_opal.sh ln_3M ln_10M ln_30M equalcov_30M
```

Or run OPAL by hand on an already-converted dataset:

```bash
B=results/bench_ln_30M
OPALPY=$OPAL_PY   # calling the env's python directly
                                                     # avoids conda-activate issues
$OPALPY $OPAL_BIN \
  -g $B/gold_standard.profile -o $B/opal \
  $B/profiles/rapdtool_screen.profile $B/profiles/rapdtool_full.profile \
  $B/profiles/kraken2_full.profile $B/profiles/kraken2_cap16.profile \
  $B/profiles/kraken2_cap8.profile $B/profiles/metaphlan.profile \
  -l "RaPDTool_screen,RaPDTool_full,Kraken2_full,Kraken2_16GB,Kraken2_8GB,MetaPhlAn4"
```

Gives per-rank **recall (Completeness), precision (Purity), F1, L1 norm, Bray–Curtis,
weighted UniFrac** in `$B/opal/results.tsv` plus an HTML report. Report and interpret at
**genus/species** (the ranks RaPDTool targets). What each metric means, and how to avoid
misreading them, is in `benchmark_rationale.md` §3.

> OPAL gives per-rank recall/precision/F1, L1 and Bray–Curtis. **For RaPDTool, use OPAL
> for ABUNDANCE only (L1, Bray–Curtis).** Its OPAL recall/precision reflect the FOCUS
> profile's false-positive tail and understate detection — detection comes from mash
> (Step 4b). The competitors' OPAL detection metrics are valid (their profiles are their
> detection).

### Step 4b — Species detection from the mash confidence table

> **RaPDTool has two species outputs.** DETECTION is the mash-screen confidence table
> (`rapdtool_confidence.tbl`), its confident species calls — the analogue of MetaPhlAn's
> marker-filtered list. ABUNDANCE is the FOCUS profile (Step 4a), which the tool itself
> flags "cautious at species level". Scoring detection from FOCUS understates RaPDTool
> badly (FOCUS ~204 species / ~184 false positives; mash: the true species, 0 FP). So
> **detection = mash, abundance = FOCUS.**

```bash
export TAXONKIT_DB=$TAXONKIT_DB
# recall / precision / F1 + reference-vs-conflictive split, from the mash table:
scripts/mash_detection.py --bench results/bench_ln_30M --gold results/bench_ln_30M/gold_standard.profile \
                    --split data/mock_genomes.list
# real data (credit a reclassified member as its gold taxon, e.g. B. spizizenii -> subtilis):
scripts/mash_detection.py --bench results/bench_zymo_even --gold results/bench_zymo_even/gold_dna.profile \
                    --extra-true 96241:1423
```

On the mocks this gives recall 1.0, precision 1.0, F1 1.0, **0 false positives** (0.97 at
3 M), and the split 10/10 reference + 10/10 conflictive — the central experiment.
RaPDTool detects the conflictive species that the competitors (0/10, absent from their
databases) cannot. **Never present the conflictive result without the census (Step 0b).**

**Abundance-threshold analysis (context for the FOCUS profile, not detection):**

```bash
scripts/threshold_sweep.py  --truth $B/gold_standard.profile -o figures/threshold_sweep
scripts/plot_f1_threshold.py                                  # -> figures/f1_threshold.*
scripts/plot_opal_depth.py                                    # -> figures/opal_depth.*  (after all 4 datasets)
```

> The FOCUS abundance profile carries a low-abundance tail; under an abundance cutoff
> applied uniformly to every tool's output it overtakes MetaPhlAn on F1 (1.0 at 0.5 %).
> This is a caveat for interpreting FOCUS composition — **RaPDTool's detection F1 is
> already 1.0 from the mash table, no threshold.** The cutoff is a post-hoc output filter,
> distinct from each tool's internal detection; full rationale in `benchmark_rationale.md`
> §3 and Phase 5.

> **Which mode is the fair comparison?** `rapdtool_screen` consumes the **same reads** as
> Kraken/MetaPhlAn — apples-to-apples, and it is where the mash confidence table lives.
> `rapdtool_full` is assembly-based and reserved for genome recovery (Step 5); the read
> classifiers cannot consume the assembly.

### Step 5 — Genome recovery (RaPDTool only)

Read classifiers produce no bins; this axis is a capability, not a contest.

Bin completeness/redundancy come from the RaPDTool output (miComplete, in
`results/bench_<dataset>/rapdtool.rep1/workfmbm/outmicomplete/miCompleteOut_*.tab`);
cross-check with CheckM if desired. Bins are in `.../rapdtool.rep1/species_bins/`.
Summarise and plot completeness vs contamination across datasets:

```bash
scripts/plot_micomplete.py      # -> figures/micomplete.{svg,png,pdf} (reads results/bench_ln_30M + equalcov)
```

The two chimeric bins (*Fusobacterium massiliense*, *Corallococcus praedator*) are
correctly flagged by miComplete's redundancy — report them (Table 5), do not hide them.

> ⚠ **AMBER is not runnable on this mock.** It needs `gold_standard_binning.tsv`, a
> contig→genome truth table, and `make_mock.sh` produces only a *composition* gold
> standard (`mock_composition.tsv`). So this step measures **bin quality, not bin
> correctness** — do not claim binning accuracy from it. To close the gap: either use a
> dataset that ships a binning gold standard (CAMI II), or extend `make_mock.sh` to
> track each simulated read's source genome and derive the contig→genome truth from the
> assembly (~a day's work, keeps everything in-house).

### Step 5b — Out-of-domain behaviour: rank resolution vs distance (mirror experiment)

Measures whether the rank RaPDTool resolves an organism to **tracks its genomic
distance** from the database — species only when genuinely close, genus at moderate
distance, abstention when far. Graceful degradation is the safe, desirable behaviour;
a confident species call for a distant genome would be the failure. See
`benchmark_rationale.md` §4b for the rationale.

```bash
export TAXONKIT_DB=$TAXONKIT_DB

# 1. distance-stratified selection: sample across novelty tiers, download, and MEASURE
#    each genome's minimum Mash distance to RaPDTool's database (the x-axis).
scripts/make_mirror_distance.py --per-tier 6 --seed 42 -o data/mirror_distance.tsv
#    -> then hand-pick ~2-3 per distance band into data/mirror_dist_genomes.list, spanning
#       100 % (positive controls, species IS in DB) down to ~70 % identity. The shipped
#       list already does this for 14 genomes.

# 2. equal-coverage vector (non-detection must be distance, not coverage)
scripts/make_abundance.py -l data/mirror_dist_genomes.list -o data/mirror_dist_abundance.txt \
    --equal-coverage --depths 15

# 3. build the mock. NOTE --acc-map: these genomes are NOT in focus_build, so the
#    default accession->taxid table cannot resolve them and the gold standard would
#    come out empty. The acc-map is taxid-per-genome, built alongside the list.
scripts/make_mock.sh -o $MOCK_ROOT/mock_mirrordist_15M -t 16 -n 15000000 \
  -b data/mirror_dist_abundance.txt --seed 42 -l data/mirror_dist_genomes.list \
  --acc-map data/acc_taxid_mirrordist.tsv

# 4. run RaPDTool full+screen, with Kraken2 as positive control (all are in its DB)
M=$MOCK_ROOT/mock_mirrordist_15M
BENCH_ENV_BIN=$BENCH_ENV_BIN \
RAPDTOOL=$RAPDTOOL \
THREADS=16 REPEATS=1 OUTDIR=$PWD/results/bench_mirrordist \
RUN_KRAKEN_CAP16=0 RUN_KRAKEN_CAP8=0 RUN_METAPHLAN=0 \
ASSEMBLY=$M/asm/final.contigs.fasta \
READS_R1=$M/reads_R1.fastq READS_R2=$M/reads_R2.fastq \
KRAKEN_DB_FULL=$KRAKEN_DB_FULL \
scripts/run_benchmark.sh

# 5. crosswalk resolved rank against measured distance, and plot
scripts/plot_mirror_distance.py     # -> figures/mirror_distance.{svg,png,pdf}
```

**What to read is the resolved rank vs distance, not an accuracy score.** For each input
genome (with its measured Mash identity), the finest rank RaPDTool resolved comes from:

- `results/bench_mirrordist/rapdtool_screen.rep1/rapdtool_confidence.tbl` — the mash-screen
  **species** calls (only above the ~95 % threshold)
- `.../profilesfmbm/*/output_All_levels.csv` — the FOCUS profile; a genome's **genus**
  appearing here (but not in the confidence table) is a genus-level call
- `results/bench_mirrordist/rapdtool.rep1/species_bins/` — bin names carry species where
  resolved and genus-only where not

**Result on the shipped 14-genome set (all behaved correctly):**

| Mash identity to nearest DB genome | rank resolved |
|---|---|
| 100 % (in DB, positive control) | species |
| 97–99 % | species (nearest congener) |
| 92–95 % | **genus only** (mash stops calling species at 95 %) |
| 82–91 % | genus (FOCUS profile) |
| 70–76 % | **abstains — not reported** |

No genome below the 95 % threshold received a species call, and nothing below ~80 %
identity was reported at all — i.e. **graceful degradation, a safety property**, not
silent misassignment. State it as "below ~80 % identity RaPDTool abstains" rather than a
precise cutoff: Mash distance saturates near 70–75 % identity, so the abstained genomes'
exact distances are not meaningful.

### Step 5c — Real-data validation (ZymoBIOMICS)

Real Illumina data for a third-party-defined community; the only dataset here not
simulated from genomes we chose. Full provenance, composition table and the two
disclosure requirements are in `$ZYMO_DIR/README.md` — **read it before
building the gold standard**, particularly the sequence-abundance vs cell-abundance
distinction, which biases the comparison if ignored.

```bash
cd $ZYMO_DIR && md5sum -c md5sums.txt     # always verify first

# even community (D6300); repeat with ERR2935805 for the log community
Z=$ZYMO_DIR
BENCH_ENV_BIN=$BENCH_ENV_BIN \
RAPDTOOL=$RAPDTOOL \
THREADS=16 REPEATS=1 OUTDIR=$PWD/results/bench_zymo_even \
READS_R1=$Z/ERR2984773_1.fastq.gz READS_R2=$Z/ERR2984773_2.fastq.gz \
KRAKEN_DB_FULL=$KRAKEN_DB_FULL \
KRAKEN_DB_CAP16=$KRAKEN_DB_CAP16 \
KRAKEN_DB_CAP8=$KRAKEN_DB_CAP8 \
MPA_DB=$MPA_DB \
RUN_RAPDTOOL_FULL=0 \
scripts/run_benchmark.sh
```

> RaPDTool `full` is disabled above because no assembly exists for these reads yet.
> To include it, assemble first (`megahit -1 … -2 … -o asm`) and pass
> `ASSEMBLY=asm/final.contigs.fasta` — remembering that FOCUS needs `.fasta`, not
> MEGAHIT's `.fa`.

Gold standards are already built: `zymo_composition_dna.tsv` (for Kraken2/Bracken) and
`zymo_composition_genomecopy.tsv` (for MetaPhlAn). Convert whichever applies with
`profile2cami.py -f truth -s mock`.

### Step 5d — MAG recovery vs a dedicated pipeline (MetaWRAP)

Shows a dedicated ensemble binner does not beat RaPDTool on recovery (rationale §4d).
MetaWRAP is given the **same assembly** RaPDTool consumed, so only binning/refinement
differs; the script then re-scores **both** bin sets with the same evaluator
(miComplete/Bact105, inside RaPDTool's SIF) so the numbers are comparable. Containerised —
nothing to compile.

```bash
# one-time setup (paths overridable via env; see config.sh: METAWRAP_SIF, CHECKM_DB, RAPDTOOL_SIF)
apptainer pull "$METAWRAP_SIF" docker://quay.io/biocontainers/metawrap-mg:1.3.0--hdfd78af_1
#   CheckM DB -> $CHECKM_DB (binning + bin_refinement, 1.4 GB):
#   wget .../CheckM_databases/checkm_data_2015_01_16.tar.gz && tar -xzf ... -C "$CHECKM_DB"

# one script per dataset: binning + bin_refinement + miComplete re-scoring, all containerised
METAWRAP_SIF=$METAWRAP_SIF CHECKM_DB=$CHECKM_DB RAPDTOOL_SIF=$RAPDTOOL_SIF \
  scripts/run_metawrap.sh ln_30M $MOCK_ROOT/mock_ln_30M
METAWRAP_SIF=$METAWRAP_SIF CHECKM_DB=$CHECKM_DB RAPDTOOL_SIF=$RAPDTOOL_SIF \
  scripts/run_metawrap.sh zymo   $ZYMO_DIR/asm_dir     # a dir with asm/final.contigs.fasta + reads_R{1,2}.fastq
# -> results/metawrap/<ds>.checkm.stats   (MetaWRAP's native CheckM numbers)
#    results/metawrap/<ds>.micomplete.tab (common-evaluator re-scoring; the comparable numbers)
#    results/metawrap/summary.csv         (wall time, peak RSS per step)
```

> Only `binning`/`bin_refinement` are run (CheckM DB, 1.4 GB). MetaWRAP's read taxonomy
> module *is* Kraken2 (§4c gap applies by identity); naming its bins would need
> `classify_bins` (NCBI_nt, +71 GB) — not run, but the footprint is reported in §4d.

### Step 6 — Fill the tables

Transfer the medians (Step 1), OPAL metrics (Step 4), miComplete + bin counts (Step 5),
the MetaWRAP MAG comparison (Step 5d), and DB sizes into your manuscript results tables.
The genome-recovery / type-material rows (genomes recovered,
completeness/redundancy, type-material placement, novel-taxon resolution,
per-species FASTA) are the columns only RaPDTool fills — the core argument.

---

## Reproducibility checklist (keep with the analysis)

- [ ] `taxdb.version` (pinned dump source + SHA-256)
- [ ] `data/census_full.tsv` + the census figure
- [ ] `data/mock_genomes.list`, `data/mock_abundance.txt`, `data/mock_abundance_equalcov.txt`
      (the fixed vectors — without these the mocks are not reproducible)
- [ ] `summary.csv` + `db_sizes.csv` **per dataset** (`bench_<dataset>/`)
- [ ] all `*.profile` CAMI files + each dataset's `mock_composition.tsv`
- [ ] tool versions and exact command lines
- [ ] hardware spec (CPU, cores, RAM) and `THREADS`
- [ ] OPAL output directory per dataset (AMBER not applicable — see Step 5)
- [ ] which dataset the reported resource medians came from, and at what depth

## Environment recorded for this run

| | |
|---|---|
| Hardware | Intel Core i9-14900, 24 physical cores / 32 threads, 125 GB RAM, NVMe |
| Threads | 16 |
| RaPDTool | v2.3.0 (Apptainer SIF; bundles FOCUS, MetaBAT2, Binning_refiner, miComplete/Bact105, Mash, KronaTools) |
| Kraken2 | 2.17.1 · DBs: standard / 16 GB / 8 GB, build 2026-06-26 |
| Bracken | 3.0.1 |
| MetaPhlAn | 4.2.5 · DB vJan26_CHOCOPhlAnSGB |
| MEGAHIT | 1.2.9 |
| InSilicoSeq | 2.0.1 |
| OPAL | 1.0.14 |
| MetaWRAP | 1.3.0 (biocontainer `metawrap-mg:1.3.0--hdfd78af_1`; bundles CheckM, DB `checkm_data_2015_01_16`) |
| Taxonomy | NCBI taxdump 2026-07-10, sha256 `c1b91199…` (the dump the FOCUS DB was built from) |

## Known gaps

Documented in full in `benchmark_rationale.md` §4. In short: **no independent dataset**
(every genome comes from RaPDTool's own curated set, including the control half),
**no binning-accuracy metric** (AMBER needs a contig→genome truth table that
`make_mock.sh` does not produce), and a **20-genome community** with no strain
variation. Adding one or two CAMI II samples would close the first two at once — that
decision is deferred, not dismissed.

## License

Dual-licensed by artifact type:

- **Code** — everything under `scripts/` — is under the **MIT License** (`LICENSE`).
- **Data, result tables and figures** — `data/`, `results/`, `figures/` — are under
  **CC-BY-4.0** (`LICENSE-CC-BY-4.0.md`).

## Citing

If you use this kit, please cite the accompanying publication and this repository
(archived at Zenodo, DOI: 10.5281/zenodo.21528297):

> The RaPDTool authors. RaPDTool: type-material–anchored, genome-resolved
> metagenomics on a laptop. *Bioinformatics* (under review). Benchmark kit:
> https://github.com/BioTools-Dev/rapdtool-benchmark
