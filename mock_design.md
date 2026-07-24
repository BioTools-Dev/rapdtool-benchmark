# Mock community design (novel-taxon benchmark)

A 20-genome synthetic community of **type-material** genomes drawn from the
`focus_build` DB, split into two functional halves so the benchmark separates a
*fair-fight baseline* from the *capability argument*:

- **10 reference** — species present in **all** competitor DBs (Kraken2 standard,
  standard16, standard8, **and** MetaPhlAn4). All tools should classify these to
  species; they anchor the profiling-accuracy comparison so any RaPDTool win/loss
  on the conflictive half is credible. Composed of **5 clinical** (recognizable
  pathogens: *E. coli, S. aureus, P. aeruginosa, M. marinum, H. cinaedi*) + **5
  phylum-diversifiers** (*Gloeothece, Borreliella, Deinococcus, Fusobacterium,
  Marinitoga*) so the baseline is broad, not just a clinical panel.
- **10 conflictive** — species **absent from Kraken2 (full/16/8) and from
  MetaPhlAn4**, but present in RaPDTool's curated type-material DB. Read/marker
  classifiers must stop at a higher rank (LCA) or miss them; RaPDTool places them
  to species. This is the novel-taxon experiment.

The 20 genomes span **15 phyla**, including one archaeon (*Haloarcula sediminis*),
so the mock resembles a real diverse community rather than a clinical panel.

## The population these genomes were sampled from (answers the obvious objection)

The obvious reviewer objection is circularity: *you chose genomes you already knew
were missing from the competitor databases, then reported that the competitors miss
them.* The defence is that those genomes are not exotic picks — they are a
**stratified sample of a large, quantified population**, and the census that
establishes it is shipped with this kit.

Censusing **all 30,209** type-material genomes in `focus_build/db`
(`check_representation.py --sample 30209 --seed 42 -o data/census_full.tsv`):

| Absent from | genomes | % of type material |
|---|---:|---:|
| Kraken2 standard (103.7 GB) | 16,854 | **55.8 %** |
| Kraken2 capped-8 (8.1 GB) | 16,856 | 55.8 % |
| MetaPhlAn4 (59.8 GB) | 5,226 | 17.3 % |
| **Both Kraken2 standard and MetaPhlAn4** | **3,423** | **11.3 %** |

Two findings stand on their own, independent of which 20 genomes went into the mock:

1. **Kraken2's 103.7 GB standard database does not contain 55.8 % of characterised
   type material.** Capping it to 8.1 GB changes that number by <0.01 pp — capping
   removes minimizers, not taxa.
2. The gap is **taxonomically pervasive**, not confined to obscure lineages: 20.9 %
   in Chloroflexota, 19.1 % in Cyanobacteriota, 15.6 % in Actinomycetota, and 9.8 %
   even in Pseudomonadota, the best-represented phylum (see `figures/census.png`).

So the 10 conflictive genomes are a **phylum-stratified sample (one per phylum) of
the 3,423-genome absent-from-both population**, and the 10 reference genomes are a
control drawn from the 11,552 present in both. That is standard sampling
methodology, not cherry-picking — and it is what should be written in the Methods,
rather than "we selected genomes absent from the competitor databases".

Figure: `figures/census.png` / `.pdf`, from `plot_census.py`.

## How the split was decided (reproducible)

Not guessed — measured with **`check_representation.py`**, which reads each DB's
own taxonomy and reports, per genome's species taxid, presence in each competitor:

```bash
# the census that defines the population
python3 scripts/check_representation.py --sample 30209 --seed 42 -o data/census_full.tsv
# the smaller run the 20 genomes were actually drawn from
python3 scripts/check_representation.py --sample 1500 --seed 42 -o data/bench_out_repcheck.tsv
```

- **Kraken2**: species taxid ∈ `inspect.txt` (col 5) of `standard/`, `standard16/`,
  `standard8/`. The three DBs share almost the same taxonomy (59098 / 58821 / 58614
  taxids), so *absence from Kraken is judged on the full DB* and holds for all three;
  the full-vs-8 gap is a run-time **sensitivity** effect (fewer minimizers), not a
  taxid-presence effect, and is what the OPAL curve will actually quantify.
- **MetaPhlAn4**: species taxid ∈ the `.pkl` `taxonomy` values' NCBI taxid path
  (last element), with the `s__` species name as a fallback (26864 taxids / 69477
  names indexed).
- Genome → species taxid → name/phylum came from `acc_taxid_strain.tsv` +
  `taxid_lineage.tsv` (one shared taxid namespace — verified 100% overlap); each
  genome's GCA was resolved back from its taxid via `acc_taxid_strain.tsv`. Same
  pinned dump the FOCUS DB was built from.

Pool from the 1500-genome sample: 504 reference, 198 conflictive-both,
88 conflictive-mpa-only, 710 conflictive-kraken-only. The 10 conflictive picks are
all `conflictive-both` (0/0/0 in Kraken, N in MetaPhlAn), one per phylum. The 5
clinical references were added deliberately (the random sample under-covers them)
and each verified reference; the 5 diversifiers were taken from reference-verdict
genomes in phyla **not** already present in the conflictive set.

## The 20 genomes

| # | Species | taxid | Phylum | Kraken f/16/8 | MetaPhlAn | Role |
|---|---------|-------|--------|:---:|:---:|------|
| 1 | Escherichia coli | 562 | Pseudomonadota | Y/Y/Y | Y | reference · clinical |
| 2 | Staphylococcus aureus | 1280 | Bacillota | Y/Y/Y | Y | reference · clinical |
| 3 | Pseudomonas aeruginosa | 287 | Pseudomonadota | Y/Y/Y | Y | reference · clinical |
| 4 | Mycobacterium marinum | 1781 | Actinomycetota | Y/Y/Y | Y | reference · clinical |
| 5 | Helicobacter cinaedi | 213 | Campylobacterota | Y/Y/Y | Y | reference · clinical |
| 6 | Gloeothece verrucosa | 497965 | Cyanobacteriota | Y/Y/Y | Y | reference · diversity |
| 7 | Borreliella andersonii | 42109 | Spirochaetota | Y/Y/Y | Y | reference · diversity |
| 8 | Deinococcus aetherius | 200252 | Deinococcota | Y/Y/Y | Y | reference · diversity |
| 9 | Fusobacterium massiliense | 1852365 | Fusobacteriota | Y/Y/Y | Y | reference · diversity |
| 10 | Marinitoga piezophila | 443254 | Thermotogota | Y/Y/Y | Y | reference · diversity |
| 11 | Streptomyces coeruleoprunus | 285563 | Actinomycetota | N/N/N | N | conflictive |
| 12 | Chitinimonas violaceus | 3459088 | Pseudomonadota | N/N/N | N | conflictive |
| 13 | Lactimicrobium massiliense | 2161814 | Bacillota | N/N/N | N | conflictive |
| 14 | Paraniabella aurantiaca | 3393740 | Bacteroidota | N/N/N | N | conflictive |
| 15 | Deferribacter thermophilus | 53573 | Deferribacterota | N/N/N | N | conflictive |
| 16 | Ktedonobacter racemifer | 485913 | Chloroflexota | N/N/N | N | conflictive |
| 17 | Haloarcula sediminis | 3111777 | Archaea (Methanobacteriota) | N/N/N | N | conflictive |
| 18 | Helicobacter zhangjianzhongii | 2974574 | Campylobacterota | N/N/N | N | conflictive |
| 19 | Corallococcus praedator | 2316724 | Myxococcota | N/N/N | N | conflictive |
| 20 | Rosettibacter firmus | 3111522 | Ignavibacteriota | N/N/N | N | conflictive |

> Rows 5 & 18 are both *Helicobacter*: one all tools resolve, one only RaPDTool —
> a within-genus illustration that Kraken/MetaPhlAn would collapse the novel species
> to the genus while RaPDTool separates it. Likewise rows 4 & 11 (Actinomycetota).

## Honesty notes (write these into the paper, don't hide them)

- **State the limitation in one explicit sentence, before a reviewer does:** *the
  conflictive taxa come from RaPDTool's own curated database, so their detection is
  expected by construction; the experiment demonstrates database coverage of type
  material, not algorithmic superiority.* With the census above, that sentence costs
  nothing and buys credibility — the claim being made is about database design, and
  a coverage gap cannot be measured with taxa that are present in every database.
- The reference half is the control: all tools should detect it, and at adequate depth
  RaPDTool does — **10/10 reference at 3/10/30 M** (OPAL, `mock_ln_*`). An earlier
  "RaPDTool loses on the reference half, 5/10" note came from the pathological 2 M mock
  (one genome held 61 % of reads, starving the rest) and **is superseded** — do not
  repeat it. The credibility of the conflictive result now rests on the census and the
  symmetric mirror experiment, not on a reference-half loss.
- The conflictive genomes come from RaPDTool's own curated DB — **that is the point**
  (a small type-material DB covers taxa the big general DBs miss), not a trick. State
  it plainly: the advantage is *curation for type material*, not a larger index.
- `make_mock.sh` simulates reads (InSilicoSeq) and **re-assembles** (MegaHit), so
  RaPDTool classifies reconstructed contigs, not the reference file it was built from
  — this blunts the "it's literally its own input" objection but does not remove that
  these organisms are in its DB.
- Complement with an **external CAMI II dataset** (organisms RaPDTool did *not*
  curate) for the fair-fight accuracy claim; use this in-house mock for the
  footprint + novel-taxon capability claims.

## Reproduce

```bash
# 1. build the mock (reads + contigs + gold standard)
scripts/make_mock.sh -o mock -t 16 -l data/mock_genomes.list        # mock lives inside the kit
# 2. set run_benchmark.sh CONFIG: ASSEMBLY/READS_R1/READS_R2 from ./mock,
#    KRAKEN_DB_FULL=$KRAKEN_DB_FULL (and CAP16/CAP8 — from config.sh),
#    MPA_DB=$MPA_DB
# 3. gold standard -> CAMI:
#    python3 scripts/profile2cami.py mock/mock_composition.tsv -f truth -o results/bench_out/gold_standard.profile -s gold
```
