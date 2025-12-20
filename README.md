# BiTUGA: Scalable Prevalence-Based Unitig Association Testing for binary traits

**BiTUGA** is a scalable workflow designed for **prevalence-based unitig association testing**, specifically optimized for **binary traits** in large and highly repetitive genomes.

Unlike traditional reference-free methods that rely on computationally expensive k-mer-by-sample matrices, BiTUGA shifts the statistical focus from raw abundance to **population-level prevalence**. It leverages **KMC** and **BCALM** to assemble differentially distributed k-mers into unitigs, allowing for the identification of trait-associated genomic regions (e.g., Sex-Determining Regions) without requiring a reference genome.

## Citation

If you use BiTUGA in your research, please cite:

> **Mittelbach, J., Kersten, B., & Kurtz, S. (2025).** *BiTUGA: scalable prevalence-based unitig association testing for binary traits in large genomes.* Manuscript submitted to **G3: Genes|Genomes|Genetics**.

BibTeX:
```bibtex
@unpublished{BiTUGA2025,
  author = {Mittelbach, Jannes and Kersten, Birgit and Kurtz, Stefan},
  title  = {BiTUGA: scalable prevalence-based unitig association testing for binary traits in large genomes},
  year   = {2025},
  note   = {Manuscript submitted to G3: Genes|Genomes|Genetics},
  url    = {https://github.com/JMittelbach/BiTUGA},
  version = {v1.0.0}
}
```

## Quick Start
```bash
# 1) Clone with submodules
git clone --recursive https://github.com/JMittelbach/BiTUGA.git
cd BiTUGA

# 2) Build native components
make

# 3) Install plotting dep (required for plots)
pip install matplotlib

# 4) Inspect options
./BiTUGA.sh --help

# Minimal execution example
./BiTUGA.sh \
  --input-dir data/reads \
  --trait-info trait_info/trait_info_populus.tsv \
  --out-dir runs/run1 \
  --threads 32 \
  --mem-gb 100

# Compiler Definitions

By default, the code is compiled using clang and clang++. You can
specify your favorite compiler by setting the environment variables
CC and CXX, like for example

export CC=/opt/homebrew/bin/gcc-15
export CXX=/opt/homebrew/bin/g++-15
```

## Conceptual Overview
- Steps: count k-mers → filter by prevalence → assemble unitigs → map to a input reads → Fisher tests + multiple testing correction → FASTA with significant unitigs. Supports raw (untrimmed) FASTQs; single-end and paired-end reads are accepted.

## Repository Structure
```
BiTUGA/
├── BiTUGA.sh          # main orchestration script
├── bin/               # helpers and compiled utilities
├── src/               # C++ sources (rebuild via make); contains external/ (KMC3/BCALM2)
├── trait_info/        # example trait tables (tsv)
├── utils/             # study-specific helper scripts (not core pipeline)
```

## Key parameters (BiTUGA.sh) — see `./BiTUGA.sh --help+` for all
- `--input-dir`        FASTQ directories  
- `--trait-info`       trait table (tsv/csv/txt; sample ID + trait)  
- `--out-dir`          output/work dir (writes `metadata.tsv` if missing)  
- `--k`                k-mer length (default 31)  
- `--ci` / `--cx`      KMC min/max counter (defaults 2 / 4,294,967,295); higher `ci` reduces noise, `cx` caps abundance  
- `--threads`          CPU threads (default 12)  
- `--mem-gb`           RAM limit (default 16 GB)  
- `--group-contrast`   presence in ≥50% of one group AND ≤30% in the other  
- Prevalence defaults   auto-adaptive: `--prev-min/--prev-max`=0.25/0.75; `--group-min/--group-max`=0.0/1.0 unless group-contrast is set  
- `--fdr` / `--fdr-alpha`   default `--fdr auto` (BH, alpha 0.05); `--raw-p-threshold` default 0.01 when FDR is off  
- `--query-split-size` MiniMatcher batch size (default 262144; use powers of two 65536–524288)  
- `--max-replicate-ref` / `--max-replicate-qry` cap replicated seeds (Pass IV)  
- `--start-pass` / `--end-pass` rerun subsets; restart from earliest affected pass after param changes  
Trait table: two columns (sample ID, trait/binary label); see `trait_info/` examples.

## Pipeline passes

- **Pass I — k-mer aggregation**  
  - Does: count pooled k-mers.  
  - Takes: FASTQs; `--k`, `--ci/--cx`.  
  - Writes: KMC DBs (`results/kmc_*`).  
  - Tuning: OOM/disk → reduce threads or raise `--ci`.  
  - Scales with: total reads and k.

- **Pass II — k-mer filtering**  
  - Does: prevalence filtering (global/adaptive, group contrast).  
  - Takes: KMC outputs; trait table `--trait-info`.  
  - Knobs: `--prev-min/--prev-max`, `--group-min/--group-max`, `--group-contrast`.  
  - Writes: filtered k-mers/FASTA (`results/kmers*`), merged stats via `merge2stats`.  
  - Tuning: over-filtering → relax thresholds.  
  - Scales with: number of unqiue k-mers.

- **Pass III — unitig assembly**  
  - Does: assemble unitigs (BCALM), optional filtering.  
  - Takes: filtered k-mers; `--k`; length/entropy filters.  
  - Writes: unitigs FASTA (`results/unitigs*.fa`).  
  - Scales with: retained k-mers.

- **Pass IV — unitig presence verification**  
  - Does: map unitigs (MiniMatcher) and compute coverage tables.  
  - Takes: unitigs FASTA, reference FASTA; `--nt-threads`, `--query-split-size` (reads/batch).  
  - Knobs: `--query-split-size` (powers of two), seed/rev-comp options, `--min-mem-length` (default k+1), replicate caps.  
  - Writes: PAF, coverage TSVs, top-window summaries, plots.  
  - Tuning: OOM/time → lower `--query-split-size`.  
  - Scales with: unitig count and reference size.

- **Pass V — Fisher's tests**  
  - Does: association testing and FASTA export of significant unitigs.  
  - Takes: coverage tables; trait table.  
  - Knobs: `--fdr` (`auto`, `bh`, `storey`, `none`), `--fdr-alpha`, `--raw-p-threshold`, prevalence deltas.  
  - Writes: results TSVs (p/q-values, counts, prevalences, bias labels), FASTA of significant unitigs, plots.  
  - Tuning: sparse tests → weak FDR; use `--fdr none` + `--raw-p-threshold` if exploratory.  
  - Scales with: number of tested unitigs.

## FDR and significance
- Modes: `--fdr auto` (default BH), `bh`, `storey`, `none` (skip FDR → use `--raw-p-threshold`).  
- Applied after prevalence/contrast filtering (Pass V).  
- Thresholds: `--fdr-alpha`; non-FDR: `--raw-p-threshold` (default 0.01).  
- Implementation: Fisher tests and BH/Storey are implemented in the pipeline code (no external Python stats deps required for Pass V).

## Prevalence / contrast filtering
- Global: `--prev-min/--prev-max` (fraction across all samples).  
- Grouped: `--group-min/--group-max` enforce presence/absence asymmetry.  
- Shortcut: `--group-contrast` sets asymmetric defaults (present in ≥50% of one group and ≤30% of the other).  
- Adaptive defaults adjust for unbalanced groups; tighten thresholds for specificity, loosen for sensitivity.

## Rerun semantics
- `--start-pass p` assumes passes < p are valid and rebuilds ≥ p with current params; downstream outputs are overwritten.  
- If earlier outputs are stale, clean them before resuming.  
- `--end-pass` can stop early after a given pass.

## Outputs (under `--out-dir`)
- KMC DBs, filtered k-mers, unitigs (`*.fa`), MiniMatcher TSV, plots (`png`/`pdf`), logs (`logging/`). 
- Results TSVs include p-value, q-value (if FDR), counts and prevalence per group, bias labels, and per-window coverage where applicable.  
- Logs capture command lines and key settings (threads, query_split_size, thresholds) for reproducibility.  
- A `metadata.tsv` is written to `--out-dir` if not provided; trait-info can be tsv/csv/txt with sample ID and trait columns.

## Dependencies
- Core pipeline: bash + compiled tools (KMC, BCALM, MiniMatcher, merge2stats). 
- Plots: `matplotlib`.  

## License
This project is licensed under the MIT License. See [LICENSE](LICENSE).
