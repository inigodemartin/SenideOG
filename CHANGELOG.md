# Changelog

## [v0.8.7] — 2026-09-17

### Fixed
- `find_orthofinder_results` only looked at direct `Results_*` children of
  the given directory. Every "-b" resume of an interrupted Module 5 run
  makes OrthoFinder nest its next attempt one level deeper
  (`Results_X/WorkingDirectory/OrthoFinder/Results_Y/...`), so after any
  resume the function kept returning the outer, incomplete shell instead of
  the real completed `Orthogroups/` output. That broken `core_results` was
  then passed to Module 6's `--core`, which OrthoFinder didn't recognize as
  valid prior orthogroups — instead of a fast assignment, it silently
  redid the full de novo clustering (BLAST, MCL, MSA/gene trees) for all
  species combined. `find_orthofinder_results` now searches recursively
  and requires `Orthogroups/` to exist; a new `_find_latest_orthofinder_attempt`
  helper (any depth, complete or not) is used separately to locate a
  `WorkingDirectory` to resume from. `run_module6`'s before/after directory
  diff (added in v0.8.6) now also searches recursively, rooted at the
  stable `of_core_dir` instead of `core_results.parent`.

## [v0.8.6] — 2026-09-17

### Fixed
- `run_module6` passed `-o of_assign_dir` to `orthofinder --assign`, which
  OrthoFinder rejects (`-o` is only valid with a fresh `-f` run) — Module 6
  failed after Module 5's 11h core inference had already completed.
  `--assign` always writes its `Results_*` dir as a sibling of `--core`
  instead; the fix snapshots that directory before/after the run and moves
  the new one under `of_assign/`, so downstream lookups are unaffected.

## [v0.8.5] — 2026-09-14

### Fixed
- `envs/senideog.yaml` listed `treepl`, which isn't packaged on
  bioconda/conda-forge — `conda env create` failed to solve on every
  fresh install. It's only needed later by `senideog_phylo.py` (fase 2)
  for divergence dating, not by the current M1-M8 pipeline; removed from
  the yaml with a note to build it from source
  (https://github.com/blackrim/treePL) when that script exists.

## [v0.8.4] — 2026-09-14

### Fixed
- `run_module5`/`run_module6` still hit OrthoFinder's `ERROR: non-default
  output directory already exists` if a *previous* attempt had been
  interrupted mid-run, leaving `of_core`/`of_assign` on disk with no
  `Orthogroups/` and no resumable `Blast*.txt.gz` (v0.8.3 only fixed the
  unconditional pre-creation on a clean first run). Both functions now
  remove the leftover directory before launching a fresh OrthoFinder run.

## [v0.8.3] — 2026-09-11

### Fixed
- M5 (`run_module5`) pre-created `of_dir` (`workdir/of_core/`) before
  launching a fresh OrthoFinder run, but OrthoFinder's `-f` mode requires
  its `-o` target to not exist yet (it creates it itself) — every first
  attempt failed immediately with `ERROR: non-default output directory
  already exists`. `of_dir` is no longer pre-created; it was never
  needed for anything else in the function (`core_proteomes`, the actual
  input, lives under `of_dir.parent`).

## [v0.8.2] — 2026-09-11

### Fixed
- M2's duplicate-id uniquification (v0.6.0) could reintroduce the exact
  collision it was meant to fix: a repeated header like `gene1` was
  renamed to `gene1_2` without checking whether the raw proteome already
  had its own distinct `gene1_2` entry — in that case the renamed and the
  pre-existing header collided again, and BUSCO still hard-errors on the
  duplicate id (`BuscoError: Duplicate of sequence`), aborting the whole
  BUSCO run for that species (recorded as `NO_BUSCO_RESULT`, with the
  output directory removed entirely, so no summary file is left behind).
  The chosen suffix now skips any candidate that collides with an id
  already present in the raw proteome or already assigned.

## [v0.8.1] — 2026-09-11

### Fixed
- M3's BUSCO checkpoint (`run_busco`) only checked that
  `short_summary.*.txt` exists and is non-empty. BUSCO writes the header
  first and the `***** Results: *****` block last, so a run killed
  mid-execution (OOM, timeout, scheduler preemption under
  `--busco-jobs` parallelism) leaves a non-empty but incomplete summary
  file, which the checkpoint mistook for a finished run — permanently
  recording that species as `NO_BUSCO_RESULT` with no retry. The
  checkpoint now also requires the `C:<pct>%` results token to be present
  before trusting an existing summary; species with a truncated summary
  are rerun automatically on the next launch.

## [v0.8.0] — 2026-09-11

### Added
- `mod03_qc.tsv` now includes every species from the M1 manifest, not just
  the ones that made it to M2/M3 — species dropped for `no_gca_found`,
  `no_proteome_found`, `ambiguous_gca_skipped` or `skipped_by_user` are
  added back in with a new `Status` column carrying that M1 reason and
  `BUSCO_status="N/A"`, so the full initially-selected species list can be
  monitored in one table instead of only the survivors.

## [v0.7.1] — 2026-09-11

### Fixed
- M3's BUSCO results-line parsing depended on the brackets closing right
  after `M:<pct>%`. BUSCO 6.0.0 moved the closing bracket to right after
  `D:<pct>%` and put `F`/`M`/`n` after it
  (`C:98.1%[S:57.3%,D:40.8%],F:1.0%,M:1.0%,n:822` vs. the old
  `C:98.4%[S:97.2%,D:1.2%,F:0.7%,M:0.9%,n:425]`), so the v0.7.0 regex fix
  still didn't match and every species kept coming out `NO_BUSCO_RESULT`
  despite the `short_summary.*.txt` files being present and correct.
  Each field (`C`/`S`/`D`/`F`/`M`) is now looked up independently by its
  `<letter>:<pct>%` token instead of relying on bracket position, so it no
  longer depends on which BUSCO version produced the summary.

## [v0.7.0] — 2026-09-11

### Fixed
- M3's BUSCO short-summary regex required the `]` to close immediately
  after `M:<pct>%`, but real BUSCO output is `M:<pct>%,n:<count>]` — the
  regex never matched, so every species (not just ones with a missing
  summary file) was recorded as `NO_BUSCO_RESULT` in `mod03_qc.tsv` even
  after a clean, successful BUSCO run. Regex now allows the optional
  `,n:<count>` before the closing bracket.

### Added
- `mod03_qc.tsv` gained `N_proteins` (carried over from M2's proteome
  stats) and `Mean_protein_length` (mean aa length of the cleaned
  proteome) columns, between `Code5` and `BUSCO_C`.

## [v0.6.0] — 2026-09-11

### Fixed
- M2 (`clean_and_prefix_fasta`) now uniquifies duplicate sequence ids within
  a single species' raw proteome (`_2`, `_3`, ... suffix on repeats) —
  BUSCO hard-errors on a repeated id (`BuscoError: Duplicate of sequence`),
  and OrthoFinder would otherwise silently conflate the two proteins.
  Real-world trigger: a proteome with a duplicated organelle gene entry.
- `run_busco` no longer aborts the whole run when BUSCO fails for one
  species (previously `_run()`'s `sys.exit(1)` on a non-zero exit code
  took every other species down with it — including all the ones already
  running under `--busco-jobs`). The failing species is now recorded as
  `NO_BUSCO_RESULT` and the run continues; `_run()` gained a `check=False`
  escape hatch for this, every other caller (seqkit/OrthoFinder/AGAT)
  still fails fast as before. Its BUSCO output directory is removed on
  failure too (nothing worth keeping without a summary file).

## [v0.5.0] — 2026-09-11

### Added
- `--busco-jobs` (default 1): runs M3's BUSCO calls for that many species in
  parallel (each still uses `--threads` internally, so total CPU used ~=
  `--busco-jobs` x `--threads`). The lineage dataset is pre-downloaded once
  via `busco --download` before dispatching, avoiding a race where several
  parallel `busco` processes try to populate the shared `busco_downloads/`
  cache at once on a first run. If one species' BUSCO run fails, jobs not
  yet started are cancelled instead of working through the full species
  list before the failure surfaces.

### Changed
- Each species' BUSCO output directory (`workdir/busco/<Code5>/`) is pruned
  down to just the `short_summary.*.txt` file right after it's parsed —
  `hmmer_output/`, `busco_sequences/`, `logs/`, etc. are large and fully
  reproducible by rerunning BUSCO, so they no longer pile up per species.
  The kept summary file is also what the existing per-species checkpoint
  looks for, so already-completed species are still recognized and skipped
  on a resumed run.

## [v0.4.0] — 2026-09-11

### Changed
- `data/species_lineage.tsv` and `data/species_taxid.tsv` are now bundled
  inside this repo instead of being read from a sibling `PCA/` checkout —
  `--lineage-file` and `--species-taxid` default to `data/` here, so
  `--kingdom` filtering and M4 core selection work out of the box after a
  plain `git clone`, with no other repo required alongside it.

## [v0.3.0] — 2026-09-11

### Added
- `--kingdom` (default `Viridiplantae`): with `--fasta-dir`, restricts the
  flat proteome folder to species whose `kingdom` column in `--lineage-file`
  matches — the folder may hold non-target species (e.g. fungi) alongside
  Viridiplantae ones. Prints the number of selected species to stdout.

## [v0.2.0] — 2026-09-10

### Added
- `--fasta-dir`: alternative to `--root` for M1 when proteomes are already
  collected in one flat folder (no `00_GenomeSource/<GCA>/` tree, no GO
  tables). `senideog_common.collect_proteome_files_from_dir` derives
  `Species` from the filename and reuses the `PROTEOME_TIERS` patterns to
  tag `Source`. GO-dependent M8 has nothing to annotate with in this mode
  (pass `--skip_mod8`).

## [v0.1.0] — 2026-09-10

### Added
- `scripts/senideog.py` (M1-M8): inventory + collection, normalization
  (one isoform/gene, seqkit cleaning, `<Code5>|` id prefix), QC (BUSCO,
  FASTA<->GO id match rate, disk check), stratified core-species selection
  (APG IV order/family), OrthoFinder de novo inference on the core,
  OrthoFinder `--assign` for the rest, species x orthogroup matrix, and
  GO-consensus orthogroup annotation with an IC-bridge file compatible with
  `PCA/scripts/general_pca_abundance.py --ic-file`.
- `scripts/senideog_common.py`: shared logging/subprocess infra, proteome
  discovery, FASTA cleaning, QC, core selection, and OrthoFinder-output
  parsing.
- `envs/senideog.yaml` conda environment.
- `test/`: 6-species synthetic FANTASIA_project-shaped fixture
  (`make_test_data.py`, seed 42) and `test_senideog.py` self-check covering
  M1/M7/M8 without requiring OrthoFinder/BUSCO/AGAT/seqkit installed.
- FAIR compliance files: `LICENSE` (MIT), `CITATION.cff`.

### Not included yet
- `scripts/senideog_phylo.py` (fase 2: IQ-TREE/ASTRAL species tree, treePL
  ultrametric calibration, CAFE5 gene-family expansion/contraction) —
  deferred until M5/M6 produce real outputs to build and test against.
