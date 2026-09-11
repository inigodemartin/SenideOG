# Changelog

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
