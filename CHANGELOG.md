# Changelog

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
