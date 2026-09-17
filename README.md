<img src="https://img.shields.io/badge/version-v0.8.9-teal"/>
<img src="https://img.shields.io/badge/python-3.10%2B-blue"/>
<img src="https://img.shields.io/badge/platform-Linux%20%7C%20macOS-lightgrey"/>

See [CHANGELOG.md](CHANGELOG.md) for release notes.

## SenideOG

Builds a **species x orthogroup matrix** for ~500 Viridiplantae reference
proteomes, replacing the GO-term unit used by the existing comparative
analysis (`PCA/`) with orthology inferred directly from the sequences —
copy-number, phylogenomics, gene-family dynamics, all without the embedding-
similarity annotation bias GO terms carry (documented in `PCA/scripts/go_phylostrata.py`).

**Name**: `senide` is Basque for "relative / kin" (no borrowed etymology) —
an orthogroup groups genes descended from a common ancestor, the same way
`senide` groups relatives. `OG` = orthogroup, per this workspace's naming
convention.

### Overview

| Module | Tool | Output |
|---|---|---|
| M1 Inventory | filesystem scan (`--root` tree or flat `--fasta-dir`) | `results/mod01_manifest.tsv` |
| M2 Normalization | AGAT, seqkit | `workdir/proteomes_clean/<Code5>.fa`, `results/mod02_proteome_stats.tsv` |
| M3 QC | BUSCO | `results/mod03_qc.tsv`, `results/mod03_busco.{pdf,png}` |
| M4 Core selection | stratified sampling (APG IV) | `results/mod04_core_species.tsv` |
| M5 De novo inference | OrthoFinder (core) | `workdir/of_core/Results_*/` |
| M6 Assignment | OrthoFinder `--assign` (SHOOT) | `workdir/of_assign/Results_*/` |
| M7 OG matrix | — | `results/mod07_og_matrix.tsv`, `results/mod07_species_stats.tsv` |
| M8 Annotation | GO consensus over FANTASIA tables | `results/mod08_og_annotation.tsv`, `results/mod08_og_desc_ic.tsv` |

Fase 2 (`scripts/senideog_phylo.py`, not yet written — needs M5/M6 outputs to
build and test against): IQ-TREE/ASTRAL species tree, treePL ultrametric
calibration, CAFE5 gene-family expansion/contraction.

### Requirements

```bash
conda env create -f envs/senideog.yaml
conda activate senideog
```

| Tool | Used by |
|---|---|
| pandas | all modules |
| AGAT | M2 (longest-isoform collapse when the source proteome isn't already `longest_isof`) |
| seqkit | M2 (sequence cleaning) |
| BUSCO >=5 | M3 |
| OrthoFinder >=3.1.5 | M5, M6 |
| diamond | M5, M6 (via OrthoFinder) |
| codecarbon (optional) | carbon footprint tracking |
| trimal, IQ-TREE, ASTRAL, treePL, CAFE5 | fase 2 (`senideog_phylo.py`) |

### Usage

```bash
scripts/senideog.py --root /data/users/demartini/FANTASIA_project --output run01 [options]

# or, with proteomes already collected in one flat folder (no GO tables => skip M8):
scripts/senideog.py --fasta-dir /path/to/proteomes_fasta --output run01 --skip_mod8 [options]
```

| Option | Default | Description |
|---|---|---|
| `--root` | — | FANTASIA_project root (required unless `--fasta-dir` or `--skip_mod1`) |
| `--fasta-dir` | — | Flat folder with one proteome FASTA per species, no GCA/GO tree — alternative to `--root`. No GO tables in this mode, so pass `--skip_mod8` |
| `--output` | — | run output directory (required) |
| `--threads` | 8 | threads for seqkit/BUSCO/OrthoFinder |
| `--species-taxid` | `data/species_taxid.tsv` (bundled) | Species -> TaxID lookup |
| `--non_interactive` | off | skip (rather than prompt for) species with >1 genome assembly |
| `--kingdom` | `Viridiplantae` | With `--fasta-dir`: keep only species whose `kingdom` column in `--lineage-file` matches this (pass `''` to disable) |
| `--min-protein-len` | 30 | M2: drop proteins shorter than this |
| `--busco-lineage` | `viridiplantae_odb12` | M3 BUSCO dataset |
| `--busco-jobs` | 1 | M3: BUSCO runs to execute in parallel (each still uses `--threads`; total CPU ~= `--busco-jobs` x `--threads`) |
| `--busco-c-pass` / `--busco-c-flag` | 0.85 / 0.80 | M3 BUSCO C%% thresholds (PASS / FLAG / FAIL) |
| `--id-match-threshold` | 0.5 | M3: min. FASTA<->GO-table id match rate before aborting |
| `--disk-estimate-gb` | 300 | M3: disk space the core all-vs-all + gene trees are expected to need |
| `--core-size` | 64 | M4 target core size |
| `--lineage-file` | `data/species_lineage.tsv` (bundled) | M4 order/family (APG IV) lookup for stratified selection; also used by M1's `--kingdom` filter |
| `--core-species-override` | — | TSV with a `Species` column to use as the core directly (bypasses M4) |
| `--min-species` | 4 | M7: drop orthogroups present in fewer species than this |
| `--go-consensus-threshold` | 0.5 | M8: min. fraction of annotated OG members sharing a GO to report it |
| `--skip_mod1` .. `--skip_mod8` | off | skip a module, reusing its checkpointed output |
| `--force` | off | rerun all steps regardless of existing outputs |
| `--dry_run` | off | validate inputs, print planned steps, exit |
| `--disable_co2_tracking` | off | disable codecarbon even if installed |

### Output layout

```
{output}/
├── results/     mod01_manifest.tsv ... mod08_og_desc_ic.tsv, {prefix}.run_summary.json
├── workdir/     proteomes_clean/, of_core/, of_assign/ (safe to delete after success)
└── logs/        Run_SenideOG_{date}.log, {prefix}.emissions.csv
```

Every run writes `logs/Run_SenideOG_{YYYYMMDD}[_N].log` (append mode — a
resumed multi-day run keeps its full history) and
`results/{prefix}.run_summary.json` (tool versions via VERSION, species
counts per module, resource usage, carbon footprint if tracked).

### Recommended workflow

```bash
# 1. dry run
scripts/senideog.py --root /data/.../FANTASIA_project --output run01 --dry_run

# 2. inventory + QC only, and REVIEW results/mod03_qc.tsv BY HAND before continuing
scripts/senideog.py --root /data/.../FANTASIA_project --output run01 \
    --skip_mod5 --skip_mod6 --skip_mod7 --skip_mod8

# 3. the long job (setsid, not tmux -- a dead tmux session kills the job)
setsid nohup scripts/senideog.py --root /data/.../FANTASIA_project --output run01 \
    --threads 64 >> run01/logs/nohup.out 2>&1 &

# 4. resume after a crash = relaunch the SAME command (checkpoints + OrthoFinder's own -b/--core resume)

# 5. feed the matrix into the existing PCA, unmodified:
cd ../PCA
scripts/general_pca_abundance.py \
    -m  ../SenideOG/run01/results/mod07_og_matrix.tsv \
    --species-stats     ../SenideOG/run01/results/mod07_species_stats.tsv \
    --taxonomy           merged_taxons_belen.tsv \
    --ic-file            ../SenideOG/run01/results/mod08_og_desc_ic.tsv \
    --metazoa-taxonomy   /dev/null \
    --output og_pca_viridiplantae.html

# 6. fase 2: scripts/senideog_phylo.py (once written)
```

Checkpoint granularity matches rerun cost: M1-M4 are minutes, **M5 and M6 are
the only ones that matter (days)** and both delegate their own resume logic
to OrthoFinder's `-b`/`--core`.

### Third-party tools

OrthoFinder ([Emms & Kelly 2019](https://doi.org/10.1186/s13059-019-1832-y)),
DIAMOND ([Buchfink et al. 2021](https://doi.org/10.1038/s41592-021-01101-x)),
BUSCO ([Manni et al. 2021](https://doi.org/10.1093/molbev/msab199)),
AGAT ([Dainat 2022](https://doi.org/10.5281/zenodo.3552717)), seqkit
([Shen et al. 2016](https://doi.org/10.1371/journal.pone.0163962)).
