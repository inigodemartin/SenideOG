# SenideOG test data

`make_test_data.py` (seed 42) generates `test_senideog_root/`, a tiny
FANTASIA_project-shaped tree with 6 synthetic Viridiplantae species (~25
genes each), exercising the three real proteome sources the pipeline sees on
the server:

| Species | Source | Exercises |
|---|---|---|
| Arabidopsis_thaliana | `longest_isof` | baseline, already-clean proteome |
| Physcomitrium_patens | `longest_isof` | baseline |
| Chlamydomonas_reinhardtii | `longest_isof` | baseline, outgroup analog |
| Amborella_trichopoda | `longest_isof` | baseline |
| Oryza_sativa | `cdhit100` (no GFF) | `Source=cdhit100` covariable, no isoform collapse |
| Solanum_lycopersicum | `cdhit100` + GFF | AGAT longest-isoform collapse path (2 isoforms/gene in fixture) |

Each species also gets a `04_FunctionalAnnotation/FANTASIA_2025/*_GOs_merged.tsv`
GO table, and a fabricated OrthoFinder output (`test_senideog_root_of/Results_test/
Orthogroups/{Orthogroups.tsv,Orthogroups.GeneCount.tsv}`) standing in for what
M5/M6 would produce, so M7/M8 are testable without OrthoFinder installed.

## Running the checks

```bash
conda activate senideog

# pure-Python self-check (M1 discovery, species_code5, M7 matrix, M8
# annotation + IC-bridge-file compatibility with PCA/) -- no external
# bioinformatics tools required:
python3 test/test_senideog.py

# M1 only, against the real CLI:
python3 scripts/senideog.py --root test/test_senideog_root --output test_run/ \
    --skip_mod2 --skip_mod3 --skip_mod4 --skip_mod5 --skip_mod6 --skip_mod7 --skip_mod8

# Full quicktest once seqkit/agat/busco/orthofinder are on PATH (M3/M4 skipped
# since BUSCO lineage DBs and stratified-selection lineage data aren't part
# of this tiny fixture -- --core-species-override stands in for M4):
python3 scripts/senideog.py --root test/test_senideog_root --output test_run/ \
    --skip_mod3 --skip_mod4 \
    --core-species-override test/core_override.tsv \
    --threads 2
```

`core_override.tsv` (not bundled — create it with a `Species` column listing
3 of the 6 fixture species) stands in for M4's stratified selection, which
needs real BUSCO scores and `PCA/data/species_lineage.tsv` lineage rows this
synthetic fixture doesn't have.

Expected: <2 min end to end, `results/mod07_og_matrix.tsv` is 6 species x
(<=12) orthogroups, `results/mod08_og_desc_ic.tsv` parses cleanly with
`PCA/scripts/general_pca_common.load_go_ic_and_descriptions`.
