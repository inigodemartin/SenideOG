#!/usr/bin/env python3
"""
Self-check for the parts of SenideOG runnable without OrthoFinder/BUSCO/AGAT/
seqkit installed (M1, M7, M8 pure-Python logic, and the species_code5 helper)
-- run with: python3 test/test_senideog.py
Regenerates test_senideog_root/ and test_senideog_root_of/ via make_test_data.py first.
"""

import shutil
import subprocess
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent / "scripts"))
import pandas as pd
import senideog_common as C

TEST_DIR = Path(__file__).parent
ROOT = TEST_DIR / "test_senideog_root"
OF_RESULTS = TEST_DIR / "test_senideog_root_of" / "Results_test"


def test_species_code5():
    taken = set()
    assert C.species_code5("Solanum_lycopersicum", taken) == "Slyco"
    assert C.species_code5("Solanum_lycopersicoides", taken) == "Slyc2"  # collision handled, stays 5 chars
    assert len(taken) == 2


def test_module1_inventory():
    manifest = C.collect_proteome_files(ROOT, non_interactive=True, taxid_map={})
    assert len(manifest) == 6, f"expected 6 species, got {len(manifest)}"
    assert (manifest["Status"] == "resolved").all(), manifest[manifest["Status"] != "resolved"]
    tomato = manifest[manifest["Species"] == "Solanum_lycopersicum"].iloc[0]
    assert tomato["Source"] == "cdhit100"
    assert tomato["GFF"] is not None and Path(tomato["GFF"]).exists()
    rice = manifest[manifest["Species"] == "Oryza_sativa"].iloc[0]
    assert rice["Source"] == "cdhit100"
    arabidopsis = manifest[manifest["Species"] == "Arabidopsis_thaliana"].iloc[0]
    assert arabidopsis["Source"] == "longest_isof"


def test_module1_flat_dir():
    """--fasta-dir mode: a flat folder of proteome FASTAs, no GCA/GO tree."""
    ara_fasta = next(ROOT.glob("Arabidopsis_thaliana/00_GenomeSource/*/*_genomic_longest_isof_proteins.fasta"))
    rice_fasta = next(ROOT.glob("Oryza_sativa/00_GenomeSource/*/*_uniq_cdhit100_5k_removed.pep"))
    flat_dir = TEST_DIR / "_tmp_flat_proteomes"
    flat_dir.mkdir(exist_ok=True)
    try:
        shutil.copy2(ara_fasta, flat_dir / ara_fasta.name)
        shutil.copy2(rice_fasta, flat_dir / rice_fasta.name)

        manifest = C.collect_proteome_files_from_dir(flat_dir, taxid_map={})
        assert len(manifest) == 2, manifest
        assert (manifest["Status"] == "resolved").all()
        assert manifest["GOFile"].isna().all() and manifest["GFF"].isna().all()

        ara_row = manifest[manifest["Species"].str.startswith("Arabidopsis_thaliana")].iloc[0]
        assert ara_row["Source"] == "longest_isof"
        rice_row = manifest[manifest["Species"].str.startswith("Oryza_sativa")].iloc[0]
        assert rice_row["Source"] == "cdhit100"
    finally:
        shutil.rmtree(flat_dir)


def test_module1_flat_dir_kingdom_filter():
    """--kingdom filter (used with --fasta-dir): only species tagged with
    the requested kingdom in --lineage-file are kept."""
    ara_fasta = next(ROOT.glob("Arabidopsis_thaliana/00_GenomeSource/*/*_genomic_longest_isof_proteins.fasta"))
    rice_fasta = next(ROOT.glob("Oryza_sativa/00_GenomeSource/*/*_uniq_cdhit100_5k_removed.pep"))
    flat_dir = TEST_DIR / "_tmp_flat_proteomes_kingdom"
    flat_dir.mkdir(exist_ok=True)
    lineage_tsv = TEST_DIR / "_tmp_lineage.tsv"
    try:
        shutil.copy2(ara_fasta, flat_dir / "Arabidopsis_thaliana.fasta")
        shutil.copy2(rice_fasta, flat_dir / "Oryza_sativa.fasta")
        pd.DataFrame([
            {"Species": "Arabidopsis_thaliana", "kingdom": "Viridiplantae"},
            {"Species": "Oryza_sativa", "kingdom": "Fungi"},
        ]).to_csv(lineage_tsv, sep="\t", index=False)

        kingdom_species = C.load_kingdom_species(lineage_tsv, "Viridiplantae")
        manifest = C.collect_proteome_files_from_dir(flat_dir, taxid_map={}, kingdom_species=kingdom_species)
        assert len(manifest) == 1, manifest
        assert manifest.iloc[0]["Species"] == "Arabidopsis_thaliana"
    finally:
        shutil.rmtree(flat_dir)
        lineage_tsv.unlink(missing_ok=True)


def test_clean_and_prefix_dedup():
    """Duplicate headers in the raw proteome must be uniquified, not silently
    kept as-is -- BUSCO hard-errors on a repeated sequence id."""
    workdir = TEST_DIR / "_tmp_clean_workdir"
    workdir.mkdir(exist_ok=True)
    raw = workdir / "raw.fa"
    raw.write_text(">gene1\nMAAA\n>gene1\nMBBB\n>gene2\nMCCC\n")
    out = workdir / "Test3.fa"
    orig_require_tool, orig_run = C._require_tool, C._run
    try:
        C._require_tool = lambda name: name
        C._run = lambda cmd, **k: (shutil.copy2(raw, Path(cmd[-1])), subprocess.CompletedProcess(cmd, 0))[1]
        stats = C.clean_and_prefix_fasta(raw, out, "Test3", min_len=1)
        assert stats["n_kept"] == 3
        headers = [h for h, _ in C.iter_fasta(out)]
        assert headers == ["Test3|gene1", "Test3|gene1_2", "Test3|gene2"], headers
    finally:
        C._require_tool, C._run = orig_require_tool, orig_run
        shutil.rmtree(workdir, ignore_errors=True)


def test_run_busco_nonfatal_failure():
    """A BUSCO failure for one species must not raise -- callers (M3's
    parallel dispatch) rely on this so other species aren't aborted."""
    workdir = TEST_DIR / "_tmp_busco_workdir"
    workdir.mkdir(exist_ok=True)
    orig_require_tool, orig_run = C._require_tool, C._run
    try:
        C._require_tool = lambda name: name
        C._run = lambda *a, **k: subprocess.CompletedProcess(args=[], returncode=1, stdout="", stderr="boom")
        result = C.run_busco(TEST_DIR / "nonexistent.fa", "viridiplantae_odb12", workdir, 4, "Test2")
        assert result == {"C": None, "S": None, "D": None, "F": None, "M": None}
    finally:
        C._require_tool, C._run = orig_require_tool, orig_run
        shutil.rmtree(workdir, ignore_errors=True)


def test_busco_cleanup():
    """After a BUSCO run, only the short_summary file should survive per species."""
    out_dir = TEST_DIR / "_tmp_busco_out"
    out_dir.mkdir(exist_ok=True)
    try:
        keep = out_dir / "short_summary.specific.viridiplantae_odb12.Test1.txt"
        keep.write_text("C:95.0%[S:90.0%,D:5.0%,F:2.0%,M:3.0%]\n")
        (out_dir / "logs").mkdir()
        (out_dir / "logs" / "busco.log").write_text("...")
        (out_dir / "hmmer_output").mkdir()
        (out_dir / "run_viridiplantae_odb12.json").write_text("{}")

        C._cleanup_busco_run(out_dir, keep=keep)

        assert list(out_dir.iterdir()) == [keep], list(out_dir.iterdir())
    finally:
        shutil.rmtree(out_dir)


def test_module7_matrix(tmp_matrix, tmp_stats):
    matrix = C.run_module7(OF_RESULTS, None, tmp_matrix, tmp_stats, min_species=1, force=True)
    assert matrix.shape == (6, 12), matrix.shape
    assert (matrix.values >= 0).all()
    return matrix


def test_module8_annotation(matrix, code5_to_go, tmp_annot, tmp_bridge):
    annot = C.run_module8(list(matrix.index), code5_to_go, OF_RESULTS / "Orthogroups" / "Orthogroups.tsv",
                           tmp_annot, tmp_bridge, threshold=0.5, force=True)
    assert len(annot) == 12
    assert annot["GO_consensus"].str.len().gt(0).any(), "expected at least one OG with GO consensus"

    # the bridge file must stay parseable by PCA/scripts/general_pca_common.load_go_ic_and_descriptions
    sys.path.insert(0, str(TEST_DIR.parent.parent / "PCA" / "scripts"))
    from general_pca_common import load_go_ic_and_descriptions
    ic, desc = load_go_ic_and_descriptions(tmp_bridge)
    assert len(desc) == 12
    assert all(og in desc for og in annot["Orthogroup"])


def main():
    subprocess.run([sys.executable, str(TEST_DIR / "make_test_data.py")], check=True)

    test_species_code5()
    test_module1_inventory()
    test_module1_flat_dir()
    test_module1_flat_dir_kingdom_filter()
    test_clean_and_prefix_dedup()
    test_run_busco_nonfatal_failure()
    test_busco_cleanup()

    tmp_matrix = TEST_DIR / "_tmp_mod07_matrix.tsv"
    tmp_stats = TEST_DIR / "_tmp_mod07_stats.tsv"
    matrix = test_module7_matrix(tmp_matrix, tmp_stats)

    manifest = C.collect_proteome_files(ROOT, non_interactive=True, taxid_map={})
    code_by_species_prefix = {"Arabidopsis_thaliana": "Atha", "Solanum_lycopersicum": "Slyco",
                               "Oryza_sativa": "Osati", "Physcomitrium_patens": "Ppate",
                               "Chlamydomonas_reinhardtii": "Crein", "Amborella_trichopoda": "Atric"}
    species_to_go = dict(zip(manifest["Species"], manifest["GOFile"]))
    code5_to_go = {code: species_to_go[species] for species, code in code_by_species_prefix.items()}

    tmp_annot = TEST_DIR / "_tmp_mod08_annot.tsv"
    tmp_bridge = TEST_DIR / "_tmp_mod08_bridge.tsv"
    test_module8_annotation(matrix, code5_to_go, tmp_annot, tmp_bridge)

    for f in (tmp_matrix, tmp_stats, tmp_annot, tmp_bridge):
        f.unlink(missing_ok=True)

    print("OK — all SenideOG self-checks passed")


if __name__ == "__main__":
    main()
