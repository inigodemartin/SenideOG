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


def test_busco_parses_real_summary_format():
    """The results line's bracket punctuation moved between BUSCO versions:
    v5 keeps F/M/n inside the brackets ('...,M:<pct>%,n:<count>]'), v6 closes
    the bracket after D and puts F/M/n after it ('...,D:<pct>%],F:<pct>%,M:
    <pct>%,n:<count>'). Parsing must not depend on where the brackets sit,
    or every species ends up NO_BUSCO_RESULT despite a clean BUSCO run."""
    formats = {
        "v5": "C:98.4%[S:97.2%,D:1.2%,F:0.7%,M:0.9%,n:425]\n",
        "v6": "C:98.1%[S:57.3%,D:40.8%],F:1.0%,M:1.0%,n:822\n",
    }
    workdir = TEST_DIR / "_tmp_busco_parse"
    orig_require_tool = C._require_tool
    try:
        C._require_tool = lambda name: name
        for tag, line in formats.items():
            out_dir = workdir / "busco" / f"Test3{tag}"
            out_dir.mkdir(parents=True, exist_ok=True)
            summary = out_dir / f"short_summary.specific.viridiplantae_odb12.Test3{tag}.txt"
            summary.write_text(line)
            result = C.run_busco(TEST_DIR / "nonexistent.fa", "viridiplantae_odb12", workdir, 4, f"Test3{tag}")
            assert result["C"] is not None, (tag, result)
    finally:
        C._require_tool = orig_require_tool
        shutil.rmtree(workdir, ignore_errors=True)


def test_busco_reruns_truncated_summary():
    """A summary file that exists but has no results line (BUSCO killed
    mid-run, e.g. OOM under --busco-jobs) must not be trusted as a
    checkpoint hit -- it has to be rerun, not permanently recorded as
    NO_BUSCO_RESULT."""
    workdir = TEST_DIR / "_tmp_busco_truncated"
    out_dir = workdir / "busco" / "Test6"
    out_dir.mkdir(parents=True, exist_ok=True)
    summary = out_dir / "short_summary.specific.viridiplantae_odb12.Test6.txt"
    summary.write_text("# BUSCO version is: 6.0.0\n# The lineage dataset is: viridiplantae_odb12\n")
    orig_require_tool, orig_run = C._require_tool, C._run
    try:
        C._require_tool = lambda name: name
        C._run = lambda cmd, **k: (
            summary.write_text("C:91.0%[S:88.0%,D:3.0%],F:5.0%,M:4.0%,n:400\n"),
            subprocess.CompletedProcess(cmd, 0))[1]
        result = C.run_busco(TEST_DIR / "nonexistent.fa", "viridiplantae_odb12", workdir, 4, "Test6")
        assert result["C"] == 91.0, result
    finally:
        C._require_tool, C._run = orig_require_tool, orig_run
        shutil.rmtree(workdir, ignore_errors=True)


def test_module3_tracks_unresolved_species():
    """Species that never made it past M1 (no GCA, no proteome file, ...)
    must still show up in mod03_qc.tsv -- otherwise the QC table can't be
    used to monitor every initially-selected species, only the survivors."""
    workdir = TEST_DIR / "_tmp_mod03_workdir"
    clean_dir = workdir / "proteomes_clean"
    clean_dir.mkdir(parents=True, exist_ok=True)
    results_path = TEST_DIR / "_tmp_mod03_qc.tsv"
    try:
        (clean_dir / "Test4.fa").write_text(">Test4|gene1\nMAAAA\n")
        manifest_df = pd.DataFrame([
            {"Species": "Test4_resolved", "GOFile": None, "Status": "resolved"},
            {"Species": "Test5_dropped", "GOFile": None, "Status": "no_gca_found"},
        ])
        proteome_stats = pd.DataFrame([
            {"Species": "Test4_resolved", "Code5": "Test4", "N_proteins": 1, "Source": "asis"},
        ])
        df = C.run_module3(proteome_stats, manifest_df, clean_dir, workdir, results_path,
                            lineage="viridiplantae_odb12", busco_pass=0.9, busco_flag=0.5, id_threshold=0.0,
                            threads=1, skip_busco=True, force=True)
        assert len(df) == 2, df
        resolved = df[df["Species"] == "Test4_resolved"].iloc[0]
        assert resolved["Status"] == "resolved" and resolved["BUSCO_status"] == "SKIPPED", resolved
        dropped = df[df["Species"] == "Test5_dropped"].iloc[0]
        assert dropped["Status"] == "no_gca_found" and dropped["BUSCO_status"] == "N/A", dropped
        assert pd.isna(dropped["N_proteins"]), dropped
    finally:
        shutil.rmtree(workdir, ignore_errors=True)
        results_path.unlink(missing_ok=True)


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
    test_busco_parses_real_summary_format()
    test_busco_reruns_truncated_summary()
    test_module3_tracks_unresolved_species()

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
