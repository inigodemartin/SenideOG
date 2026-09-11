#!/usr/bin/env python3
"""
SenideOG ("senide" = pariente, euskera + OG) -- builds a species x orthogroup
matrix for ~500 Viridiplantae proteomes, replacing the GO-term unit (biased
by embedding-similarity transfer, see go_phylostrata.py) with orthogroups
inferred from the sequences themselves.

M1 Inventory      -- scan <root>/<Species>/00_GenomeSource/<GCA>/ for proteomes + GO tables,
                      or a flat --fasta-dir of one proteome FASTA per species (no GO tables)
M2 Normalization   -- one isoform/gene, seqkit cleaning, <Code5>| ID prefix
M3 QC              -- BUSCO completeness, FASTA<->GO id match rate, disk check
M4 Core selection  -- ~64 species, stratified by order/family (APG IV)
M5 De novo         -- OrthoFinder on the core (diamond_ultra_sens, MSA trees)
M6 Assignment      -- OrthoFinder --assign for the remaining species (SHOOT)
M7 OG matrix       -- species x orthogroup counts, same shape as the GO matrix
M8 Annotation      -- GO consensus per orthogroup + IC/description bridge file
                       for PCA/scripts/general_pca_abundance.py --ic-file

Phylogeny/CAFE5 (fase 2) is a separate script, senideog_phylo.py, written
once M5/M6 outputs exist to build on -- see the design doc.
"""

VERSION = "v0.7.0"

import argparse
import getpass
import json
import os
import platform
import resource
import sys
import time
from datetime import datetime
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
import senideog_common as C  # noqa: E402


def _dated_log_path(logs_dir: Path, base_name: str) -> Path:
    date_str = datetime.now().strftime("%Y%m%d")
    candidate = logs_dir / f"{base_name}_{date_str}.log"
    if not candidate.exists():
        return candidate
    n = 2
    while (logs_dir / f"{base_name}_{date_str}_{n}.log").exists():
        n += 1
    return logs_dir / f"{base_name}_{date_str}_{n}.log"


def _validate_inputs(pairs: list) -> None:
    ok = True
    for flag, path in pairs:
        if path is not None and not path.exists():
            print(f"ERROR: {flag} not found: {path}", file=sys.stderr)
            ok = False
    if not ok:
        sys.exit(1)


def parse_args():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--root", type=Path, default=None,
                     help="FANTASIA_project root, e.g. /data/users/demartini/FANTASIA_project "
                          "(required unless --fasta-dir or --skip_mod1)")
    ap.add_argument("--fasta-dir", type=Path, default=None,
                     help="Folder with one proteome FASTA per species, flat (no 00_GenomeSource/<GCA>/ tree) "
                          "-- alternative to --root. No GFF/GO tables in this mode, so M8 has nothing to "
                          "annotate with (pass --skip_mod8)")
    ap.add_argument("--output", type=Path, required=True, help="Run output directory")
    ap.add_argument("--threads", type=int, default=8, help="Threads for seqkit/BUSCO/OrthoFinder (default: 8)")
    ap.add_argument("--species-taxid", type=Path, default=C.DEFAULT_SPECIES_TAXID_PATH,
                     help="Species -> TaxID lookup TSV (default: data/species_taxid.tsv, bundled in this repo)")
    ap.add_argument("--non_interactive", action="store_true",
                     help="Skip species with >1 genome assembly (GCA) instead of prompting")
    ap.add_argument("--kingdom", default="Viridiplantae",
                     help="With --fasta-dir, keep only species whose 'kingdom' column in --lineage-file matches "
                          "this value (default: Viridiplantae). Pass '' to disable the filter.")

    ap.add_argument("--min-protein-len", type=int, default=30,
                     help="Drop proteins shorter than this many aa in M2 (default: 30)")

    ap.add_argument("--busco-lineage", default="viridiplantae_odb12", help="BUSCO lineage dataset (default: viridiplantae_odb12)")
    ap.add_argument("--busco-jobs", type=int, default=1,
                     help="Number of BUSCO runs to execute in parallel; each still uses --threads internally, "
                          "so total CPU used ~= --busco-jobs x --threads (default: 1, sequential)")
    ap.add_argument("--busco-c-pass", type=float, default=0.85, help="BUSCO C%% >= this -> PASS (default: 0.85)")
    ap.add_argument("--busco-c-flag", type=float, default=0.80, help="BUSCO C%% >= this -> FLAG, below -> FAIL (default: 0.80)")
    ap.add_argument("--id-match-threshold", type=float, default=0.5,
                     help="Minimum mean FASTA<->GO-table id match rate; M3 aborts below this (default: 0.5)")
    ap.add_argument("--disk-estimate-gb", type=float, default=300.0,
                     help="Estimated disk needed for the core all-vs-all + gene trees (default: 300 GB)")

    ap.add_argument("--core-size", type=int, default=64, help="Target core size (default: 64)")
    ap.add_argument("--lineage-file", type=Path, default=C.DEFAULT_LINEAGE_PATH,
                     help="Species -> order/family (APG IV) TSV for stratified core selection "
                          "(default: data/species_lineage.tsv, bundled in this repo)")
    ap.add_argument("--core-species-override", type=Path, default=None,
                     help="TSV with a 'Species' column to use as the core directly, bypassing stratified "
                          "selection (manual curation, or local testing without real BUSCO/lineage data)")

    ap.add_argument("--min-species", type=int, default=4,
                     help="Drop orthogroups present in fewer than this many species from the M7 matrix (default: 4)")
    ap.add_argument("--go-consensus-threshold", type=float, default=0.5,
                     help="Fraction of an orthogroup's annotated members that must carry a GO for it to be "
                          "reported as consensus (default: 0.5)")

    for i in range(1, 9):
        ap.add_argument(f"--skip_mod{i}", action="store_true", help=f"Skip Module {i}")
    ap.add_argument("--force", action="store_true",
                     help="Rerun all steps from scratch even if intermediate outputs exist in workdir/")
    ap.add_argument("--dry_run", action="store_true",
                     help="Validate inputs and print the steps that would run, then exit without executing anything")
    ap.add_argument("--disable_co2_tracking", action="store_true",
                     help="Disable carbon footprint tracking even if codecarbon is installed")
    ap.add_argument("--version", action="version", version=f"%(prog)s {VERSION}")
    return ap.parse_args()


def main():
    args = parse_args()

    args.output = args.output.resolve()
    if args.root:
        args.root = args.root.resolve()
    if args.fasta_dir:
        args.fasta_dir = args.fasta_dir.resolve()
    args.species_taxid = args.species_taxid.resolve()
    args.lineage_file = args.lineage_file.resolve()
    if args.core_species_override:
        args.core_species_override = args.core_species_override.resolve()

    if args.root is not None and args.fasta_dir is not None:
        print("ERROR: pass only one of --root / --fasta-dir", file=sys.stderr)
        sys.exit(1)
    if not args.skip_mod1 and args.root is None and args.fasta_dir is None:
        print("ERROR: --root or --fasta-dir is required (unless --skip_mod1)", file=sys.stderr)
        sys.exit(1)
    _validate_inputs([("--root", args.root), ("--fasta-dir", args.fasta_dir)])

    results = args.output / "results"
    workdir = args.output / "workdir"
    logs_dir = args.output / "logs"
    for d in (results, workdir, logs_dir):
        d.mkdir(parents=True, exist_ok=True)

    C._LOG_FH = open(_dated_log_path(logs_dir, "Run_SenideOG"), "a")
    sep = "=" * 62
    C._LOG_FH.write(f"{sep}\n  SenideOG {VERSION}  —  Run Log\n{sep}\n")
    C._LOG_FH.write(f"Date      : {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
    C._LOG_FH.write(f"User      : {getpass.getuser()}\n")
    C._LOG_FH.write(f"Server    : {platform.node()}\n")
    C._LOG_FH.write(f"OS        : {platform.system()} {platform.release()} ({platform.machine()})\n")
    C._LOG_FH.write(f"Directory : {os.getcwd()}\n")
    C._LOG_FH.write(f"Command   : {' '.join(sys.argv)}\n")
    C._LOG_FH.write(f"{sep}\n\n")
    C._LOG_FH.flush()

    if args.dry_run:
        C._banner("Dry run — no steps will be executed")
        if args.fasta_dir:
            C._log(f"  Fasta dir   : {args.fasta_dir}")
        else:
            C._log(f"  Root        : {args.root}")
        C._log(f"  Output      : {args.output}/")
        C._log("  Steps that would run:")
        labels = [
            "M1 Inventory      → results/mod01_manifest.tsv",
            "M2 Normalization  → workdir/proteomes_clean/, results/mod02_proteome_stats.tsv",
            "M3 QC (BUSCO + id match + disk) → results/mod03_qc.tsv",
            "M4 Core selection → results/mod04_core_species.tsv",
            "M5 De novo inference on core (OrthoFinder) → workdir/of_core/",
            "M6 Assignment of the rest (OrthoFinder --assign) → workdir/of_assign/",
            "M7 Species x orthogroup matrix → results/mod07_og_matrix.tsv",
            "M8 OG annotation + IC bridge file → results/mod08_og_annotation.tsv",
        ]
        for i, label in enumerate(labels, 1):
            if not getattr(args, f"skip_mod{i}"):
                C._log(f"    [{i}] {label}")
        C._log("  Exiting (--dry_run).")
        sys.exit(0)

    if args.force:
        C._log("--force set: all steps will rerun regardless of existing outputs")
    elif workdir.exists() and any(workdir.iterdir()):
        C._log("Existing workdir found — resuming from checkpoints (use --force to rerun all steps from scratch)")

    _tracker = None
    if args.disable_co2_tracking:
        C._log("  Carbon footprint tracking disabled (--disable_co2_tracking)")
    else:
        try:
            from codecarbon import EmissionsTracker
            _tracker = EmissionsTracker(output_dir=str(logs_dir), output_file=f"{args.output.name}.emissions.csv",
                                        project_name="SenideOG", log_level="warning")
            _tracker.start()
            C._log("  codecarbon tracker started")
        except ImportError:
            C._log("  codecarbon not installed — carbon tracking skipped (conda install -c conda-forge codecarbon)")

    t_start = time.monotonic()
    stats_json = {}

    manifest_path = results / "mod01_manifest.tsv"
    C._banner("Module 1 — Inventory and collection")
    if args.skip_mod1:
        if not manifest_path.exists():
            print(f"ERROR: --skip_mod1 set but {manifest_path} doesn't exist yet — run Module 1 at least once first", file=sys.stderr)
            sys.exit(1)
        C._log("  --skip_mod1 set: loading existing manifest")
        import pandas as pd
        manifest_df = pd.read_csv(manifest_path, sep="\t")
    else:
        manifest_df = C.run_module1(args.root, args.fasta_dir, manifest_path, args.species_taxid, args.non_interactive,
                                     args.force, args.lineage_file, args.kingdom or None)
    stats_json["n_species_manifest"] = len(manifest_df)
    stats_json["n_species_resolved"] = int((manifest_df["Status"] == "resolved").sum())

    proteome_stats_path = results / "mod02_proteome_stats.tsv"
    clean_dir = workdir / "proteomes_clean"
    if not args.skip_mod2:
        C._banner("Module 2 — Normalization")
        proteome_stats = C.run_module2(manifest_df, workdir, proteome_stats_path, args.min_protein_len, args.force)
    elif proteome_stats_path.exists():
        import pandas as pd
        proteome_stats = pd.read_csv(proteome_stats_path, sep="\t")
    else:
        proteome_stats = None

    qc_path = results / "mod03_qc.tsv"
    if not args.skip_mod3:
        C._banner("Module 3 — QC (BUSCO, id match, disk)")
        C.check_disk(workdir, args.disk_estimate_gb)
        qc_df = C.run_module3(proteome_stats, manifest_df, clean_dir, workdir, qc_path,
                               args.busco_lineage, args.busco_c_pass, args.busco_c_flag,
                               args.id_match_threshold, args.threads, skip_busco=False, force=args.force,
                               busco_jobs=args.busco_jobs)
    elif qc_path.exists():
        import pandas as pd
        qc_df = pd.read_csv(qc_path, sep="\t")
    else:
        qc_df = None

    core_path = results / "mod04_core_species.tsv"
    if not args.skip_mod4:
        C._banner("Module 4 — Core selection")
        core_df = C.run_module4(qc_df, args.lineage_file, core_path, args.core_size,
                                 args.core_species_override, args.force)
    elif core_path.exists():
        import pandas as pd
        core_df = pd.read_csv(core_path, sep="\t")
    else:
        core_df = None
    if core_df is not None:
        stats_json["n_core_species"] = len(core_df)

    of_core_dir = workdir / "of_core"
    if not args.skip_mod5:
        C._banner("Module 5 — De novo inference on the core")
        core_codes = core_df["Code5"].tolist()
        core_results = C.run_module5(core_codes, clean_dir, of_core_dir, args.threads, args.force)
    else:
        core_results = C.find_orthofinder_results(of_core_dir)

    of_assign_dir = workdir / "of_assign"
    if not args.skip_mod6:
        C._banner("Module 6 — Assignment of the remaining species")
        core_species = set(core_df["Species"])
        rest = proteome_stats[~proteome_stats["Species"].isin(core_species)]
        rest_codes = rest["Code5"].tolist()
        if rest_codes:
            assign_results = C.run_module6(rest_codes, clean_dir, core_results, of_assign_dir, args.threads, args.force)
        else:
            C._log("  no non-core species to assign")
            assign_results = None
    else:
        assign_results = C.find_orthofinder_results(of_assign_dir)

    matrix_path = results / "mod07_og_matrix.tsv"
    stats_path = results / "mod07_species_stats.tsv"
    if not args.skip_mod7:
        C._banner("Module 7 — Species x orthogroup matrix")
        matrix = C.run_module7(core_results, assign_results, matrix_path, stats_path, args.min_species, args.force)
    elif matrix_path.exists():
        import pandas as pd
        matrix = pd.read_csv(matrix_path, sep="\t", index_col=0)
    else:
        matrix = None
    if matrix is not None:
        stats_json["og_matrix_shape"] = list(matrix.shape)

    if not args.skip_mod8:
        C._banner("Module 8 — Orthogroup annotation")
        og_source = assign_results or core_results
        orthogroups_tsv = og_source / "Orthogroups" / "Orthogroups.tsv"
        code_to_species = dict(zip(proteome_stats["Code5"], proteome_stats["Species"]))
        species_to_go = dict(zip(manifest_df["Species"], manifest_df["GOFile"]))
        code5_to_go = {code: species_to_go.get(species) for code, species in code_to_species.items()}
        C.run_module8(list(matrix.index), code5_to_go, orthogroups_tsv,
                       results / "mod08_og_annotation.tsv", results / "mod08_og_desc_ic.tsv",
                       args.go_consensus_threshold, args.force)

    elapsed_s = time.monotonic() - t_start
    ru = resource.getrusage(resource.RUSAGE_SELF)
    peak_mem_mb = (ru.ru_maxrss / (1024 * 1024) if platform.system() == "Darwin" else ru.ru_maxrss / 1024)

    emissions_kg = None
    if _tracker is not None:
        try:
            emissions_kg = _tracker.stop()
        except Exception:
            pass

    summary = {
        "date": datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
        "version": VERSION,
        "root": str(args.root) if args.root else None,
        "fasta_dir": str(args.fasta_dir) if args.fasta_dir else None,
        **stats_json,
        "parameters": {
            "core_size": args.core_size,
            "min_species": args.min_species,
            "go_consensus_threshold": args.go_consensus_threshold,
            "busco_lineage": args.busco_lineage,
        },
        "resource_usage": {
            "wall_clock_s": round(elapsed_s, 1),
            "peak_mem_mb": round(peak_mem_mb, 1),
            "emissions_kg_CO2eq": emissions_kg,
        },
    }
    with open(results / f"{args.output.name}.run_summary.json", "w") as fh:
        json.dump(summary, fh, indent=2)
        fh.write("\n")

    C._banner("Done")
    C._log(f"  Wall-clock time : {elapsed_s:.1f}s")
    C._log(f"  Peak memory     : {peak_mem_mb:.1f} MB")

    if C._LOG_FH is not None:
        C._LOG_FH.close()


if __name__ == "__main__":
    main()
