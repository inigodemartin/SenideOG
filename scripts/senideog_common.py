#!/usr/bin/env python3
"""
Shared helpers for senideog.py (M1-M8): logging/subprocess infra (same shape
as every other tool in this workspace, kept here once instead of duplicated
between senideog.py and the future senideog_phylo.py), proteome discovery,
FASTA cleaning, QC, core-species selection, and OrthoFinder-output parsing
into the species x orthogroup matrix and the OG annotation bridge file.
"""

import concurrent.futures
import fnmatch
import gzip
import re
import shutil
import subprocess
import sys
from collections import Counter
from datetime import datetime
from pathlib import Path

import pandas as pd

DEFAULT_LINEAGE_PATH = Path(__file__).parent.parent / "data" / "species_lineage.tsv"
DEFAULT_SPECIES_TAXID_PATH = Path(__file__).parent.parent / "data" / "species_taxid.tsv"

# ------------------------------------------------------------------ logging
_LOG_FH = None


def _log(msg: str) -> None:
    ts = datetime.now().strftime("%H:%M:%S")
    line = f"[{ts}] {msg}"
    print(line, file=sys.stderr)
    if _LOG_FH is not None:
        print(line, file=_LOG_FH, flush=True)


def _banner(title: str) -> None:
    bar = "─" * (len(title) + 4)
    _log(f"┌{bar}┐")
    _log(f"│  {title}  │")
    _log(f"└{bar}┘")


def _checkpoint(path: Path, label: str, force: bool) -> bool:
    if not force and path.exists() and path.stat().st_size > 0:
        _log(f"  [checkpoint] {label} — {path.name} already exists, skipping")
        return True
    return False


def _require_tool(name: str) -> str:
    tool = shutil.which(name)
    if tool is None:
        print(f"ERROR: '{name}' not found in PATH.\n"
              f"       Install with:  conda install -c bioconda {name}",
              file=sys.stderr)
        sys.exit(1)
    return tool


def _run(cmd: list, capture_stdout: bool = False, env: dict = None, cwd: Path = None,
         check: bool = True) -> subprocess.CompletedProcess:
    """check=False lets the caller handle a non-zero exit itself (used for
    BUSCO: one species' bad input shouldn't sys.exit() the whole run out
    from under the other species still processing in parallel)."""
    _log(f"  $ {' '.join(str(c) for c in cmd)}")
    result = subprocess.run(cmd,
                             stdout=subprocess.PIPE if capture_stdout else None,
                             stderr=subprocess.PIPE, text=True,
                             env=env, cwd=cwd)
    if check and result.returncode != 0:
        print(f"ERROR: command failed (exit {result.returncode}):\n"
              f"{result.stderr[-3000:]}", file=sys.stderr)
        sys.exit(1)
    return result


# --------------------------------------------------------------------- FASTA
def iter_fasta(path: Path):
    """Yield (id, seq) for every record. id = header up to first whitespace."""
    opener = gzip.open if str(path).endswith(".gz") else open
    header, chunks = None, []
    with opener(path, "rt") as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(chunks)
                header = line[1:].split()[0]
                chunks = []
            else:
                chunks.append(line)
        if header is not None:
            yield header, "".join(chunks)


def fasta_ids(path: Path) -> set:
    return {h for h, _ in iter_fasta(path)}


def count_fasta(path: Path) -> int:
    return sum(1 for _ in iter_fasta(path))


def mean_protein_length(path: Path) -> float:
    lengths = [len(seq) for _, seq in iter_fasta(path)]
    return sum(lengths) / len(lengths) if lengths else float("nan")


# ------------------------------------------------------- M1: inventory / collection
# Fixed preference order (per species, recorded as "Source"): a species'
# proteome file is picked from the first tier that has a match -- see the
# design doc's rationale (CD-HIT100 collapses recent real gene duplicates,
# so it's used only when nothing better is on disk).
PROTEOME_TIERS = [
    ("longest_isof", ["*_genomic_longest_isof_proteins.fasta"]),
    ("cdhit100", ["*_uniq_cdhit100_5k_removed.pep", "*_uniq_cdhit100.pep"]),
    ("rmdup", ["*_rmdup.fasta"]),
]

FANTASIA_DIR_RE = re.compile(r"^FANTASIA_2025(_.+)?$")
GOS_MERGED_RE = re.compile(r"^(\w+)_GOs_merged\.tsv$")


def find_proteome(gca_dir: Path):
    """(path, source_tag) for the first matching tier, or (None, None)."""
    for source, patterns in PROTEOME_TIERS:
        for pattern in patterns:
            hits = sorted(gca_dir.glob(pattern))
            if hits:
                return hits[0], source
    return None, None


def find_gff(gca_dir: Path):
    for pattern in ("*.gff3", "*.gff"):
        hits = sorted(gca_dir.glob(pattern))
        if hits:
            return hits[0]
    return None


def find_go_file(species_dir: Path):
    """First `<prefix>_GOs_merged.tsv` under 04_FunctionalAnnotation/FANTASIA_2025*/.
    Ambiguity (multiple valid annotation runs, e.g. different gene models) is
    only resolved by M8 at consensus time -- here we just pick the first
    match and log so it's visible in the manifest, not silently arbitrary."""
    func_dir = species_dir / "04_FunctionalAnnotation"
    if not func_dir.is_dir():
        return None
    hits = []
    for candidate in sorted(func_dir.iterdir()):
        if candidate.is_dir() and FANTASIA_DIR_RE.match(candidate.name):
            hits += sorted(candidate.glob("*_GOs_merged.tsv"))
    if not hits:
        return None
    if len(hits) > 1:
        _log(f"  [WARN] {species_dir.name}: {len(hits)} GO annotation files found, using {hits[0].name}")
    return hits[0]


def collect_gcas(species_dir: Path) -> list:
    genome_dir = species_dir / "00_GenomeSource"
    if not genome_dir.is_dir():
        return []
    return sorted(p for p in genome_dir.iterdir() if p.is_dir())


def prompt_choice_gca(species: str, gcas: list):
    """Interactively ask which GCA assembly to use. Returns a Path or None to skip."""
    print(f"\n[ambiguous] {species} has {len(gcas)} genome assemblies (GCA):")
    for i, gca in enumerate(gcas, 1):
        print(f"  {i}. {gca.name}")
    print("  0. skip this species")
    while True:
        choice = input(f"  Choose 1-{len(gcas)} (0 = skip): ").strip()
        if choice == "0":
            return None
        if choice.isdigit() and 1 <= int(choice) <= len(gcas):
            return gcas[int(choice) - 1]
        print("  invalid choice, try again")


def collect_proteome_files(root: Path, non_interactive: bool, taxid_map: dict) -> pd.DataFrame:
    """One row per species under <root>: Species, TaxID, Accession, ProteomeFile, Source, GFF, GOFile."""
    rows = []
    n_no_gca, n_no_proteome, n_multi = 0, 0, 0
    for species_dir in sorted(p for p in root.iterdir() if p.is_dir()):
        species = species_dir.name
        gcas = collect_gcas(species_dir)
        if not gcas:
            n_no_gca += 1
            rows.append({"Species": species, "TaxID": taxid_map.get(species), "Accession": None,
                         "ProteomeFile": None, "Source": None, "GFF": None,
                         "GOFile": find_go_file(species_dir), "Status": "no_gca_found"})
            continue

        gca_dir = gcas[0]
        if len(gcas) > 1:
            n_multi += 1
            if non_interactive:
                rows.append({"Species": species, "TaxID": taxid_map.get(species), "Accession": None,
                             "ProteomeFile": None, "Source": None, "GFF": None,
                             "GOFile": find_go_file(species_dir), "Status": "ambiguous_gca_skipped"})
                continue
            chosen = prompt_choice_gca(species, gcas)
            if chosen is None:
                rows.append({"Species": species, "TaxID": taxid_map.get(species), "Accession": None,
                             "ProteomeFile": None, "Source": None, "GFF": None,
                             "GOFile": find_go_file(species_dir), "Status": "skipped_by_user"})
                continue
            gca_dir = chosen

        proteome, source = find_proteome(gca_dir)
        go_file = find_go_file(species_dir)
        if proteome is None:
            n_no_proteome += 1
            rows.append({"Species": species, "TaxID": taxid_map.get(species), "Accession": gca_dir.name,
                         "ProteomeFile": None, "Source": None, "GFF": find_gff(gca_dir), "GOFile": go_file,
                         "Status": "no_proteome_found"})
            continue

        rows.append({"Species": species, "TaxID": taxid_map.get(species), "Accession": gca_dir.name,
                     "ProteomeFile": str(proteome), "Source": source, "GFF": find_gff(gca_dir),
                     "GOFile": str(go_file) if go_file else None, "Status": "resolved"})

    _log(f"  {root}: {len(rows)} species scanned — "
         f"{n_multi} with >1 GCA, {n_no_gca} with no GCA, {n_no_proteome} with no matching proteome file")
    cols = ["Species", "TaxID", "Accession", "ProteomeFile", "Source", "GFF", "GOFile", "Status"]
    return pd.DataFrame(rows, columns=cols)


FASTA_EXTS = (".fa", ".fasta", ".pep", ".fa.gz", ".fasta.gz", ".pep.gz")


def load_kingdom_species(lineage_path: Path, kingdom: str) -> set:
    """Species names from lineage_path whose 'kingdom' column matches (e.g. Viridiplantae)."""
    lineage_df = pd.read_csv(lineage_path, sep="\t")
    return set(lineage_df.loc[lineage_df["kingdom"] == kingdom, "Species"])


def collect_proteome_files_from_dir(fasta_dir: Path, taxid_map: dict, kingdom_species: set = None) -> pd.DataFrame:
    """One row per FASTA file directly under fasta_dir (flat layout: no
    00_GenomeSource/<GCA>/ tree, no GO tables). Species = filename with its
    FASTA extension stripped. Source is tagged from the same PROTEOME_TIERS
    patterns the GCA-tree layout uses, so M2 still skips the AGAT isoform
    collapse when a name already says 'longest_isof'; falls back to 'asis'.

    kingdom_species, if given (from --lineage-file, filtered to --kingdom),
    restricts the folder -- which may hold non-target species, e.g. fungi --
    to just those species names."""
    files = sorted(p for p in fasta_dir.iterdir() if p.is_file() and p.name.lower().endswith(FASTA_EXTS))
    rows = []
    n_skipped = 0
    for f in files:
        species = f.name
        for ext in sorted(FASTA_EXTS, key=len, reverse=True):
            if species.lower().endswith(ext):
                species = species[: -len(ext)]
                break
        if kingdom_species is not None and species not in kingdom_species:
            n_skipped += 1
            continue
        source = next((tag for tag, patterns in PROTEOME_TIERS
                        if any(fnmatch.fnmatch(f.name, pat) for pat in patterns)), "asis")
        rows.append({"Species": species, "TaxID": taxid_map.get(species), "Accession": None,
                     "ProteomeFile": str(f), "Source": source, "GFF": None, "GOFile": None,
                     "Status": "resolved"})
    if kingdom_species is not None:
        _log(f"  {fasta_dir}: {len(files)} FASTA files found, {len(rows)} match the lineage-file kingdom filter "
             f"({n_skipped} skipped)")
        print(f"Selected {len(rows)} species", file=sys.stdout)
    else:
        _log(f"  {fasta_dir}: {len(rows)} FASTA files found")
    cols = ["Species", "TaxID", "Accession", "ProteomeFile", "Source", "GFF", "GOFile", "Status"]
    return pd.DataFrame(rows, columns=cols)


def run_module1(root: Path, fasta_dir: Path, manifest_path: Path, species_taxid: Path, non_interactive: bool,
                 force: bool, lineage_path: Path = None, kingdom: str = None) -> pd.DataFrame:
    if _checkpoint(manifest_path, "manifest", force):
        return pd.read_csv(manifest_path, sep="\t")
    taxid_map = {}
    if species_taxid and species_taxid.exists():
        taxid_df = pd.read_csv(species_taxid, sep="\t")
        taxid_map = dict(zip(taxid_df["Species"], taxid_df["TaxID"]))
    else:
        _log(f"  [WARN] --species-taxid not found ({species_taxid}) — TaxID column left empty")

    if fasta_dir:
        kingdom_species = None
        if kingdom:
            if lineage_path and lineage_path.exists():
                kingdom_species = load_kingdom_species(lineage_path, kingdom)
                _log(f"  --kingdom={kingdom}: {len(kingdom_species)} species in {lineage_path.name}")
            else:
                _log(f"  [WARN] --lineage-file not found ({lineage_path}) — --kingdom={kingdom} filter skipped")
        df = collect_proteome_files_from_dir(fasta_dir, taxid_map, kingdom_species)
    else:
        df = collect_proteome_files(root, non_interactive, taxid_map)
    df.to_csv(manifest_path, sep="\t", index=False)
    _log(f"  wrote {manifest_path.name} ({len(df)} species, {(df['Status'] == 'resolved').sum()} resolved)")
    return df


# ------------------------------------------------------------- M2: normalization
def species_code5(species_name: str, taken: set) -> str:
    """<Genus>_<epithet> -> 5-char OrthoFinder-style code, e.g. Solanum_lycopersicum -> Slyco."""
    genus, _, epithet = species_name.partition("_")
    base = ((genus[:1] + epithet[:4]) or species_name[:5]).capitalize()
    code = base
    n = 1
    while code in taken:
        n += 1
        suffix = str(n)
        code = (base[: max(1, 5 - len(suffix))] + suffix)[:5]
    taken.add(code)
    return code


def gff_longest_isoform_ids(filtered_gff: Path) -> set:
    """mRNA/protein IDs surviving `agat_sp_keep_longest_isoform.pl` (parses GFF3 attribute ID=...)."""
    ids = set()
    id_re = re.compile(r"ID=([^;]+)")
    with open(filtered_gff) as fh:
        for line in fh:
            if line.startswith("#") or "\t" not in line:
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9 or fields[2].lower() not in ("mrna", "cds"):
                continue
            m = id_re.search(fields[8])
            if m:
                ids.add(m.group(1))
    return ids


def keep_longest_isoform(fasta_in: Path, gff_in: Path, workdir: Path, code5: str) -> Path:
    agat = _require_tool("agat_sp_keep_longest_isoform.pl")
    filtered_gff = workdir / f"{code5}.longest_isof.gff"
    if not _checkpoint(filtered_gff, f"{code5} longest-isoform GFF", False):
        _run([agat, "--gff", str(gff_in), "-o", str(filtered_gff)], cwd=workdir)
    keep_ids = gff_longest_isoform_ids(filtered_gff)
    n_matched = sum(1 for h, _ in iter_fasta(fasta_in) if any(h.startswith(k) or k.startswith(h) for k in keep_ids))
    if keep_ids and n_matched / max(1, count_fasta(fasta_in)) < 0.5:
        _log(f"  [WARN] {code5}: longest-isoform GFF IDs match only "
             f"{n_matched}/{count_fasta(fasta_in)} FASTA records — gene-model IDs may not correspond, check manually")
    filtered_fasta = workdir / f"{code5}.longest_isof.fa"
    with open(filtered_fasta, "w") as out:
        for h, seq in iter_fasta(fasta_in):
            if any(h.startswith(k) or k.startswith(h) for k in keep_ids):
                out.write(f">{h}\n{seq}\n")
    return filtered_fasta


def clean_and_prefix_fasta(fasta_in: Path, out_path: Path, code5: str, min_len: int) -> dict:
    """seqkit for the real cleaning (strip trailing '*', U/J/Z/B -> X, drop
    < min_len aa); prefixing with <Code5>| happens in the same pass since
    seqkit has no rename-with-prefix mode that also strips *. Duplicate
    headers in the raw proteome (seen in the wild, e.g. repeated organelle
    genes) are uniquified with a '_N' suffix -- BUSCO hard-errors on a
    duplicate id and OrthoFinder would silently conflate the two proteins.
    First occurrence keeps its id as-is; the 2nd/3rd/... get a '_2'/'_3'/...
    suffix, skipping any candidate that collides with an id that already
    exists in the raw proteome (e.g. 'gene1' repeated would naively become
    'gene1_2', but if the proteome separately already has a distinct
    'gene1_2' entry that just reintroduces the duplicate BUSCO chokes on)."""
    seqkit = _require_tool("seqkit")
    tmp = out_path.with_suffix(".tmp.fa")
    _run([seqkit, "seq", "-g", "-M", "999999", "-m", str(min_len), str(fasta_in), "-o", str(tmp)])
    n_kept = 0
    all_headers = {h for h, _ in iter_fasta(tmp)}
    seen = Counter()
    used = set()
    with open(tmp) as fin, open(out_path, "w") as fout:
        for h, seq in iter_fasta(Path(tmp)):
            seq = seq.rstrip("*").upper()
            seq = re.sub(r"[UJZB]", "X", seq)
            seen[h] += 1
            if seen[h] > 1:
                n = seen[h]
                candidate = f"{h}_{n}"
                while candidate in all_headers or candidate in used:
                    n += 1
                    candidate = f"{h}_{n}"
                h = candidate
            used.add(h)
            fout.write(f">{code5}|{h}\n{seq}\n")
            n_kept += 1
    tmp.unlink()
    n_dup = sum(c - 1 for c in seen.values() if c > 1)
    if n_dup:
        _log(f"  [WARN] {code5}: {n_dup} duplicate sequence ids in the raw proteome, uniquified with '_N'")
    return {"n_kept": n_kept}


def run_module2(manifest_df: pd.DataFrame, workdir: Path, results_path: Path, min_len: int, force: bool) -> pd.DataFrame:
    clean_dir = workdir / "proteomes_clean"
    clean_dir.mkdir(parents=True, exist_ok=True)
    if _checkpoint(results_path, "proteome stats", force):
        return pd.read_csv(results_path, sep="\t")

    resolved = manifest_df[manifest_df["Status"] == "resolved"]
    taken_codes = set()
    rows = []
    for _, row in resolved.iterrows():
        species, source = row["Species"], row["Source"]
        code5 = species_code5(species, taken_codes)
        raw_fasta = Path(row["ProteomeFile"])
        out_fasta = clean_dir / f"{code5}.fa"
        if _checkpoint(out_fasta, f"{species} ({code5})", force):
            rows.append({"Species": species, "Code5": code5, "N_proteins": count_fasta(out_fasta), "Source": source})
            continue

        fasta_for_cleaning = raw_fasta
        effective_source = source
        if source != "longest_isof":
            gff = Path(row["GFF"]) if pd.notna(row.get("GFF")) else None
            if gff and gff.exists():
                _log(f"  {species}: collapsing to one isoform/gene via AGAT ({gff.name})")
                fasta_for_cleaning = keep_longest_isoform(raw_fasta, gff, workdir, code5)
            else:
                _log(f"  {species}: no GFF available — keeping proteome as-is (Source stays '{source}')")

        stats = clean_and_prefix_fasta(fasta_for_cleaning, out_fasta, code5, min_len)
        rows.append({"Species": species, "Code5": code5, "N_proteins": stats["n_kept"], "Source": effective_source})

    df = pd.DataFrame(rows, columns=["Species", "Code5", "N_proteins", "Source"])
    df.to_csv(results_path, sep="\t", index=False)
    _log(f"  wrote {results_path.name} ({len(df)} species cleaned)")
    return df


# ----------------------------------------------------------------- M3: QC
def ensure_busco_dataset(lineage: str, workdir: Path) -> None:
    """Pre-download+extract the BUSCO lineage dataset once, sequentially,
    before any parallel `busco` runs start -- otherwise N concurrent BUSCO
    processes race to populate the same shared busco_downloads/ cache on the
    very first run. No-op (fast) if already downloaded."""
    busco = _require_tool("busco")
    busco_dir = workdir / "busco"
    busco_dir.mkdir(parents=True, exist_ok=True)
    _run([busco, "--download", lineage], cwd=busco_dir)


def _cleanup_busco_run(out_dir: Path, keep: Path) -> None:
    """Delete everything BUSCO wrote under out_dir except `keep` (the parsed
    short_summary file) -- hmmer_output/, busco_sequences/, logs/, etc. are
    large and fully reproducible by rerunning BUSCO, so there's no reason to
    keep them per species once the C/S/D/F/M numbers are parsed out."""
    for p in out_dir.iterdir():
        if p != keep:
            shutil.rmtree(p) if p.is_dir() else p.unlink()


def _busco_summary_valid(summary: Path) -> bool:
    """A summary file can exist (checkpoint-visible) yet be incomplete --
    BUSCO writes the header first and the '***** Results: *****' block
    last, so a process killed mid-run (OOM, timeout, scheduler preemption)
    leaves a non-empty file with no results line. Only trust it once the
    'C:<pct>%' token is actually present."""
    return (summary.exists() and summary.stat().st_size > 0
            and re.search(r"C:[\d.]+%", summary.read_text()) is not None)


def run_busco(fasta: Path, lineage: str, workdir: Path, threads: int, code5: str) -> dict:
    busco = _require_tool("busco")
    out_dir = workdir / "busco" / code5
    summary = out_dir / f"short_summary.specific.{lineage}.{code5}.txt"
    if _busco_summary_valid(summary):
        _log(f"  [checkpoint] {code5} BUSCO — {summary.name} already exists, skipping")
    else:
        if summary.exists():
            _log(f"  [WARN] {code5}: existing BUSCO summary has no results line "
                 f"(likely killed mid-run) -- rerunning")
        (workdir / "busco").mkdir(parents=True, exist_ok=True)
        result = _run([busco, "-i", str(fasta), "-m", "proteins", "-l", lineage, "-o", code5,
                       "-c", str(threads), "-f"], cwd=workdir / "busco", check=False)
        if result.returncode != 0:
            _log(f"  [WARN] {code5}: BUSCO failed (exit {result.returncode}) -- recorded as "
                 f"NO_BUSCO_RESULT, other species continue. stderr tail:\n{(result.stderr or '')[-1500:]}")
    if not _busco_summary_valid(summary):
        if out_dir.exists():
            shutil.rmtree(out_dir)
        return {"C": None, "S": None, "D": None, "F": None, "M": None}
    _cleanup_busco_run(out_dir, keep=summary)
    text = summary.read_text()
    # Field positions relative to the brackets moved between BUSCO versions
    # (v5: "C:_%[S:_%,D:_%,F:_%,M:_%,n:_]"; v6: "C:_%[S:_%,D:_%],F:_%,M:_%,n:_")
    # -- look up each "<letter>:<pct>%" token independently instead of
    # depending on where the brackets close.
    fields = {}
    for key in ("C", "S", "D", "F", "M"):
        m = re.search(rf"{key}:([\d.]+)%", text)
        fields[key] = float(m.group(1)) if m else None
    if fields["C"] is None:
        return {"C": None, "S": None, "D": None, "F": None, "M": None}
    return fields


def id_match_rate(fasta_path: Path, go_file: Path) -> float:
    """Fraction of GO-annotated protein IDs that also appear in the cleaned
    FASTA -- catches the #1 project risk (gene models don't correspond)."""
    prot_ids = fasta_ids(fasta_path)
    prot_ids_bare = {h.split("|", 1)[-1] for h in prot_ids}
    go_ids = set()
    with open(go_file) as fh:
        for line in fh:
            pid = line.split("\t", 1)[0].strip()
            if pid:
                go_ids.add(pid)
    if not go_ids:
        return float("nan")
    return len(go_ids & prot_ids_bare) / len(go_ids)


def run_module3(proteome_stats: pd.DataFrame, manifest_df: pd.DataFrame, clean_dir: Path, workdir: Path,
                results_path: Path, lineage: str, busco_pass: float, busco_flag: float, id_threshold: float,
                threads: int, skip_busco: bool, force: bool, busco_jobs: int = 1) -> pd.DataFrame:
    if _checkpoint(results_path, "QC table", force):
        return pd.read_csv(results_path, sep="\t")

    go_by_species = dict(zip(manifest_df["Species"], manifest_df["GOFile"]))
    species_rows = list(proteome_stats.iterrows())

    def qc_row(row):
        species, code5 = row["Species"], row["Code5"]
        fasta = clean_dir / f"{code5}.fa"
        busco = {"C": None} if skip_busco else run_busco(fasta, lineage, workdir, threads, code5)
        c = busco.get("C")
        if c is None:
            busco_status = "SKIPPED" if skip_busco else "NO_BUSCO_RESULT"
        elif c >= busco_pass * 100:
            busco_status = "PASS"
        elif c >= busco_flag * 100:
            busco_status = "FLAG"
        else:
            busco_status = "FAIL"
        go_file = go_by_species.get(species)
        match_rate = (id_match_rate(fasta, Path(go_file))
                      if go_file and pd.notna(go_file) and Path(go_file).exists() else float("nan"))
        return {"Species": species, "Code5": code5, "Status": "resolved", "N_proteins": row["N_proteins"],
                "Mean_protein_length": round(mean_protein_length(fasta), 1),
                "BUSCO_C": c, "BUSCO_status": busco_status, "GO_ID_match_rate": match_rate}

    if not skip_busco and busco_jobs > 1:
        n_cached = sum(1 for _, row in species_rows
                       if (workdir / "busco" / row["Code5"] /
                           f"short_summary.specific.{lineage}.{row['Code5']}.txt").exists())
        _log(f"  running BUSCO for {len(species_rows)} species ({n_cached} already done, reused as-is), "
             f"--busco-jobs={busco_jobs} in parallel (--threads={threads} each)")
        ensure_busco_dataset(lineage, workdir)
        with concurrent.futures.ThreadPoolExecutor(max_workers=busco_jobs) as pool:
            futures = {pool.submit(qc_row, row): i for i, (_, row) in enumerate(species_rows)}
            rows = [None] * len(species_rows)
            try:
                for fut in concurrent.futures.as_completed(futures):
                    rows[futures[fut]] = fut.result()
            except BaseException:
                # a species' BUSCO run failed (_run() exits) -- drop jobs not yet started instead of
                # burning through the rest of the list before the failure surfaces
                pool.shutdown(wait=False, cancel_futures=True)
                raise
    else:
        rows = [qc_row(row) for _, row in species_rows]

    df = pd.DataFrame(rows)

    # Species dropped before M2/M3 (no GCA, no proteome file, ambiguous GCA,
    # user-skipped) never reach qc_row -- keep them in the table anyway so
    # every initially-selected species can be tracked, not just the ones
    # that made it this far.
    not_resolved = manifest_df[~manifest_df["Species"].isin(df["Species"])]
    if len(not_resolved):
        extra = pd.DataFrame({"Species": not_resolved["Species"], "Code5": None,
                               "Status": not_resolved["Status"], "N_proteins": None,
                               "Mean_protein_length": None, "BUSCO_C": None,
                               "BUSCO_status": "N/A", "GO_ID_match_rate": None})
        df = pd.concat([df, extra], ignore_index=True)

    low_match = df["GO_ID_match_rate"].dropna()
    low_match = low_match[low_match < id_threshold]
    if len(low_match):
        print(f"ERROR: {len(low_match)} species have FASTA<->GO id match rate below {id_threshold:.0%} "
              f"(risk #1 of this project — gene models don't correspond). See {results_path}.", file=sys.stderr)
        df.to_csv(results_path, sep="\t", index=False)
        sys.exit(1)

    df.to_csv(results_path, sep="\t", index=False)
    n_pass = (df["BUSCO_status"] == "PASS").sum()
    n_flag = (df["BUSCO_status"] == "FLAG").sum()
    n_fail = (df["BUSCO_status"] == "FAIL").sum()
    n_not_resolved = (df["Status"] != "resolved").sum()
    _log(f"  wrote {results_path.name} — BUSCO PASS={n_pass} FLAG={n_flag} FAIL={n_fail}, "
         f"{n_not_resolved} species never resolved past M1")
    return df


def check_disk(path: Path, estimate_gb: float) -> None:
    usage = shutil.disk_usage(path)
    free_gb = usage.free / (1024 ** 3)
    if free_gb < estimate_gb:
        print(f"ERROR: only {free_gb:.0f} GB free under {path}, estimated need is {estimate_gb:.0f} GB "
              f"(OrthoFinder blast + gene-tree output). Free space or lower --disk-estimate-gb to override.",
              file=sys.stderr)
        sys.exit(1)
    _log(f"  disk check: {free_gb:.0f} GB free under {path} (estimate: {estimate_gb:.0f} GB) — OK")


# --------------------------------------------------------- M4: core selection
def select_core(qc_df: pd.DataFrame, lineage_df: pd.DataFrame, core_size: int) -> pd.DataFrame:
    """Stratified by order/family (APG IV via species_lineage.tsv), best-BUSCO
    within each stratum -- not top-N BUSCO overall, which would concentrate
    the core in a few well-studied families (Brassicaceae, Poaceae)."""
    merged = qc_df.merge(lineage_df, on="Species", how="left")
    merged = merged[merged["BUSCO_status"].isin(["PASS", "FLAG", "SKIPPED"])].copy()
    merged["stratum"] = merged["order"].fillna(merged["family"]).fillna("unknown")
    merged["_rank_score"] = merged["BUSCO_C"].fillna(-1)

    chosen = []
    for stratum, group in merged.groupby("stratum"):
        best = group.sort_values("_rank_score", ascending=False).iloc[0]
        chosen.append(best)
    core = pd.DataFrame(chosen).sort_values("_rank_score", ascending=False)

    if len(core) > core_size:
        core = core.head(core_size)
    elif len(core) < core_size:
        remaining = merged[~merged["Species"].isin(core["Species"])].sort_values("_rank_score", ascending=False)
        core = pd.concat([core, remaining.head(core_size - len(core))])

    core = core[["Species", "Code5", "order", "family", "BUSCO_C", "stratum"]].rename(
        columns={"stratum": "Selection_reason"})
    return core.reset_index(drop=True)


def run_module4(qc_df: pd.DataFrame, lineage_path: Path, results_path: Path, core_size: int,
                override_path, force: bool) -> pd.DataFrame:
    if _checkpoint(results_path, "core species", force):
        return pd.read_csv(results_path, sep="\t")

    if override_path is not None:
        core = pd.read_csv(override_path, sep="\t")
        _log(f"  --core-species-override set: using {len(core)} manually-specified species (bypasses stratified selection)")
        core.to_csv(results_path, sep="\t", index=False)
        return core

    if not lineage_path.exists():
        print(f"ERROR: --lineage-file not found: {lineage_path} (needed for APG IV stratified core selection, "
              f"or pass --core-species-override for a manual list)", file=sys.stderr)
        sys.exit(1)
    lineage_df = pd.read_csv(lineage_path, sep="\t")[["Species", "order", "family"]]
    core = select_core(qc_df, lineage_df, core_size)
    core.to_csv(results_path, sep="\t", index=False)
    _log(f"  wrote {results_path.name} ({len(core)} core species, {core['Selection_reason'].nunique()} strata)")
    return core


# --------------------------------------------------- M5/M6: OrthoFinder wrappers
def find_orthofinder_results(of_dir: Path):
    """Deepest Results_* dir under of_dir with a completed Orthogroups/ output.

    Every "-b" resume of an interrupted run makes OrthoFinder nest its next
    attempt one level deeper (Results_X/WorkingDirectory/OrthoFinder/Results_Y/...),
    so the real output is not necessarily a direct child of of_dir — search
    recursively and require Orthogroups/ to rule out empty/interrupted shells.
    """
    hits = [d for d in of_dir.rglob("Results_*") if (d / "Orthogroups").is_dir()]
    if not hits:
        return None
    return max(hits, key=lambda d: len(d.parts))


def _find_latest_orthofinder_attempt(of_dir: Path):
    """Deepest Results_* dir under of_dir, complete or not — used to locate a
    WorkingDirectory with reusable Blast*.txt.gz for a "-b" resume."""
    hits = list(of_dir.rglob("Results_*"))
    if not hits:
        return None
    return max(hits, key=lambda d: len(d.parts))


def _flatten_orthofinder_results(of_dir: Path, result: Path) -> Path:
    """Relocate a nested "-b"-resume Results_* dir to be a direct child of
    of_dir. OrthoFinder's --assign/--core expects a plain, non-nested
    Results dir; repeated interrupted resumes otherwise leave the real
    output buried several "WorkingDirectory/OrthoFinder/Results_*" levels
    deep, which OrthoFinder's own --core path resolution can't follow
    ("Couldn't find previous orthogroups")."""
    if result is None or result.parent == of_dir:
        return result
    flat = of_dir / result.name
    if flat.exists() and flat != result:
        shutil.rmtree(flat)
    _log(f"  flattening nested {result} -> {flat}")
    shutil.move(str(result), str(flat))
    return flat


def run_module5(core_codes: list, clean_dir: Path, of_dir: Path, threads: int, force: bool) -> Path:
    existing = find_orthofinder_results(of_dir)
    if existing and not force:
        existing = _flatten_orthofinder_results(of_dir, existing)
        _log(f"  [checkpoint] core inference — {existing} already complete, skipping")
        return existing

    orthofinder = _require_tool("orthofinder")
    working_dir = None
    latest_attempt = _find_latest_orthofinder_attempt(of_dir)
    if latest_attempt:
        wd = latest_attempt / "WorkingDirectory"
        blast_files = list(wd.glob("Blast*.txt.gz")) if wd.is_dir() else []
        if blast_files:
            working_dir = wd

    # of_dir itself must NOT be pre-created: OrthoFinder's "-f" (fresh run)
    # mode requires -o to not exist yet, it creates it. On a "-b" resume,
    # of_dir already exists (that's how `existing` was found above).
    core_proteomes = of_dir.parent / "core_proteomes"
    core_proteomes.mkdir(parents=True, exist_ok=True)
    for code in core_codes:
        dest = core_proteomes / f"{code}.fa"
        if not dest.exists():
            shutil.copy2(clean_dir / f"{code}.fa", dest)

    cmd = ["orthofinder", "-t", str(threads), "-a", str(max(1, threads // 4)),
           "-S", "diamond_ultra_sens", "-M", "msa", "-A", "famsa", "-T", "fasttree",
           "-o", str(of_dir)]
    if working_dir:
        _log(f"  resuming from existing all-vs-all in {working_dir}")
        cmd = ["orthofinder", "-b", str(working_dir)] + cmd[1:-2]
    else:
        if of_dir.exists():
            _log(f"  removing incomplete {of_dir} from a previous interrupted run")
            shutil.rmtree(of_dir)
        cmd = ["orthofinder", "-f", str(core_proteomes)] + cmd[1:]
    _run(cmd)
    return _flatten_orthofinder_results(of_dir, find_orthofinder_results(of_dir))


def run_module6(rest_codes: list, clean_dir: Path, core_results: Path, of_core_dir: Path,
                 assign_dir: Path, threads: int, force: bool) -> Path:
    existing = find_orthofinder_results(assign_dir)
    if existing and not force:
        _log(f"  [checkpoint] assignment — {existing} already complete, skipping")
        return existing

    _require_tool("orthofinder")
    rest_proteomes = assign_dir.parent / "rest_proteomes"
    rest_proteomes.mkdir(parents=True, exist_ok=True)
    for code in rest_codes:
        dest = rest_proteomes / f"{code}.fa"
        if not dest.exists():
            shutil.copy2(clean_dir / f"{code}.fa", dest)

    if assign_dir.exists():
        _log(f"  removing incomplete {assign_dir} from a previous interrupted run")
        shutil.rmtree(assign_dir)
    assign_dir.mkdir(parents=True, exist_ok=True)

    # --assign rejects -o ("only with -f"); it writes its Results_* dir
    # somewhere under --core's tree instead (nesting depth depends on how
    # deep core_results itself is). Snapshot before/after under the stable
    # of_core_dir root to find it regardless of depth, then move it under
    # assign_dir so find_orthofinder_results(assign_dir) still works.
    before = set(of_core_dir.rglob("Results_*"))
    _run(["orthofinder", "--assign", str(rest_proteomes), "--core", str(core_results),
          "-t", str(threads)])
    new_dirs = set(of_core_dir.rglob("Results_*")) - before
    if not new_dirs:
        print(f"ERROR: OrthoFinder --assign produced no new Results_* directory under {of_core_dir}",
              file=sys.stderr)
        sys.exit(1)
    new_dir = max(new_dirs, key=lambda d: len(d.parts))
    shutil.move(str(new_dir), str(assign_dir / new_dir.name))

    result = find_orthofinder_results(assign_dir)
    if result is None:
        # OrthoFinder can exit 0 even after an internal "ERROR:" (e.g.
        # "Couldn't find previous orthogroups"), so _run()'s returncode
        # check alone won't catch this — the new dir exists but has no
        # Orthogroups/, meaning the assignment silently failed.
        print(f"ERROR: orthofinder --assign finished without producing Orthogroups/ in "
              f"{assign_dir / new_dir.name} — check the OrthoFinder output above for an "
              f"'ERROR:' line; the run failed internally despite exiting 0.", file=sys.stderr)
        sys.exit(1)
    return result


# ------------------------------------------------------- M7: species x OG matrix
def build_og_matrix(genecount_tsv: Path, min_species: int) -> pd.DataFrame:
    df = pd.read_csv(genecount_tsv, sep="\t", index_col=0)
    if "Total" in df.columns:
        df = df.drop(columns=["Total"])
    matrix = df.T  # species x OG
    n_species_per_og = (matrix > 0).sum(axis=0)
    keep = n_species_per_og >= min_species
    n_dropped = (~keep).sum()
    matrix = matrix.loc[:, keep]
    _log(f"  --min-species={min_species}: dropped {n_dropped}/{len(keep)} species-specific-ish OGs, "
         f"{matrix.shape[1]} remain")
    matrix.index.name = "Species"
    return matrix


def run_module7(core_results: Path, assign_results: Path, results_matrix: Path, results_stats: Path,
                 min_species: int, force: bool) -> pd.DataFrame:
    if _checkpoint(results_matrix, "OG matrix", force):
        return pd.read_csv(results_matrix, sep="\t", index_col=0)

    genecount = assign_results / "Orthogroups" / "Orthogroups.GeneCount.tsv" if assign_results else None
    if not genecount or not genecount.exists():
        genecount = core_results / "Orthogroups" / "Orthogroups.GeneCount.tsv"
    matrix = build_og_matrix(genecount, min_species)
    matrix.to_csv(results_matrix, sep="\t")

    species_stats = pd.DataFrame({
        "Species": matrix.index,
        "Total_prots": matrix.sum(axis=1).values,
        "N_orthogroups": (matrix > 0).sum(axis=1).values,
    })
    species_stats.to_csv(results_stats, sep="\t", index=False)
    _log(f"  wrote {results_matrix.name} ({matrix.shape[0]} species x {matrix.shape[1]} OGs) and {results_stats.name}")
    return matrix


# ---------------------------------------------------------- M8: OG annotation
def parse_orthogroup_members(orthogroups_tsv: Path) -> dict:
    """OG id -> list of protein ids, from Orthogroups.tsv (OG, species1, species2, ...; cells = comma-space-separated ids)."""
    df = pd.read_csv(orthogroups_tsv, sep="\t", index_col=0)
    members = {}
    for og, row in df.iterrows():
        ids = []
        for cell in row.dropna():
            ids += [p.strip() for p in str(cell).split(",") if p.strip()]
        members[og] = ids
    return members


def parse_go_sets(go_file: Path) -> dict:
    """protein_id -> set(GO ids), from a `<prefix>_GOs_merged.tsv` topgo file."""
    out = {}
    with open(go_file) as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 2 or not parts[1].strip():
                continue
            out[parts[0].strip()] = {g.strip() for g in parts[1].split(",") if g.strip()}
    return out


def consensus_annotation(og_members: dict, protein_go: dict, threshold: float) -> pd.DataFrame:
    """A GO is assigned to an OG if >= threshold of its members (that have any
    GO annotation at all) carry it -- consensus over the orthogroup dampens
    the protein-to-protein spurious-transfer bias documented in go_phylostrata.py."""
    rows = []
    for og, members in og_members.items():
        bare_ids = [m.split("|", 1)[-1] for m in members]
        annotated = [protein_go[i] for i in bare_ids if i in protein_go]
        n_genes = len(members)
        n_species = len({m.split("|", 1)[0] for m in members})
        if not annotated:
            rows.append({"Orthogroup": og, "N_species": n_species, "N_genes": n_genes, "GO_consensus": ""})
            continue
        go_counts = Counter(g for gos in annotated for g in gos)
        n_annot = len(annotated)
        consensus = sorted(
            ((g, c / n_annot) for g, c in go_counts.items() if c / n_annot >= threshold),
            key=lambda x: -x[1],
        )
        # "|" (not ":") separates id from fraction -- GO ids themselves contain a colon (GO:0016787).
        rows.append({"Orthogroup": og, "N_species": n_species, "N_genes": n_genes,
                     "GO_consensus": ";".join(f"{g}|{frac:.2f}" for g, frac in consensus)})
    return pd.DataFrame(rows)


def write_ic_bridge_file(annotation_df: pd.DataFrame, go_desc: dict, out_path: Path) -> None:
    """Headerless TSV shaped exactly like data/All_GOs_ic.tsv so general_pca_common.load_go_ic_and_descriptions
    (parts[0]=id, parts[4]=IC-or-'-', parts[5]=description) reads it unchanged -- see design doc M8."""
    with open(out_path, "w") as fh:
        for _, row in annotation_df.iterrows():
            terms = [t for t in row["GO_consensus"].split(";") if t]
            descs = []
            for t in terms[:5]:
                go_id, frac = t.split("|")
                pct = round(float(frac) * 100)
                descs.append(f"{go_desc.get(go_id, go_id)} [{go_id}, {pct}% sp]")
            desc = "; ".join(descs) + f" ({row['N_species']} sp / {row['N_genes']} genes)" if descs else \
                f"unannotated orthogroup ({row['N_species']} sp / {row['N_genes']} genes)"
            fh.write(f"{row['Orthogroup']}\torthogroup\t{row['N_species']}\t{row['N_genes']}\t-\t{desc}\t\n")


def run_module8(matrix_codes: list, code5_to_go: dict, orthogroups_tsv: Path,
                 results_annotation: Path, results_bridge: Path, threshold: float, force: bool) -> pd.DataFrame:
    """matrix_codes: the Code5 values that are matrix.index (OrthoFinder uses
    Code5 as column/species names, not the full Species name) -- code5_to_go
    maps each one to its GOFile, since GO tables are keyed by original
    protein id, not the <Code5>| prefixed one."""
    if _checkpoint(results_annotation, "OG annotation", force):
        return pd.read_csv(results_annotation, sep="\t")

    protein_go = {}
    n_with_go = 0
    for code in matrix_codes:
        go_file = code5_to_go.get(code)
        if go_file and pd.notna(go_file) and Path(go_file).exists():
            protein_go.update(parse_go_sets(Path(go_file)))
            n_with_go += 1
    _log(f"  loaded GO annotations for {n_with_go}/{len(matrix_codes)} species ({len(protein_go)} annotated proteins)")

    og_members = parse_orthogroup_members(orthogroups_tsv)
    annotation_df = consensus_annotation(og_members, protein_go, threshold)
    annotation_df.to_csv(results_annotation, sep="\t", index=False)

    go_desc = {}
    for gos in protein_go.values():
        for g in gos:
            go_desc.setdefault(g, g)
    write_ic_bridge_file(annotation_df, go_desc, results_bridge)
    _log(f"  wrote {results_annotation.name} ({len(annotation_df)} orthogroups) and {results_bridge.name}")
    return annotation_df
