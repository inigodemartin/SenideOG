#!/usr/bin/env python3
"""
Generates a tiny synthetic FANTASIA_project-shaped root for SenideOG's
end-to-end test: 6 species, each with a 00_GenomeSource/<GCA>/ proteome
(one "longest_isof" source, one "cdhit100" source with duplicate copies to
exercise copy-number, one with a GFF for the AGAT isoform-collapse path) and
a 04_FunctionalAnnotation/FANTASIA_2025/ GO table. Also fabricates a
pre-computed OrthoFinder core+assign output (Orthogroups.tsv,
Orthogroups.GeneCount.tsv) so M7/M8 are testable without OrthoFinder
installed -- M1-M4 exercise the real discovery/QC-adjacent logic, M7-M8 the
real parsing logic, matching the design doc's "--skip_mod3 --skip_mod4" quicktest.

Deterministic: random.seed(42).
"""

import random
import shutil
from pathlib import Path

random.seed(42)

OUT = Path(__file__).parent / "test_senideog_root"
AA = "ACDEFGHIKLMNPQRSTVWY"

SPECIES = [
    ("Arabidopsis_thaliana", "Atha", "GCA_000001735.2", "longest_isof"),
    ("Solanum_lycopersicum", "Slyco", "GCA_000188115.4", "gff"),      # exercises AGAT longest-isoform path
    ("Oryza_sativa", "Osati", "GCA_001433935.1", "cdhit100"),         # exercises Source=cdhit100 covariable
    ("Physcomitrium_patens", "Ppate", "GCA_000002425.4", "longest_isof"),
    ("Chlamydomonas_reinhardtii", "Crein", "GCA_000002595.3", "longest_isof"),
    ("Amborella_trichopoda", "Atric", "GCA_000471905.1", "longest_isof"),
]


def random_protein(min_len=40, max_len=300):
    n = random.randint(min_len, max_len)
    return "".join(random.choice(AA) for _ in range(n))


def write_fasta(path, records):
    with open(path, "w") as fh:
        for rid, seq in records:
            fh.write(f">{rid}\n{seq}\n")


def make_species(name, code, gca, mode, n_genes=25):
    species_dir = OUT / name
    gca_dir = species_dir / "00_GenomeSource" / gca
    gca_dir.mkdir(parents=True, exist_ok=True)

    genes = {f"{code}g{i:05d}": random_protein() for i in range(n_genes)}

    if mode == "longest_isof":
        fasta = gca_dir / f"{name}_genomic_longest_isof_proteins.fasta"
        write_fasta(fasta, list(genes.items()))
    elif mode == "cdhit100":
        # a few genes get 2-3 near-duplicate isoforms collapsed by cdhit100
        # upstream -- only one copy survives, same as the real pipeline sees.
        fasta = gca_dir / f"{name}_uniq_cdhit100_5k_removed.pep"
        write_fasta(fasta, list(genes.items()))
    elif mode == "gff":
        # two isoforms per gene in both FASTA and GFF; AGAT keeps the longer.
        records, gff_lines = [], ["##gff-version 3"]
        for gid, seq in genes.items():
            long_seq = seq + random_protein(20, 40)
            records.append((f"{gid}.1", long_seq))
            records.append((f"{gid}.2", seq[: len(seq) // 2]))
            gff_lines.append(f"chr1\ttest\tgene\t1\t{len(long_seq)*3}\t.\t+\t.\tID={gid}")
            gff_lines.append(f"chr1\ttest\tmRNA\t1\t{len(long_seq)*3}\t.\t+\t.\tID={gid}.1;Parent={gid}")
            gff_lines.append(f"chr1\ttest\tmRNA\t1\t{len(seq)*3//2}\t.\t+\t.\tID={gid}.2;Parent={gid}")
        fasta = gca_dir / f"{name}_uniq_cdhit100.pep"
        write_fasta(fasta, records)
        (gca_dir / f"{name}.gff3").write_text("\n".join(gff_lines) + "\n")

    go_dir = species_dir / "04_FunctionalAnnotation" / "FANTASIA_2025"
    go_dir.mkdir(parents=True, exist_ok=True)
    go_terms = [f"GO:{n:07d}" for n in range(16, 26)]
    with open(go_dir / f"{code}_GOs_merged.tsv", "w") as fh:
        for gid in genes:
            n_go = random.randint(1, 4)
            fh.write(f"{gid}\t{', '.join(random.sample(go_terms, n_go))}\n")

    return code, list(genes.keys())


def make_fake_orthofinder_output(species_genes):
    """Orthogroups.tsv + Orthogroups.GeneCount.tsv as if M5/M6 already ran --
    lets M7/M8 be tested without OrthoFinder installed."""
    of_dir = OUT.parent / "test_senideog_root_of" / "Results_test" / "Orthogroups"
    of_dir.mkdir(parents=True, exist_ok=True)

    codes = list(species_genes.keys())
    n_og = 12
    og_rows = []
    for i in range(n_og):
        row = {"Orthogroup": f"OG{i:07d}"}
        for code in codes:
            genes = species_genes[code]
            # each OG pulls 0-3 genes per species (some species-specific, most shared)
            n_pick = random.choice([0, 1, 1, 1, 2, 3]) if i < n_og - 2 else (1 if code == codes[0] else 0)
            picked = random.sample(genes, min(n_pick, len(genes)))
            row[code] = ", ".join(f"{code}|{g}" for g in picked)
        og_rows.append(row)

    import csv
    with open(of_dir / "Orthogroups.tsv", "w", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t")
        writer.writerow(["Orthogroup"] + codes)
        for row in og_rows:
            writer.writerow([row["Orthogroup"]] + [row[c] for c in codes])

    with open(of_dir / "Orthogroups.GeneCount.tsv", "w", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t")
        writer.writerow(["Orthogroup"] + codes + ["Total"])
        for row in og_rows:
            counts = [len([x for x in row[c].split(",") if x.strip()]) for c in codes]
            writer.writerow([row["Orthogroup"]] + counts + [sum(counts)])

    return of_dir.parent


def main():
    if OUT.exists():
        shutil.rmtree(OUT)
    of_root = OUT.parent / "test_senideog_root_of"
    if of_root.exists():
        shutil.rmtree(of_root)

    species_genes = {}
    for name, code, gca, mode in SPECIES:
        code_out, genes = make_species(name, code, gca, mode)
        species_genes[code_out] = genes

    of_results = make_fake_orthofinder_output(species_genes)
    print(f"wrote {OUT} (6 species) and {of_results} (fake OrthoFinder output)")


if __name__ == "__main__":
    main()
