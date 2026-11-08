"""Generate deterministic, entirely synthetic inputs for the workflow smoke test."""
from __future__ import annotations

import csv
import gzip
import hashlib
import json
import random
import textwrap
from pathlib import Path

AA = "ACDEFGHIKLMNPQRSTVWY"
CODON = dict(zip(AA, "GCT TGT GAT GAA TTT GGT CAT ATT AAA TTA ATG AAT CCT CAA CGT CTG ACT GTT TGG TAT".split()))
ASSEMBLIES = ("GCA_999999991.1", "GCA_999999992.1")
ROOT = Path("tests/_work")
DATA = ROOT / "inputs"
CATALOG = DATA / "uniprot"
GENOMES = DATA / "ncbi" / "data"


def tsv(path, headers, rows):
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(headers)
        writer.writerows(rows)


def fasta(records):
    return "".join(f">{name}\n" + "\n".join(textwrap.wrap(seq, 80)) + "\n" for name, seq in records)


def main():
    rng = random.Random(20261108)
    CATALOG.mkdir(parents=True, exist_ok=True)
    p1 = "M" + "".join(rng.choices(AA, k=149))
    p1_variant = p1[:70] + ("A" if p1[70] != "A" else "V") + p1[71:]
    p2 = "M" + "".join(rng.choices(AA, k=149))
    p3 = "M" + "".join(rng.choices(AA, k=149))
    proteins = {"TEST001": p1, "TEST002": p2, "TEST003": p3}
    (CATALOG / "surface_candidates.faa").write_text(fasta(proteins.items()), encoding="utf-8")
    tsv(CATALOG / "surface_candidates.tsv",
        ["accession", "protein_name", "organism_name", "organism_id", "reviewed", "evidence_group",
         "length", "sequence_length", "sequence_sha256", "matched_categories"],
        [[key, f"Synthetic protein {key}", "Candida auris", "498019", "unreviewed",
          "cell_wall_annotated", len(seq), len(seq), hashlib.sha256(seq.encode()).hexdigest(), "synthetic"]
         for key, seq in proteins.items()])
    (CATALOG / "run_metadata.json").write_text('{"synthetic":true,"categories_failed":0}\n', encoding="utf-8")
    (DATA / "human.faa").write_text(fasta([("SYNTHETIC_REFERENCE", p1)]), encoding="utf-8")
    (DATA / "selected.txt").write_text("\n".join(ASSEMBLIES) + "\n", encoding="utf-8")
    tsv(DATA / "clades.tsv", ["assembly", "assigned_clade", "assignment_status", "known_control_clade"],
        [(ASSEMBLIES[0], "I", "synthetic", "I"), (ASSEMBLIES[1], "II", "synthetic", "II")])
    tsv(DATA / "comparisons.tsv",
        ["assembly", "ref_clade", "reference_accession", "ANI_pct", "AF_ref_pct", "AF_query_pct"],
        [(a, cl, f"SYNTHETIC_{cl}", 100 if idx == ref else 99, 95, 94)
         for idx, a in enumerate(ASSEMBLIES) for ref, cl in enumerate(("I", "II"))])

    background = "".join(rng.choices("ACGT", k=260_000))
    for index, (assembly, protein) in enumerate(zip(ASSEMBLIES, (p1, p1_variant))):
        folder = GENOMES / assembly
        folder.mkdir(parents=True, exist_ok=True)
        genome = list(background)
        for pos in rng.sample(range(20_000, len(genome)), k=index * 200):
            genome[pos] = rng.choice("ACGT")
        for start, seq in ((10_000, protein), (15_000, p2)):
            coding = "".join(CODON[aa] for aa in seq) + "TAA"
            genome[start:start + len(coding)] = coding
        content = fasta([(f"contig_{index + 1}", "".join(genome))])
        if index:
            (folder / f"{assembly}_genomic.fna.gz").write_bytes(gzip.compress(content.encode(), mtime=0))
        else:
            (folder / f"{assembly}_genomic.fna").write_text(content, encoding="utf-8")
        (folder / "protein.faa").write_text(fasta([("mock_protein_1", protein), ("mock_protein_2", p2)]), encoding="utf-8")

    for name, pairs in (("ani.txt", (99.7, 99.7)), ("af.txt", (97.5, 98.5))):
        (DATA / name).write_text(f"2\n{ASSEMBLIES[0]}\t100\t{pairs[0]}\n{ASSEMBLIES[1]}\t{pairs[1]}\t100\n", encoding="utf-8")

    config = {
        "results": str(ROOT / "results"),
        "catalog": {"directory": str(CATALOG)},
        "inputs": {
            "selected": str(DATA / "selected.txt"),
            "genomes": str(GENOMES),
            "human_proteome": str(DATA / "human.faa"),
            "clade_assignments": str(DATA / "clades.tsv"),
            "clade_comparisons": str(DATA / "comparisons.tsv"),
            "isolates": str(GENOMES),
        },
        "threads": {"search": 2, "mapping": 1, "mafft": 2},
        "analysis": {"min_clade_n": 2},
    }
    (ROOT / "config.json").write_text(json.dumps(config, indent=2) + "\n", encoding="utf-8")
    config["results"] = str(ROOT / "results-archived")
    config["inputs"].update(ani_matrix=str(DATA / "ani.txt"), af_matrix=str(DATA / "af.txt"))
    (ROOT / "config-archived.json").write_text(json.dumps(config, indent=2) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
