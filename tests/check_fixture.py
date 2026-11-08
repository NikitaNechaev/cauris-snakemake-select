"""Check scientific output contracts on the synthetic fixture."""
import argparse
import csv
import hashlib
import json
import zipfile
from pathlib import Path


def rows(path, key):
    with path.open(encoding="utf-8", newline="") as handle:
        data = list(csv.DictReader(handle, delimiter="\t"))
    index = {row[key]: row for row in data}
    assert len(index) == len(data), path
    return index


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("results", type=Path)
    parser.add_argument("--archived", action="store_true")
    parser.add_argument("--legacy", action="store_true")
    args = parser.parse_args()
    root = args.results
    final = root / "12_final_report"
    summary = json.loads((final / "final_summary.json").read_text(encoding="utf-8"))
    assert summary["n_candidates"] == 3
    assert summary["n_genomes"] == 2
    assert summary["clade_counts"]["I"] == summary["clade_counts"]["II"] == 1
    assert summary["n_human_reported_hits"] == 1
    candidates = rows(final / "integrated_candidates.tsv", "accession")
    assert set(candidates) == {"TEST001", "TEST002", "TEST003"}
    assert [candidates[c]["n_QC_extracted"] for c in sorted(candidates)] == ["2", "2", "0"]
    assert candidates["TEST001"]["human_has_reported_hit"] == "True"
    aligned = rows(root / "10_mafft/alignment_summary.tsv", "candidate")
    assert [aligned[c]["status"] for c in sorted(aligned)] == [
        "aligned", "single_unique_allele", "no_extracted_sequences",
    ]
    assert aligned["TEST001"]["n_unique_alleles"] == "2"

    population = root / "08_population/results"
    comparisons = rows(population / "ani_population/all_estimated_pairs.tsv", "assembly_1")
    assert len(comparisons) == 1
    pair = next(iter(comparisons.values()))
    assert float(pair["ANI_pct"]) > 99
    if args.archived:
        assert float(pair["AF_min_pct"]) == 97.5
        assert {float(pair["AF_1_to_2_pct"]), float(pair["AF_2_to_1_pct"])} == {97.5, 98.5}
    clades = rows(population / "clade_typing/clade_diagnostics.tsv", "assembly")
    assert all(row["nearest_matches_control"] == "YES" for row in clades.values())
    isolates = rows(root / "07_annotated_isolates/summary/candidate_presence_summary.tsv", "candidate")
    assert [isolates[c]["n_assemblies_strict_hit"] for c in sorted(isolates)] == ["2", "2", "0"]
    if args.legacy:
        legacy = rows(root / "13_legacy_population/summary/candidate_presence_summary.tsv", "candidate")
        assert [legacy[c]["n_pass"] for c in sorted(legacy)] == ["2", "2", "0"]

    provenance = rows(final / "raw_inputs_sha256.tsv", "path")
    for name, record in provenance.items():
        source = Path(name)
        assert source.stat().st_size == int(record["bytes"])
        assert hashlib.sha256(source.read_bytes()).hexdigest() == record["sha256"]
    with zipfile.ZipFile(final / "final_results_bundle.zip") as archive:
        assert archive.testzip() is None
        names = set(archive.namelist())
        assert {"integrated_candidates.tsv", "candidate_by_clade.tsv", "raw_inputs_sha256.tsv", "RESULTS_RU.md"} <= names
        assert all(not n.endswith((".faa", ".fna", ".gz")) for n in names)
    print("Synthetic results and provenance verified.")


if __name__ == "__main__":
    main()
