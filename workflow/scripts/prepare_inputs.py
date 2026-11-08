"""Validate and stage archived inputs, candidate metadata and provenance."""
from __future__ import annotations

import argparse
import csv
import json
import math
import shutil
from collections import Counter, defaultdict
from pathlib import Path

from cauris_population_stage import digest_file, read_manifest

CLADES = {"I", "II", "III", "IV", "V", "VI", "Unassigned"}


def read_tsv(path, required):
    with path.open(encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not set(required).issubset(reader.fieldnames or ()):
            raise ValueError(f"{path}: expected columns {', '.join(required)}")
        return reader.fieldnames, list(reader)


def write_tsv(path, fields, rows):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def snapshot(args):
    args.out.mkdir(parents=True, exist_ok=True)
    for source in args.files:
        shutil.copyfile(source, args.out / source.name)


def candidates(args):
    fields, rows = read_tsv(
        args.qc / "surface_priority.tsv",
        ("accession", "organism_id", "protein_name", "reviewed", "evidence_group"),
    )
    by_accession = {row["accession"]: row for row in rows}
    if len(by_accession) != len(rows):
        raise ValueError("Duplicate accessions in the priority catalog")
    primary = {key for key, row in by_accession.items() if row["organism_id"].strip() == "498019"}
    if args.allowlist:
        allowed = [
            line.strip() for line in args.allowlist.read_text(encoding="utf-8-sig").splitlines()
            if line.strip() and not line.lstrip().startswith("#")
        ]
        if len(allowed) != len(set(allowed)):
            raise ValueError(f"Duplicate accession in {args.allowlist}")
        unknown = set(allowed) - primary
        if unknown:
            raise ValueError(f"Allowlist accessions outside the C. auris priority catalog: {sorted(unknown)}")
        primary = set(allowed)
    if not primary:
        raise ValueError("The C. auris priority candidate panel is empty")

    by_member = {}
    cluster_members = defaultdict(set)
    with args.clusters.open(encoding="utf-8", newline="") as handle:
        for number, pair in enumerate(csv.reader(handle, delimiter="\t"), 1):
            if len(pair) != 2 or not set(pair).issubset(by_accession):
                raise ValueError(f"Invalid MMseqs2 cluster row {number}: {pair}")
            representative, member = pair
            if member in by_member:
                raise ValueError(f"Repeated MMseqs2 cluster member: {member}")
            by_member[member] = representative
            cluster_members[representative].add(member)
    if set(by_member) != set(by_accession):
        raise ValueError(f"Cluster/catalog membership mismatch: {sorted(set(by_accession) - set(by_member))[:10]}")

    selected_rows = []
    for accession in sorted(primary):
        representative = by_member[accession]
        members = cluster_members[representative]
        taxids = {by_accession[m]["organism_id"] for m in members} - {""}
        selected_rows.append({
            **by_accession[accession],
            "cluster_id": representative,
            "cluster_members": len(members),
            "cluster_organism_taxids": ";".join(sorted(taxids)),
        })
    args.fasta.parent.mkdir(parents=True, exist_ok=True)
    written = set()
    with (args.qc / "surface_priority.faa").open(encoding="utf-8") as source, args.fasta.open("w", encoding="utf-8") as target:
        keep = False
        for line in source:
            if line.startswith(">"):
                accession = line[1:].split(maxsplit=1)[0].strip()
                keep = accession in primary
                if keep:
                    if accession in written:
                        raise ValueError(f"Repeated FASTA accession {accession}")
                    written.add(accession)
            if keep:
                target.write(line)
    if written != primary:
        raise ValueError(f"Missing candidate FASTA records: {sorted(primary - written)}")
    write_tsv(args.metadata, fields + ["cluster_id", "cluster_members", "cluster_organism_taxids"], selected_rows)
    print(f"Selected {len(primary)} catalog proteins across {len(cluster_members)} similarity clusters")


def clades(args):
    expected = {row["assembly"] for row in read_manifest(args.manifest)}
    fields, assignments = read_tsv(args.assignments, ("assembly", "assigned_clade", "assignment_status"))
    observed = {row["assembly"] for row in assignments}
    if len(assignments) != len(observed) or observed != expected:
        raise ValueError(f"Clade panel differs from genome manifest. Missing: {sorted(expected - observed)}; extra: {sorted(observed - expected)}")
    for row in assignments:
        row["assigned_clade"] = row["assigned_clade"].strip() or "Unassigned"
        if row["assigned_clade"] not in CLADES:
            raise ValueError(f"Unknown clade for {row['assembly']}: {row['assigned_clade']}")
    fields_cmp, comparisons = read_tsv(args.comparisons, (
        "assembly", "ref_clade", "reference_accession", "ANI_pct", "AF_ref_pct", "AF_query_pct",
    ))
    if not comparisons:
        raise ValueError("Reference comparison table is empty")
    seen = set()
    for row in comparisons:
        key = row["assembly"], row["reference_accession"]
        if key in seen or key[0] not in expected or not key[1] or row["ref_clade"] not in CLADES - {"Unassigned"}:
            raise ValueError(f"Invalid or duplicate reference comparison: {key}")
        seen.add(key)
        for field in ("ANI_pct", "AF_ref_pct", "AF_query_pct"):
            value = row[field].strip()
            if value in ("", "NA", "NaN", "nan"):
                row[field] = ""
                continue
            parsed = float(value)
            if not math.isfinite(parsed) or not 0 <= parsed <= 100:
                raise ValueError(f"{key}: {field} is outside 0–100")
    write_tsv(args.summary.parent / "clade_assignments.tsv", fields, assignments)
    write_tsv(args.summary.parent / "all_reference_comparisons.tsv", fields_cmp, comparisons)
    args.summary.write_text(json.dumps({
        "assignment_source": "curated input table",
        "clade_counts": dict(Counter(row["assigned_clade"] for row in assignments)),
        "reference_comparison_rows": len(comparisons),
        "assemblies_with_reference_comparisons": len({assembly for assembly, _ in seen}),
    }, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def provenance(args):
    root = Path.cwd()
    rows = []
    for path in sorted({Path(f).absolute() for f in args.files}):
        relative = path.relative_to(root) if path.is_relative_to(root) else path
        rows.append({"path": str(relative), "bytes": path.stat().st_size, "sha256": digest_file(path)})
    write_tsv(args.manifest, ["path", "bytes", "sha256"], rows)
    config = json.loads(args.config_json)
    args.config_out.write_text(json.dumps(config, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(required=True)
    sub = commands.add_parser("snapshot")
    sub.add_argument("--out", type=Path, required=True)
    sub.add_argument("files", nargs="+", type=Path)
    sub.set_defaults(run=snapshot)

    sub = commands.add_parser("candidates")
    for name in ("qc", "clusters", "fasta", "metadata"):
        sub.add_argument("--" + name, type=Path, required=True)
    sub.add_argument("--allowlist", type=Path)
    sub.set_defaults(run=candidates)

    sub = commands.add_parser("clades")
    for name in ("manifest", "assignments", "comparisons", "summary"):
        sub.add_argument("--" + name, type=Path, required=True)
    sub.set_defaults(run=clades)

    sub = commands.add_parser("provenance")
    for name in ("manifest", "config-out"):
        sub.add_argument("--" + name, type=Path, required=True)
    sub.add_argument("--config-json", required=True)
    sub.add_argument("files", nargs="+", type=Path)
    sub.set_defaults(run=provenance)
    args = parser.parse_args()
    args.run(args)


if __name__ == "__main__":
    main()
