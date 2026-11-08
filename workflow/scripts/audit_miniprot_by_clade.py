#!/usr/bin/env python3
"""TEN / Candida auris: step 17, QC audit of miniprot protein-to-genome alignments.

Inputs (relative to the TEN project root):
  08_population/results/candidate_genome_details.tsv
  08_population/results/clade_typing/clade_assignments.tsv
  08_population/qc/genome_manifest.tsv
  06_blastp_on_human/auris_cluster_audit.tsv (optional names/evidence)
  08_population/qc/genome_qc.tsv (optional genome-level assembly metrics)

Outputs: 08_population/results/mapping_audit/
No external Python libraries. This is a QC audit, not a test of orthology,
protein expression, cell-surface exposure, immunogenicity, or vaccine protection.

Examples:
  python scripts/audit_miniprot_by_clade.py
  python scripts/audit_miniprot_by_clade.py --edge-check --edge-window 1500
"""
from __future__ import annotations

import argparse
import csv
import gzip
import json
import math
import re
import sys
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path

CLADES = ("I", "II", "III", "IV", "V", "VI")
VALID_STATUSES = {"pass", "partial_or_divergent", "no_alignment"}
Z = 1.959963984540054


def read_tsv(path: Path, required: set[str]) -> list[dict[str, str]]:
    if not path.is_file():
        raise FileNotFoundError(f"Input file not found: {path}")
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        missing = required - set(reader.fieldnames or [])
        if missing:
            raise ValueError(f"{path}: missing columns {', '.join(sorted(missing))}")
        rows = list(reader)
    if not rows:
        raise ValueError(f"Empty table: {path}")
    return rows


def maybe_tsv(path: Path, required: set[str]) -> list[dict[str, str]]:
    return read_tsv(path, required) if path.is_file() else []


def write_tsv(path: Path, fields: list[str], rows: list[dict]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def number(value, default=None):
    try:
        if value is None or str(value).strip() == "":
            return default
        return float(value)
    except (ValueError, TypeError):
        return default


def integer(value, default=None):
    v = number(value)
    if v is None or not math.isfinite(v) or int(v) != v:
        return default
    return int(v)


def wilson(n: int, total: int) -> tuple[str, str]:
    if total == 0:
        return "NA", "NA"
    p = n / total
    z2 = Z ** 2
    mid = (p + z2 / (2 * total)) / (1 + z2 / total)
    half = (Z * math.sqrt((p * (1 - p) + z2 / (4 * total)) / total)) / (1 + z2 / total)
    return f"{100 * max(0, mid - half):.2f}", f"{100 * min(1, mid + half):.2f}"


def fasta_lengths(path: Path) -> dict[str, int]:
    """Optional contig edge check; sequence is streamed, not kept in memory."""
    result: dict[str, int] = {}
    curr = None
    length = 0
    opener = gzip.open if path.suffix.lower() == ".gz" else open
    with opener(path, "rt", encoding="utf-8") as handle:
        for raw in handle:
            if raw.startswith(">"):
                if curr is not None:
                    if curr in result:
                        raise ValueError(f"Duplicate contig {curr} in {path}")
                    result[curr] = length
                curr = raw[1:].split()[0]
                length = 0
            elif raw.strip():
                length += len(raw.strip())
        if curr is not None:
            if curr in result:
                raise ValueError(f"Duplicate contig {curr} in {path}")
            result[curr] = length
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("--root", type=Path, default=None, help="TEN project root, auto-detected as parent of scripts/")
    parser.add_argument("--details", type=Path, help="miniprot candidate_genome_details.tsv")
    parser.add_argument("--clades", type=Path, help="clade_assignments.tsv, after step 16 reclassification")
    parser.add_argument("--manifest", type=Path, help="genome_manifest.tsv")
    parser.add_argument("--metadata", type=Path, help="optional auris_cluster_audit.tsv")
    parser.add_argument("--out", type=Path, help="output directory")
    parser.add_argument("--edge-check", action="store_true", help="stream FASTA to check genomic loci close to contig ends")
    parser.add_argument("--edge-window", type=int, default=1500, help="bp from contig end to flag as uncertain")
    parser.add_argument("--min-clade-n", type=int, default=5, help="below this n, note sparse clade sampling")
    parser.add_argument("--near-identity", type=float, default=5, help="pp above mapping identity threshold considered borderline")
    parser.add_argument("--near-coverage", type=float, default=5, help="pp above mapping coverage threshold considered borderline")
    args = parser.parse_args()
    root = (args.root or Path(__file__).resolve().parent.parent).expanduser().resolve()
    pop = root / "08_population"
    details_path = args.details or pop / "results/candidate_genome_details.tsv"
    clades_path = args.clades or pop / "results/clade_typing/clade_assignments.tsv"
    manifest_path = args.manifest or pop / "qc/genome_manifest.tsv"
    metadata_path = args.metadata or root / "06_blastp_on_human/auris_cluster_audit.tsv"
    out = args.out or pop / "results/mapping_audit"
    mapping_summary_path = details_path.parent / "summary.json"
    mapping_thresholds = {"identity_pct": 40.0, "candidate_coverage_pct": 70.0}
    if mapping_summary_path.is_file():
        mapping_meta = json.loads(mapping_summary_path.read_text(encoding="utf-8"))
        mapping_thresholds.update(mapping_meta.get("thresholds", {}))
    min_identity = float(mapping_thresholds["identity_pct"])
    min_coverage = float(mapping_thresholds["candidate_coverage_pct"])
    details = read_tsv(details_path, {"assembly", "candidate", "status", "n_alignments", "n_strict_alignments",
                                     "best_identity_pct", "best_coverage_pct", "frameshifts", "stop_codons",
                                     "contig", "genome_start", "genome_end"})
    clades = read_tsv(clades_path, {"assembly", "assigned_clade", "assignment_status"})
    manifest = read_tsv(manifest_path, {"assembly", "fasta"})
    annotation = maybe_tsv(metadata_path, {"accession", "protein_name"})
    meta_lookup = {r["accession"]: r for r in annotation}
    manifest_map = {r["assembly"]: r["fasta"] for r in manifest}
    if len(manifest_map) != len(manifest):
        raise ValueError("Duplicate assembly in genome_manifest.tsv")
    clade_map = {r["assembly"]: r for r in clades}
    if len(clade_map) != len(clades):
        raise ValueError("Duplicate assembly in clade_assignments.tsv")
    assemblies = set(manifest_map)
    if set(clade_map) != assemblies:
        raise ValueError(f"Clade/manifest mismatch. No clade row: {sorted(assemblies - set(clade_map))[:5]}; "
                         f"extra clades: {sorted(set(clade_map) - assemblies)[:5]}")
    counts = Counter(r["assembly"] for r in details)
    if set(counts) != assemblies:
        raise ValueError(f"Details/manifest mismatch. Missing: {sorted(assemblies - set(counts))[:5]}; "
                         f"extra: {sorted(set(counts) - assemblies)[:5]}")
    candidates = sorted({r["candidate"] for r in details})
    expected = len(candidates) * len(assemblies)
    if len(details) != expected:
        raise ValueError(f"Expected {expected} candidate×assembly rows, found {len(details)}")
    pairs = Counter((r["assembly"], r["candidate"]) for r in details)
    duplicate_pairs = [key for key, count in pairs.items() if count != 1]
    if duplicate_pairs:
        raise ValueError(f"Non-unique candidate×assembly combinations: {duplicate_pairs[:5]}")

    # Class assignments that remain unresolved are kept in a separate denominator.
    clade_counts = Counter(r["assigned_clade"].strip() if r["assigned_clade"].strip() in CLADES else "Unassigned"
                           for r in clades)
    if clade_counts["Unassigned"]:
        print(f"[WARN] {clade_counts['Unassigned']} genomes still Unassigned; not forced into a clade", file=sys.stderr)

    edge_lengths = {}
    if args.edge_check:
        for i, (assembly, filepath) in enumerate(manifest_map.items(), 1):
            p = Path(filepath)
            if not p.is_file():
                raise FileNotFoundError(f"Missing FASTA for {assembly}: {filepath}")
            edge_lengths[assembly] = fasta_lengths(p)
            if i % 20 == 0:
                print(f"[EDGE] Processed {i}/{len(manifest_map)} genomes", file=sys.stderr)

    # We audit existing classifications; do not overwrite miniprot summaries or clade calls.
    enriched = []
    by_candidate: dict[str, list] = defaultdict(list)
    by_candidate_clade: dict[tuple, list] = defaultdict(list)
    status_counter = Counter()
    issue_counter = Counter()
    edge_missing = Counter()
    for r in details:
        acc, candidate, status = r["assembly"], r["candidate"], r["status"]
        if status not in VALID_STATUSES:
            raise ValueError(f"Unknown miniprot status {status!r} in {acc} / {candidate}")
        nhit = integer(r["n_alignments"])
        nstrict = integer(r["n_strict_alignments"])
        if nhit is None or nstrict is None or nhit < nstrict or nhit < 0 or nstrict < 0:
            raise ValueError(f"Invalid alignment counts in {acc}/{candidate}")
        if (status == "no_alignment") != (nhit == 0):
            raise ValueError(f"Contradictory no_alignment/counts in {acc}/{candidate}")
        if status == "pass" and nstrict < 1:
            raise ValueError(f"Contradictory pass/strict count in {acc}/{candidate}")
        if status == "partial_or_divergent" and (nhit == 0 or nstrict != 0):
            raise ValueError(f"Contradictory partial-or-divergent counts in {acc}/{candidate}")
        identity = number(r["best_identity_pct"])
        coverage = number(r["best_coverage_pct"])
        frameshifts = integer(r["frameshifts"], 0)
        stops = integer(r["stop_codons"], 0)
        if nhit and (identity is None or coverage is None):
            raise ValueError(f"Missing identity/coverage with reported hit: {acc}/{candidate}")
        if identity is not None and not (0 <= identity <= 100) or coverage is not None and not (0 <= coverage <= 100.0001):
            raise ValueError(f"Invalid identity/coverage: {acc}/{candidate}")
        group = clade_map[acc]["assigned_clade"].strip()
        if group not in CLADES:
            group = "Unassigned"
        flags = []
        if status == "partial_or_divergent":
            if coverage is not None and coverage < min_coverage and identity is not None and identity >= min_identity:
                flags.append("partial_coverage")
            elif identity is not None and identity < min_identity and coverage is not None and coverage >= min_coverage:
                flags.append("low_identity")
            else:
                flags.append("partial_and_or_divergent")
        elif status == "no_alignment":
            flags.append("no_reported_alignment")
        if nstrict > 1:
            flags.append("multiple_strict_loci")
        elif nhit > 1:
            flags.append("multiple_reported_alignments")
        if (frameshifts or 0) > 0:
            flags.append("frameshift_in_best")
        if (stops or 0) > 0:
            flags.append("stop_codon_in_best")
        if status == "pass" and ((identity is not None and identity < min_identity + args.near_identity) or
                                 (coverage is not None and coverage < min_coverage + args.near_coverage)):
            flags.append("near_pass_threshold")
        if group == "Unassigned":
            flags.append("clade_unassigned")
        edge_distance = ""
        if args.edge_check and status != "no_alignment":
            contig = r["contig"]
            contig_len = edge_lengths[acc].get(contig)
            st = integer(r["genome_start"])
            en = integer(r["genome_end"])
            if contig_len is None:
                edge_missing[acc] += 1
                flags.append("contig_name_not_found")
            elif st is None or en is None:
                flags.append("missing_genome_coordinates")
            elif st < 1 or en > contig_len or st > en:
                flags.append("invalid_genome_coordinates")
            else:
                edge_distance = min(st - 1, contig_len - en)
                if edge_distance <= args.edge_window:
                    flags.append("near_contig_end")
        severity = "high" if any(x in flags for x in (
            "multiple_strict_loci", "frameshift_in_best", "stop_codon_in_best", "invalid_genome_coordinates")) else (
            "review" if flags else "ok")
        ann = meta_lookup.get(candidate, {})
        item = {
            "assembly": acc, "candidate": candidate,
            "protein_name": ann.get("protein_name", ""),
            "evidence_group": ann.get("evidence_group", ""),
            "clade": group, "clade_assignment_status": clade_map[acc]["assignment_status"],
            "status": status, "identity_pct": r["best_identity_pct"],
            "coverage_pct": r["best_coverage_pct"],
            "n_alignments": nhit, "n_strict_alignments": nstrict,
            "frameshifts_best": frameshifts if nhit else "",
            "stop_codons_best": stops if nhit else "",
            "contig": r["contig"], "genome_start": r["genome_start"], "genome_end": r["genome_end"],
            "distance_to_contig_end_bp": edge_distance,
            "flags": ";".join(flags), "severity": severity,
        }
        enriched.append(item)
        by_candidate[candidate].append(item)
        by_candidate_clade[candidate, group].append(item)
        status_counter[status] += 1
        issue_counter.update(flags)
    if len(enriched) != expected:
        raise AssertionError("Unexpected row count after enrichment")
    enriched.sort(key=lambda r: (r["candidate"], r["assembly"]))
    detail_fields = list(enriched[0])
    write_tsv(out / "mapping_audit.tsv", detail_fields, enriched)
    review = [r for r in enriched if r["severity"] != "ok"]
    review.sort(key=lambda r: ({"high": 0, "review": 1}.get(r["severity"], 2), r["candidate"], r["assembly"]))
    write_tsv(out / "review_queue.tsv", detail_fields, review)

    group_order = CLADES + ("Unassigned",)
    by_clade = []
    for candidate in candidates:
        name = meta_lookup.get(candidate, {}).get("protein_name", "")
        for clade in group_order:
            n = clade_counts[clade]
            records = by_candidate_clade[candidate, clade]
            if len(records) != n:
                raise ValueError(f"Candidate {candidate}/clade {clade}: found {len(records)} rows, expected {n}")
            counts = Counter(r["status"] for r in records)
            good = counts["pass"]
            lo, hi = wilson(good, n)
            by_clade.append({
                "candidate": candidate, "protein_name": name, "clade": clade, "n_genomes": n,
                "n_pass": good, "n_partial_or_divergent": counts["partial_or_divergent"],
                "n_no_alignment": counts["no_alignment"],
                "pass_pct": f"{100 * good / n:.2f}" if n else "NA",
                "wilson95_low_pct": lo, "wilson95_high_pct": hi,
                "interpretation_note": ("no_sample" if n == 0 else "low_n" if n < args.min_clade_n else
                                        "provisional_ANI_clades" if clade != "Unassigned" else "clade_unassigned"),
            })
    by_clade_fields = list(by_clade[0])
    write_tsv(out / "candidate_by_clade.tsv", by_clade_fields, by_clade)

    cand_rows = []
    for candidate in candidates:
        rows = by_candidate[candidate]
        statuses = Counter(r["status"] for r in rows)
        flags = Counter(flag for r in rows for flag in r["flags"].split(";") if flag)
        good, total = statuses["pass"], len(rows)
        lo, hi = wilson(good, total)
        cand_rows.append({
            "candidate": candidate,
            "protein_name": meta_lookup.get(candidate, {}).get("protein_name", ""),
            "evidence_group": meta_lookup.get(candidate, {}).get("evidence_group", ""),
            "n_genomes": total, "n_pass": good,
            "n_partial_or_divergent": statuses["partial_or_divergent"],
            "n_no_alignment": statuses["no_alignment"],
            "pass_pct": f"{100 * good / total:.2f}",
            "wilson95_low_pct": lo, "wilson95_high_pct": hi,
            "n_multiple_strict_loci": flags["multiple_strict_loci"],
            "n_frameshifts_best": flags["frameshift_in_best"],
            "n_stop_codons_best": flags["stop_codon_in_best"],
            "n_near_contig_end": flags["near_contig_end"],
            "n_high_severity": sum(r["severity"] == "high" for r in rows),
            "n_any_flags": sum(r["severity"] != "ok" for r in rows),
        })
    cand_rows.sort(key=lambda r: (-r["n_high_severity"], -r["n_partial_or_divergent"], r["candidate"]))
    write_tsv(out / "candidate_qc_summary.tsv", list(cand_rows[0]), cand_rows)

    clade_rows = []
    for clade in group_order:
        n = clade_counts[clade]
        clade_rows.append({"clade": clade, "n_genomes": n,
                           "sampling_flag": "no_genomes" if n == 0 else "low_n" if n < args.min_clade_n else "sampled",
                           "clade_status": "ANI_supported_provisional" if clade != "Unassigned" else "unassigned"})
    write_tsv(out / "clade_panel_qc.tsv", list(clade_rows[0]), clade_rows)

    report = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "project_root": str(root),
        "inputs": {"mapping_details": str(details_path), "clade_assignments": str(clades_path),
                   "manifest": str(manifest_path), "metadata_optional": str(metadata_path)},
        "assemblies": len(assemblies), "candidates": len(candidates), "candidate_assembly_pairs": len(enriched),
        "status_counts": dict(status_counter), "flag_counts": dict(issue_counter),
        "mapping_thresholds": {"identity_pct": min_identity, "candidate_coverage_pct": min_coverage},
        "clade_counts": {key: clade_counts[key] for key in group_order},
        "unassigned_assemblies": clade_counts["Unassigned"],
        "review_queue_rows": len(review), "high_severity_rows": sum(r["severity"] == "high" for r in enriched),
        "contig_edge_check_enabled": args.edge_check,
        "contig_names_missing_in_edge_check": dict(edge_missing),
        "warnings": (["ANI clade assignment is provisional, not phylogenetic confirmation."] +
                     ([f"One or more clades have fewer than {args.min_clade_n} assemblies; within-clade percentages are unstable."]
                      if any(0 < clade_counts[k] < args.min_clade_n for k in CLADES) else []) +
                     (["Some genomes lack assigned clades and have been kept separately."] if clade_counts["Unassigned"] else [])),
        "interpretation": (f"Pass is based on miniprot alignment thresholds ({min_identity:g}% identity, {min_coverage:g}% candidate coverage) "
                           "from prior step. These results do not demonstrate one-to-one orthology, gene expression, "
                           "surface localization or immunogenicity. Wilson intervals show sampling uncertainty "
                           "conditional on these non-random genomes and must not be treated as population prevalence CIs."),
    }
    out.mkdir(parents=True, exist_ok=True)
    (out / "audit_summary.json").write_text(json.dumps(report, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(f"[OK] {len(assemblies)} genomes × {len(candidates)} proteins = {len(enriched)} records")
    print(f"[OK] clades: {dict(report['clade_counts'])}")
    print(f"[OK] statuses: {dict(status_counter)}; review items: {len(review)}")
    print(f"[OK] output: {out}")
    if args.edge_check and edge_missing:
        print(f"[WARN] GFF contig IDs missing in FASTA for {len(edge_missing)} genomes", file=sys.stderr)


if __name__ == "__main__":
    try:
        main()
    except (ValueError, FileNotFoundError, OSError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        sys.exit(2)
