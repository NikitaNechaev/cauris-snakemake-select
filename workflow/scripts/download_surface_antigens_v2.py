#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Download UniProtKB-annotated fungal surface/extracellular protein candidates.

Examples:
    python download_surface_antigens_v2.py --scope auris --out results
    python download_surface_antigens_v2.py --scope family --out results
    python download_surface_antigens_v2.py --scope expanded --out results
    python download_surface_antigens_v2.py --scope expanded --max-pages 1 --out test

Needs: requests (pip install requests)
Note: UniProt annotations provide candidate evidence, not verified immunogenicity.
"""

import argparse
import csv
import io
import json
import re
import sys
import time
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path
from urllib.parse import urlparse

import requests

BASE_URL = "https://rest.uniprot.org/uniprotkb/search"

TAXA = {
    "auris": "organism_id:498019",
    "family": "taxonomy_id:27319",
    # Historical Candida sensu lato is polyphyletic. Explicit species are more reliable.
    "expanded": "(taxonomy_id:27319 OR organism_id:5476 OR organism_id:5482 "
                "OR organism_id:5480 OR organism_id:5478)",
}

# Fetch narrowly annotated candidate classes, NOT complete proteomes.
# Gene Ontology terms: GO:0009277 fungal-type cell wall; GO:0005618 cell wall;
# GO:0009986 cell surface; GO:0005576 extracellular region.
CATEGORIES = {
    "annot_cell_wall": 'cc_scl_term:"Cell wall"',
    "go_fungal_cell_wall": "go:0009277",
    "go_cell_wall": "go:0005618",
    "go_cell_surface": "go:0009986",
    "name_adhesin": "protein_name:adhesin",
    "name_agglutinin": "protein_name:agglutinin",
    "name_mannoprotein": "protein_name:mannoprotein",
    "name_surface": "protein_name:surface",
    "name_gpi": "protein_name:GPI",
    "go_extracellular": "go:0005576",
}

FIELDS = [
    "accession", "protein_name", "gene_names", "organism_name", "organism_id",
    "reviewed", "length", "cc_subcellular_location", "ft_lipid",
    "ft_signal", "go_c", "keyword", "sequence",
]

OUTPUT_FIELDS = [
    "accession", "protein_name", "gene_names", "organism_name", "organism_id",
    "reviewed", "length", "cc_subcellular_location", "ft_lipid",
    "ft_signal", "go_c", "keyword", "matched_categories", "evidence_group",
    "sequence_length", "sequence_sha256",
]


def retry_get(session, url, params=None, timeout=120, retries=5):
    """Get UniProt page; retry transient failures, never hide HTTP 400."""
    last = None
    for attempt in range(retries):
        try:
            response = session.get(url, params=params, timeout=timeout)
            if response.status_code in {429, 500, 502, 503, 504}:
                delay = min(40, 2 ** (attempt + 1))
                if response.headers.get("Retry-After", "").isdigit():
                    delay = min(90, int(response.headers["Retry-After"]))
                print(f"  HTTP {response.status_code}; repeat in {delay}s "
                      f"({attempt + 1}/{retries})", file=sys.stderr)
                last = f"HTTP {response.status_code}: {response.text[:500]}"
                if attempt < retries - 1:
                    time.sleep(delay)
                    continue
            response.raise_for_status()
            return response
        except (requests.Timeout, requests.ConnectionError) as error:
            last = repr(error)
            if attempt < retries - 1:
                time.sleep(min(40, 2 ** (attempt + 1)))
                continue
            raise RuntimeError(f"Network error after {retries} attempts: {last}") from error
        except requests.HTTPError as error:
            raise RuntimeError(f"HTTP {response.status_code}: "
                               f"{response.text[:1000]}") from error
    raise RuntimeError(f"Unable to fetch data: {last}")


def pages(session, query, max_pages=0, size=500):
    url = BASE_URL
    params = {
        "query": query,
        "format": "tsv",
        "fields": ",".join(FIELDS),
        "size": size,
    }
    number = 0
    while url:
        response = retry_get(session, url, params=params)
        content_type = response.headers.get("Content-Type", "")
        if "html" in content_type.lower():
            raise ValueError("Received HTML instead of TSV")
        reader = csv.reader(io.StringIO(response.text), delimiter="\t")
        header = next(reader, None)
        if header is None or len(header) != len(FIELDS):
            raise ValueError(f"Unexpected TSV columns ({len(header or [])}): {header!r}; "
                             f"expected {len(FIELDS)}. Response: {response.text[:300]}")
        chunk = []
        for values in reader:
            if not values or all(not cell for cell in values):
                continue
            if len(values) != len(FIELDS):
                raise ValueError(f"Bad TSV row with {len(values)} columns")
            row = dict(zip(FIELDS, values))
            if row["accession"]:
                chunk.append(row)
        yield chunk, response.headers.get("x-total-results", ""), response.headers.get("x-uniprot-release", "")
        number += 1
        if max_pages and number >= max_pages:
            break
        next_link = response.links.get("next", {}).get("url")
        if next_link and urlparse(next_link).netloc != "rest.uniprot.org":
            raise ValueError(f"Unexpected pagination host: {next_link}")
        url, params = next_link, None


def evidence_group(record):
    """Descriptive annotation bins, NOT immunogenicity or localization proof."""
    name = record.get("protein_name", "").lower()
    loc = record.get("cc_subcellular_location", "").lower()
    go = record.get("go_c", "").lower()
    lipid = record.get("ft_lipid", "").lower()
    keywords = record.get("keyword", "").lower()
    matched = record.get("matched_categories", set())

    if ("cell wall" in loc or "fungal-type cell wall" in go
        or "go:0009277" in go or "go:0005618" in go
        or "annot_cell_wall" in matched or "go_fungal_cell_wall" in matched
        or "go_cell_wall" in matched):
        return "cell_wall_annotated"
    if ("gpi-anchor" in lipid or "gpi anchor" in lipid or
        "adhesin" in name or "agglutinin" in name or
        "mannoprotein" in name or "cell surface" in loc or
        "go:0009986" in go):
        return "surface_related_annotation"
    if ("extracellular" in loc or "secreted" in loc or
        "go:0005576" in go or "go_extracellular" in matched):
        return "extracellular_related_annotation"
    return "name_or_other_annotation_only"


def safe_fasta_description(record):
    desc = f"{record['organism_name']} | {record['protein_name']}"
    return re.sub(r"\s+", " ", desc).strip()[:180]


def write_fasta(records, path):
    with path.open("w", encoding="utf-8") as out:
        for accession, record in sorted(records.items()):
            seq = record["sequence"].strip().replace(" ", "")
            if not seq:
                continue
            out.write(f">{accession} {safe_fasta_description(record)}\n")
            for i in range(0, len(seq), 80):
                out.write(seq[i:i + 80] + "\n")


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--scope", choices=TAXA, default="expanded",
                        help="auris, family (Metschnikowiaceae), expanded (+ clinical Candida spp.)")
    parser.add_argument("--out", type=Path, default=Path("09_surface_antigens"))
    parser.add_argument("--max-pages", type=int, default=0,
                        help="Limit pages per category for smoke test; 0 = all pages")
    parser.add_argument("--sleep", type=float, default=0.3)
    args = parser.parse_args()
    if args.max_pages < 0:
        parser.error("--max-pages must be non-negative")

    args.out.mkdir(parents=True, exist_ok=True)
    session = requests.Session()
    session.headers.update({"User-Agent": "Academic-surface-protein-catalog/1.0 (UniProtKB REST API)",
                            "Accept": "text/tab-separated-values"})

    records = {}
    audit = []
    release = ""
    errors = 0
    scope_query = TAXA[args.scope]
    print("Scope:", args.scope, "=>", scope_query, flush=True)
    for category, filter_query in CATEGORIES.items():
        query = f"({scope_query}) AND ({filter_query})"
        print(f"\n[SEARCH] {category}\n  {query}", flush=True)
        count = 0
        chunks = 0
        status = "OK"
        reported_total = ""
        try:
            for page, total, page_release in pages(session, query, max_pages=args.max_pages):
                if page_release:
                    if release and release != page_release:
                        raise ValueError(f"UniProt release changed during the download: {release} → {page_release}")
                    release = page_release
                reported_total = total or reported_total
                chunks += 1
                count += len(page)
                for entry in page:
                    key = entry["accession"]
                    if key not in records:
                        records[key] = entry
                        records[key]["matched_categories"] = set()
                    elif (records[key].get("sequence") and entry["sequence"] and
                          records[key]["sequence"] != entry["sequence"]):
                        raise ValueError(f"Different sequences for accession {key}")
                    records[key]["matched_categories"].add(category)
                print(f"  page {chunks}: +{len(page)}; cumulative matches: {count}", flush=True)
                time.sleep(args.sleep)
        except (requests.RequestException, RuntimeError, ValueError) as error:
            status = f"FAILED: {error}"
            errors += 1
            print(f"ERROR [{category}]: {error}", file=sys.stderr)
        audit.append({"category": category, "query": query,
                      "reported_total": reported_total,
                      "retrieved_rows": count, "pages": chunks, "status": status})

    import hashlib
    for record in records.values():
        record["evidence_group"] = evidence_group(record)
        record["sequence_length"] = len(record["sequence"].strip())
        record["sequence_sha256"] = hashlib.sha256(record["sequence"].encode()).hexdigest()

    with (args.out / "surface_candidates.tsv").open("w", encoding="utf-8", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=OUTPUT_FIELDS, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        for accession, record in sorted(records.items()):
            row = dict(record)
            row["matched_categories"] = ";".join(sorted(record["matched_categories"]))
            writer.writerow(row)

    write_fasta(records, args.out / "surface_candidates.faa")
    wall = {acc: r for acc, r in records.items() if r["evidence_group"] == "cell_wall_annotated"}
    write_fasta(wall, args.out / "cell_wall_annotated.faa")

    with (args.out / "query_audit.tsv").open("w", encoding="utf-8", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=["category", "query", "reported_total", "retrieved_rows", "pages", "status"], delimiter="\t")
        writer.writeheader()
        writer.writerows(audit)

    org_counts = Counter((r["organism_id"], r["organism_name"]) for r in records.values())
    with (args.out / "by_organism.tsv").open("w", encoding="utf-8", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t")
        writer.writerow(["organism_id", "organism_name", "unique_accessions"])
        for (taxid, name), n in sorted(org_counts.items(), key=lambda kv: (-kv[1], kv[0][1])):
            writer.writerow([taxid, name, n])

    metadata = {
        "timestamp_utc": datetime.now(timezone.utc).isoformat(),
        "uniprot_release": release or None,
        "scope": args.scope,
        "scope_query": scope_query,
        "max_pages_per_category": args.max_pages,
        "total_unique_accessions": len(records),
        "cell_wall_annotated": len(wall),
        "categories_failed": errors,
        "evidence_groups": dict(Counter(r["evidence_group"] for r in records.values())),
        "queries": audit,
    }
    (args.out / "run_metadata.json").write_text(json.dumps(metadata, indent=2, ensure_ascii=False), encoding="utf-8")

    print("\nRESULTS:", args.out.resolve())
    print("  Unique accessions:", len(records))
    print("  Wall-annotated:", len(wall))
    print("  Organism taxa:", len(org_counts))
    print("  Failed categories:", errors)
    if args.max_pages:
        print("NOTE: TEST MODE. Pagination was intentionally truncated.", file=sys.stderr)
    if errors:
        print("WARNING: Catalog is INCOMPLETE. Inspect query_audit.tsv!", file=sys.stderr)
        sys.exit(2)
    if not records:
        print("WARNING: Zero records; verify taxonomy and search filters.", file=sys.stderr)
        sys.exit(3)


if __name__ == "__main__":
    main()
