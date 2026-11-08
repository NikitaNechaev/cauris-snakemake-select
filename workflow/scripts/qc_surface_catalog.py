#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""QC and classify a UniProtKB fungal surface-candidate catalog downloaded by
 download_surface_antigens_v2.py.

Usage:
  python qc_surface_catalog.py --input 09_surface_antigens --out 10_qc

No third-party Python dependencies.
Produces a descriptive catalog only: evidence_group is not evidence of
immunogenicity or vaccine efficacy.
"""
import argparse
import csv
import hashlib
import json
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path

REQUIRED_COLS = {
    'accession', 'organism_name', 'organism_id', 'protein_name',
    'reviewed', 'length', 'evidence_group', 'matched_categories'
}
PRIMARY_GROUPS = {'cell_wall_annotated', 'surface_related_annotation'}
AA_STANDARD = set('ACDEFGHIKLMNPQRSTVWY')


def iter_fasta(path):
    header, parts = None, []
    with path.open(encoding='utf-8') as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith('>'):
                if header is not None:
                    yield header.split()[0], ''.join(parts)
                header, parts = line[1:], []
            else:
                if header is None:
                    raise ValueError(f'{path}: sequence before FASTA header')
                parts.append(line)
    if header is not None:
        yield header.split()[0], ''.join(parts)


def write_tsv(path, rows, columns):
    with path.open('w', newline='', encoding='utf-8') as fh:
        writer = csv.DictWriter(fh, fieldnames=columns, delimiter='\t', extrasaction='ignore')
        writer.writeheader()
        writer.writerows(rows)


def write_fasta(path, items, seqs):
    with path.open('w', encoding='utf-8') as fh:
        for accession in sorted(items):
            seq = seqs.get(accession)
            if not seq:
                continue
            fh.write('>' + accession + '\n')
            for i in range(0, len(seq), 80):
                fh.write(seq[i:i + 80] + '\n')


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--input', type=Path, required=True)
    ap.add_argument('--out', type=Path, default=Path('10_qc'))
    args = ap.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)

    src = args.input
    tsv = src / 'surface_candidates.tsv'
    fasta = src / 'surface_candidates.faa'
    if not tsv.is_file() or not fasta.is_file():
        ap.error(f'Expected both {tsv} and {fasta}')

    issues = []
    def note(severity, issue, accession='', detail=''):
        issues.append(dict(severity=severity, issue=issue,
                           accession=accession, detail=str(detail)))

    metadata_path = src / 'run_metadata.json'
    runmeta = {}
    if metadata_path.exists():
        runmeta = json.loads(metadata_path.read_text(encoding='utf-8'))
        if runmeta.get('categories_failed', 0) != 0:
            note('ERROR', 'download_categories_failed', detail=runmeta['categories_failed'])
        if runmeta.get('max_pages_per_category', 0) != 0:
            note('WARNING', 'download_was_page_limited', detail=runmeta['max_pages_per_category'])
    else:
        note('WARNING', 'missing_download_metadata')

    audit_path = src / 'query_audit.tsv'
    if audit_path.exists():
        with audit_path.open(encoding='utf-8', newline='') as fh:
            audit = list(csv.DictReader(fh, delimiter='\t'))
        for q in audit:
            if q.get('status', '').strip() != 'OK':
                note('ERROR', 'query_failed', detail=f"{q.get('category')}: {q.get('status')}")
    else:
        note('WARNING', 'missing_query_audit')

    records = {}
    with tsv.open(encoding='utf-8', newline='') as fh:
        reader = csv.DictReader(fh, delimiter='\t')
        columns = reader.fieldnames or []
        missing = REQUIRED_COLS - set(columns)
        if missing:
            raise ValueError('Missing expected TSV columns: ' + ', '.join(sorted(missing)))
        for row in reader:
            accession = row['accession'].strip()
            if not accession:
                note('ERROR', 'empty_accession')
                continue
            if accession in records:
                note('ERROR', 'duplicate_accession_in_tsv', accession)
            records[accession] = row

    seqs = {}
    for accession, seq in iter_fasta(fasta):
        if accession in seqs:
            note('ERROR', 'duplicate_accession_in_fasta', accession)
        seqs[accession] = seq.upper()
        if not seq:
            note('ERROR', 'empty_sequence', accession)
        nonstandard = set(seq.upper()) - AA_STANDARD
        if nonstandard:
            note('WARNING', 'nonstandard_aa', accession, ''.join(sorted(nonstandard)))

    for accession in records.keys() - seqs.keys():
        note('ERROR', 'tsv_accession_missing_from_fasta', accession)
    for accession in seqs.keys() - records.keys():
        note('ERROR', 'fasta_accession_missing_from_tsv', accession)

    seq_to_acc = defaultdict(list)
    for accession, seq in seqs.items():
        if seq:
            seq_to_acc[hashlib.sha256(seq.encode()).hexdigest()].append(accession)
        row = records.get(accession)
        if row:
            if row.get('length', '').isdigit() and int(row['length']) != len(seq):
                note('WARNING', 'uniprot_length_mismatch', accession,
                     f"tsv={row['length']} fasta={len(seq)}")
            if row.get('sequence_length', '').isdigit() and int(row['sequence_length']) != len(seq):
                note('WARNING', 'download_length_mismatch', accession)
            if row.get('sequence_sha256') and row['sequence_sha256'] != hashlib.sha256(seq.encode()).hexdigest():
                note('WARNING', 'sha256_mismatch', accession)

    for _, accs in seq_to_acc.items():
        if len(accs) > 1:
            note('INFO', 'identical_aa_sequence_multiple_accessions',
                 ','.join(accs[:5]), f'n={len(accs)}')

    by_species = Counter()
    by_evidence = Counter()
    by_reviewed = Counter()
    species_evidence = Counter()
    for row in records.values():
        species = row.get('organism_name', '') or 'UNKNOWN'
        group = row.get('evidence_group', '') or 'UNKNOWN'
        by_species[species] += 1
        by_evidence[group] += 1
        by_reviewed[row.get('reviewed', '') or 'UNKNOWN'] += 1
        species_evidence[(species, group)] += 1

    write_tsv(args.out / 'counts_by_organism.tsv',
              [dict(organism_name=k, count=v) for k, v in by_species.most_common()],
              ['organism_name', 'count'])
    write_tsv(args.out / 'counts_by_evidence.tsv',
              [dict(evidence_group=k, count=v) for k, v in by_evidence.most_common()],
              ['evidence_group', 'count'])
    write_tsv(args.out / 'species_by_evidence.tsv',
              [dict(organism_name=sp, evidence_group=ev, count=n)
               for (sp, ev), n in sorted(species_evidence.items())],
              ['organism_name', 'evidence_group', 'count'])
    write_tsv(args.out / 'qc_issues.tsv', issues,
              ['severity', 'issue', 'accession', 'detail'])

    # Annotations only: priority dataset is not proof of external accessibility.
    priority = {acc for acc, r in records.items()
                if r.get('evidence_group') in PRIMARY_GROUPS and acc in seqs}
    write_fasta(args.out / 'surface_priority.faa', priority, seqs)
    priority_fields = list(columns)
    write_tsv(args.out / 'surface_priority.tsv',
              [records[a] for a in sorted(priority)], priority_fields)

    auris = {acc for acc, r in records.items()
             if ' auris' in (' ' + r.get('organism_name', '').lower()) and acc in seqs}
    write_fasta(args.out / 'cauris_all.faa', auris, seqs)
    write_tsv(args.out / 'cauris_all.tsv',
              [records[a] for a in sorted(auris)], priority_fields)

    result = {
        'timestamp_utc': datetime.now(timezone.utc).isoformat(),
        'records_tsv': len(records), 'records_fasta': len(seqs),
        'unique_sequence_hashes': len(seq_to_acc),
        'distinct_organism_labels': len(by_species),
        'auris_records': len(auris),
        'priority_records': len(priority),
        'evidence_groups': dict(by_evidence),
        'reviewed_values': dict(by_reviewed),
        'issue_severity_counts': dict(Counter(i['severity'] for i in issues)),
        'source_scope': runmeta.get('scope', 'UNKNOWN'),
    }
    (args.out / 'qc_summary.json').write_text(json.dumps(result, indent=2, ensure_ascii=False), encoding='utf-8')
    print(json.dumps(result, ensure_ascii=False, indent=2))
    print('\nGenerated:', args.out.resolve())
    print('Read qc_issues.tsv and query_audit.tsv before biological interpretation.')


if __name__ == '__main__':
    main()
