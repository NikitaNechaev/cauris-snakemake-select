#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Summarize MMseqs2 easy-cluster assignments with UniProtKB metadata.

Usage:
 python summarize_surface_clusters.py \
   --clusters 11_clusters/surface_cluster.tsv \
   --metadata 10_qc/surface_priority.tsv \
   --out 11_clusters/cluster_summary.tsv

Important: clusters are similarity groups, NOT proven orthologous groups.
"""
import argparse
import csv
from collections import defaultdict, Counter
from pathlib import Path


def approximate_species(name):
    parts = name.split()
    return ' '.join(parts[:2]) if len(parts) >= 2 else name


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--clusters', required=True, type=Path)
    ap.add_argument('--metadata', required=True, type=Path)
    ap.add_argument('--out', required=True, type=Path)
    args = ap.parse_args()

    with args.metadata.open(newline='', encoding='utf-8') as f:
        metadata = {r['accession']: r for r in csv.DictReader(f, delimiter='\t')}

    clusters = defaultdict(set)
    with args.clusters.open(encoding='utf-8') as f:
        for line_no, line in enumerate(f, 1):
            if not line.strip():
                continue
            pair = line.rstrip('\n').split('\t')
            if len(pair) < 2:
                raise ValueError(f'Bad cluster row {line_no}: {line!r}')
            representative, member = pair[0], pair[1]
            clusters[representative].add(member)

    rows, missing = [], Counter()
    for rep, members in clusters.items():
        known = [metadata[a] for a in members if a in metadata]
        for a in members:
            if a not in metadata:
                missing[a] += 1
        species = sorted(set(approximate_species(r['organism_name']) for r in known))
        contains_auris = any(' auris' in (' ' + r['organism_name'].lower()) for r in known)
        gene_names = sorted(set(r.get('protein_name', '') for r in known if r.get('protein_name')))
        n_reviewed = sum('reviewed' == r.get('reviewed', '').strip().lower() for r in known)
        rows.append({
            'representative': rep,
            'n_members': len(members),
            'n_species_labels_approx': len(species),
            'has_auris': str(contains_auris),
            'n_reviewed': n_reviewed,
            'species_labels_approx': '; '.join(species),
            'example_protein_names': '; '.join(gene_names[:5]),
            'member_accessions': ';'.join(sorted(members)),
        })

    rows.sort(key=lambda r: (-int(r['has_auris'] == 'True'),
                             -r['n_species_labels_approx'], -r['n_members'], r['representative']))
    args.out.parent.mkdir(parents=True, exist_ok=True)
    fields = ['representative', 'n_members', 'n_species_labels_approx',
              'has_auris', 'n_reviewed', 'species_labels_approx',
              'example_protein_names', 'member_accessions']
    with args.out.open('w', encoding='utf-8', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=fields, delimiter='\t')
        writer.writeheader()
        writer.writerows(rows)

    print('Clusters:', len(rows))
    print('Clusters containing C. auris:', sum(r['has_auris'] == 'True' for r in rows))
    print('Missing metadata accessions:', len(missing))
    if missing:
        print('Examples:', ', '.join(list(missing)[:10]))
    print('Report:', args.out)


if __name__ == '__main__':
    main()
