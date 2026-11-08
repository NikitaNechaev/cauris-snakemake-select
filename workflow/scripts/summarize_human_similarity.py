#!/usr/bin/env python3
"""Summarize DIAMOND blastp against human reference proteome.

Input hit format:
--outfmt 6 qseqid sseqid pident length qlen slen qstart qend sstart send evalue bitscore

Outputs every query accession, including proteins with no reported hit.
No-hit means no alignment passed the chosen search thresholds, not proof
of absence of similarity, immune cross-reactivity or vaccine safety.
"""
import argparse
import csv
from collections import Counter
from pathlib import Path

FIELDS = ['qseqid', 'sseqid', 'pident', 'length', 'qlen', 'slen',
          'qstart','qend','sstart','send','evalue','bitscore']


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--hits', type=Path, required=True)
    ap.add_argument('--metadata', type=Path, required=True)
    ap.add_argument('--out', type=Path, required=True)
    args = ap.parse_args()
    with args.metadata.open(newline='', encoding='utf-8') as f:
        metas = list(csv.DictReader(f, delimiter='\t'))
    by_acc = {r['accession']: r for r in metas}
    if len(by_acc) != len(metas):
        raise ValueError('Duplicate accession in input metadata')
    hits_count = Counter()
    best = {}
    with args.hits.open(encoding='utf-8', newline='') as f:
        reader = csv.reader(f, delimiter='\t')
        for idx, t in enumerate(reader, start=1):
            if not t:
                continue
            if len(t) != len(FIELDS):
                raise ValueError(f'Wrong DIAMOND columns in row {idx}: {len(t)} (expected 12)')
            d = dict(zip(FIELDS, t))
            q = d['qseqid']
            if q not in by_acc:
                raise ValueError(f'Unknown query ID in DIAMOND output: {q}')
            hits_count[q] += 1
            try:
                d['pident'] = float(d['pident'])
                d['bitscore'] = float(d['bitscore'])
                d['evalue'] = float(d['evalue'])
                d['qlen'] = int(d['qlen'])
                d['qstart'] = int(d['qstart'])
                d['qend'] = int(d['qend'])
                if d['qlen'] < 1:
                    raise ValueError('Zero qlen')
            except ValueError as ex:
                raise ValueError(f'Invalid numerical field in hit row {idx}') from ex
            d['qcov_pct'] = round(100 * (abs(d['qend']-d['qstart'])+1) / d['qlen'], 2)
            if q not in best or (d['bitscore'],d['qcov_pct']) > (best[q]['bitscore'],best[q]['qcov_pct']):
                best[q] = d
    cols = ['accession','protein_name','reviewed','evidence_group','cluster_id',
            'n_reported_hsps','has_reported_hit','best_human_hit',
            'best_pident_pct','best_query_coverage_pct','best_evalue','best_bitscore']
    args.out.parent.mkdir(parents=True, exist_ok=True)
    with args.out.open('w', encoding='utf-8', newline='') as f:
        w=csv.DictWriter(f, fieldnames=cols, delimiter='\t')
        w.writeheader()
        for acc,m in sorted(by_acc.items()):
            h=best.get(acc)
            w.writerow({
                'accession':acc,
                'protein_name':m.get('protein_name',''),
                'reviewed':m.get('reviewed',''),
                'evidence_group':m.get('evidence_group',''),
                'cluster_id':m.get('cluster_id',''),
                'n_reported_hsps':hits_count.get(acc,0),
                'has_reported_hit':str(h is not None),
                'best_human_hit':h['sseqid'] if h else '',
                'best_pident_pct':h['pident'] if h else '',
                'best_query_coverage_pct':h['qcov_pct'] if h else '',
                'best_evalue':h['evalue'] if h else '',
                'best_bitscore':h['bitscore'] if h else '',
            })
    print(f'Queries: {len(by_acc)}; with hits: {len(best)}; without reported hits: {len(by_acc)-len(best)}')
    print(f'Summary: {args.out}')

if __name__ == '__main__':
    main()
