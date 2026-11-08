#!/usr/bin/env python3
"""Stage 4 of C. auris genomic survey: QC genomes and summarize miniprot GFF3.
Python >=3.9; standard library only. Does not infer clade membership or vaccine efficacy.
"""
from __future__ import annotations
import argparse, csv, hashlib, json, gzip, re, statistics, sys
from collections import Counter, defaultdict
from pathlib import Path

ACCESSION = re.compile(r'GC[AF]_\d+\.\d+')
DNA_EXTS = ('.fna', '.fna.gz', '.fa', '.fa.gz', '.fasta', '.fasta.gz')


def open_txt(p):
    return gzip.open(p, 'rt') if str(p).endswith('.gz') else open(p, 'rt', encoding='utf-8')


def write_tsv(path, header, rows):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open('w', newline='', encoding='utf-8') as fh:
        w = csv.writer(fh, delimiter='\t')
        w.writerow(header)
        w.writerows(rows)


def read_expect(path):
    with open_txt(path) as fh:
        first = fh.readline()
    if '\t' in first and first.lower().startswith('assembly accession'):
        with open_txt(path) as fh:
            rows = list(csv.DictReader(fh, delimiter='\t'))
        accessions = [r['Assembly Accession'].strip() for r in rows]
    else:
        with open_txt(path) as fh:
            accessions = [line.strip().split()[0] for line in fh if line.strip()]
    if not accessions or any(not ACCESSION.fullmatch(a) for a in accessions):
        raise ValueError('Expected --selected file must contain GCA/GCF accessions')
    if len(accessions) != len(set(accessions)):
        raise ValueError('Duplicate accessions inside --selected')
    return accessions


def extract_acc(path):
    for part in reversed(path.parts):
        match = ACCESSION.search(part)
        if match:
            return match.group(0)
    return None


def genome_stats(path):
    nseq = 0; lengths = []; nbases = 0; other = 0; gc = 0; total = 0
    digest_list = []
    seq_parts = []
    def flush():
        nonlocal nseq, total, nbases, other, gc, seq_parts
        if not seq_parts: return
        s = ''.join(seq_parts).upper()
        lengths.append(len(s)); total += len(s); nseq += 1
        nbases += s.count('N'); gc += s.count('G') + s.count('C')
        other += sum(1 for c in s if c not in 'ACGTN')
        digest_list.append(hashlib.sha256(s.encode('ascii', errors='replace')).hexdigest())
        seq_parts = []
    with open_txt(path) as fh:
        for line in fh:
            if line.startswith('>'):
                flush()
            else:
                seq_parts.append(line.strip())
        flush()
    if total == 0: raise ValueError(f'Empty or invalid FASTA: {path}')
    lengths.sort(reverse=True)
    running = 0; n50 = 0
    for length in lengths:
        running += length
        if running >= total / 2:
            n50 = length; break
    digest = hashlib.sha256('\n'.join(sorted(digest_list)).encode()).hexdigest()
    return [nseq, total, n50, max(lengths), round(nbases*100/total, 4), round(gc*100/total, 4), other, digest]


def prepare(args):
    selected = read_expect(args.selected)
    expected = set(selected)
    root = Path(args.genomes).resolve(); out = Path(args.out).resolve()
    out.mkdir(parents=True, exist_ok=True)
    if not root.is_dir(): raise ValueError(f'Genome directory not found: {root}')
    all_genomes = [p for p in root.rglob('*') if p.is_file() and p.name.lower().endswith(DNA_EXTS)]
    buckets = defaultdict(list)
    for p in all_genomes:
        acc = extract_acc(p)
        if acc in expected: buckets[acc].append(p)
    issues = []; qc = []; manifest=[]; seen_hash = {}
    for acc in selected:
        choices = buckets.get(acc, [])
        if not choices:
            issues.append([acc, 'ERROR', 'Missing assembly FASTA']); continue
        if len(choices) > 1:
            issues.append([acc, 'ERROR', 'Multiple FASTA files for one assembly; resolve manually']); continue
        path = choices[0]
        try:
            stats = genome_stats(path)
        except Exception as exc:
            issues.append([acc, 'ERROR', str(exc)]); continue
        nseq, length, n50, longest, n_percent, gc_percent, other, digest = stats
        if length < 11800000 or length > 13500000: issues.append([acc, 'WARNING', f'Unexpected genome length: {length}'])
        if nseq > 50: issues.append([acc, 'WARNING', f'High fragmentation: {nseq} sequences'])
        if n_percent > 1: issues.append([acc, 'WARNING', f'High ambiguous fraction: {n_percent}% N'])
        if other: issues.append([acc, 'WARNING', f'Non-ACGTN symbols: {other}'])
        if digest in seen_hash: issues.append([acc, 'WARNING', f'Exact duplicate assembled sequences of {seen_hash[digest]}'])
        else: seen_hash[digest] = acc
        manifest.append([acc, str(path)])
        qc.append([acc, str(path), nseq, length, n50, longest, n_percent, gc_percent, other, digest])
    write_tsv(out/'genome_manifest.tsv', ['assembly','fasta'], manifest)
    write_tsv(out/'genome_qc.tsv', ['assembly','fasta','n_contigs','length_bp','N50_bp','longest_contig_bp','N_pct','GC_pct','non_ACGTN','sequence_hash'], qc)
    write_tsv(out/'genome_issues.tsv', ['assembly','severity','message'], issues)
    (out/'genome_paths.txt').write_text(''.join(row[1] + '\n' for row in manifest))
    (out/'genome_missing.txt').write_text(''.join(acc+'\n' for acc in selected if acc not in {r[0] for r in manifest}))
    counts = Counter(x[1] for x in issues)
    summary = {'expected_assemblies':len(selected), 'valid_genome_paths':len(manifest), 'n_errors':counts['ERROR'], 'n_warnings':counts['WARNING'], 'duplicate_sequence_assemblies':sum('Exact duplicate' in x[2] for x in issues)}
    (out/'qc_summary.json').write_text(json.dumps(summary, indent=2, ensure_ascii=False)+'\n')
    print(json.dumps(summary, ensure_ascii=False, indent=2))
    if summary['n_errors']: sys.exit(2)


def fasta_ids(path):
    names=[]
    with open_txt(path) as fh:
        for line in fh:
            if line.startswith('>'):
                names.append(line[1:].strip().split()[0])
    if len(names) != len(set(names)): raise ValueError('Duplicate candidate FASTA identifiers')
    return names


def parse_gff(path):
    hits=[]
    with open_txt(path) as fh:
        for line in fh:
            if not line or line[0] == '#': continue
            fields=line.rstrip('\n').split('\t')
            if len(fields)<9 or fields[2]!='mRNA': continue
            attrs={}
            for kv in fields[8].split(';'):
                if '=' in kv:
                    k,v=kv.split('=',1); attrs[k]=v
            target=attrs.get('Target','').split()
            if len(target)<3: continue
            q=target[0]; start=int(target[1]); end=int(target[2]);
            ident=float(attrs.get('Identity','0')) * 100
            frameshifts=int(attrs.get('Frameshift','0'))
            stop=int(attrs.get('StopCodon','0'))
            hits.append({'candidate':q,'contig':fields[0], 'start_genome':int(fields[3]),'end_genome':int(fields[4]), 'strand':fields[6], 'start_protein':start,'end_protein':end, 'identity_pct':ident, 'frameshifts':frameshifts, 'stop_codons':stop})
    return hits


def faa_lengths(path):
    seq_lengths={}; cur=None
    with open_txt(path) as fh:
        for line in fh:
            if line.startswith('>'):
                cur=line[1:].strip().split()[0]; seq_lengths[cur]=0
            elif cur:
                seq_lengths[cur]+=len(line.strip())
    return seq_lengths


def summarize(args):
    out=Path(args.out); out.mkdir(parents=True, exist_ok=True)
    with open(args.manifest, encoding='utf-8') as fh:
        manifest=list(csv.DictReader(fh, delimiter='\t'))
    ids = faa_lengths(args.candidates)
    if not ids: raise ValueError('No candidate sequences')
    if len(ids)!=len(set(ids)): raise ValueError('Duplicate candidate IDs')
    gffdir=Path(args.gff_dir)
    if not manifest: raise ValueError('Empty manifest')
    long=[]; matrix=[]; counts=Counter(); missing_files=[]
    assembly_order=[r['assembly'] for r in manifest]
    for rec in manifest:
        acc=rec['assembly']
        gff=gffdir/(acc+'.gff')
        done=gffdir/(acc + args.stamp_suffix)
        if not gff.is_file() or not done.is_file():
            missing_files.append(acc); continue
        grouped=defaultdict(list)
        for hit in parse_gff(gff):
            if hit['candidate'] not in ids: continue
            hit['coverage_pct']=(hit['end_protein']-hit['start_protein']+1)*100/ids[hit['candidate']]
            grouped[hit['candidate']].append(hit)
        for cid in ids:
            hs=grouped[cid]
            good=[h for h in hs if h['identity_pct'] >= args.min_identity and h['coverage_pct'] >= args.min_coverage]
            # This is best aligned match, not ortholog assignment.
            hs.sort(key=lambda h:(h['coverage_pct']*h['identity_pct'], h['coverage_pct']), reverse=True)
            best=hs[0] if hs else None
            status='pass' if good else ('partial_or_divergent' if hs else 'no_alignment')
            counts[status]+=1
            row = [acc,cid,status,len(hs),len(good)]
            if best:
                row += [best['identity_pct'],round(best['coverage_pct'],2),best['contig'],best['start_genome'],best['end_genome'],best['strand'],best['frameshifts'],best['stop_codons']]
            else: row += ['']*8
            long.append(row)
    if missing_files: raise RuntimeError('Missing/unfinished miniprot output: '+', '.join(missing_files[:12]) + (f' (+{len(missing_files)-12})' if len(missing_files)>12 else ''))
    cols=['assembly','candidate','status','n_alignments','n_strict_alignments','best_identity_pct','best_coverage_pct','contig','genome_start','genome_end','strand','frameshifts','stop_codons']
    write_tsv(out/'candidate_genome_details.tsv',cols,long)
    by_pair={(r[1],r[0]): r[2] for r in long}
    for cid in ids:
        matrix.append([cid]+[{'pass':1,'partial_or_divergent':0,'no_alignment':0}[by_pair[(cid,a)]] for a in assembly_order])
    write_tsv(out/'candidate_presence_matrix.tsv',['candidate']+assembly_order,matrix)
    summary=[]
    for cid in ids:
        statuses=[by_pair[(cid,a)] for a in assembly_order]
        summary.append([cid,len(assembly_order),statuses.count('pass'),statuses.count('partial_or_divergent'),statuses.count('no_alignment'),round(100*statuses.count('pass')/len(statuses),2)])
    summary.sort(key=lambda r:(-r[5],r[0]))
    write_tsv(out/'candidate_presence_summary.tsv',['candidate','n_genomes','n_pass','n_partial_or_divergent','n_no_alignment','pass_pct'], summary)
    meta={'genomes':len(assembly_order),'candidates':len(ids),'filters':{'identity_pct':args.min_identity,'protein_coverage_pct':args.min_coverage},'mapping_status_counts':dict(counts), 'interpretation':'Protein-to-genome matches; not orthology, expression, accessibility, or vaccine efficacy'}
    (out/'summary.json').write_text(json.dumps(meta,indent=2,ensure_ascii=False)+'\n')
    print(json.dumps(meta,indent=2,ensure_ascii=False))


def main():
    p=argparse.ArgumentParser(description=__doc__)
    sub=p.add_subparsers(dest='cmd',required=True)
    q=sub.add_parser('qc',help='Validate selected genomes, make skani and miniprot input manifest')
    q.add_argument('--genomes',required=True,help='Unpacked NCBI data directory, with accession-specific subfolders')
    q.add_argument('--selected',required=True,help='selected_100.tsv or selected_accessions.txt')
    q.add_argument('--out',default='14_population/qc')
    m=sub.add_parser('summarize',help='Summarize miniprot --gff results')
    m.add_argument('--manifest',required=True)
    m.add_argument('--gff-dir',required=True)
    m.add_argument('--stamp-suffix', default='.done', choices=('.done', '.done.json'))
    m.add_argument('--candidates',required=True)
    m.add_argument('--min-identity',type=float,default=40.0)
    m.add_argument('--min-coverage',type=float,default=70.0)
    m.add_argument('--out',default='14_population/results')
    a=p.parse_args()
    (prepare if a.cmd=='qc' else summarize)(a)

if __name__=='__main__': main()
