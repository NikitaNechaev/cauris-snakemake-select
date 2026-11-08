#!/usr/bin/env python3
"""Compare a fixed protein panel with annotated C. auris assemblies.

Workflow
--------
  1) `datasets download genome taxon 498019 --annotated --assembly-source RefSeq --include protein ...`
  2) `unzip ...`
  3) `python check_auris_isolates.py prepare --ncbi-dir 13_isolates/ncbi --out 13_isolates`
  4) Create DIAMOND database from the 51 query proteins and align combined annotated proteins to it.
  5) `python check_auris_isolates.py summarize --hits 13_isolates/diamond_hits.tsv --candidates ...`

Outputs describe *assembly-level similarity*, not biological demonstration of
antigen conservation, orthology, surface accessibility or vaccine efficacy.
"""
import argparse
import csv
import json
import re
import sys
from collections import Counter, defaultdict
from pathlib import Path

FIELDNAMES = ["qseqid","sseqid","pident","length","qlen","slen","qstart","qend","sstart","send","evalue","bitscore"]
ASSEMBLY_RE = re.compile(r"^GC[AF]_\d+\.\d+$")
AA_ALPHABET = set("ACDEFGHIKLMNPQRSTVWYXUBZOJ*-")


def read_fasta(path):
    identifier = None
    seq_parts = []
    with open(path, encoding="utf-8") as f:
        for line in f:
            if line.startswith(">"):
                if identifier is not None:
                    yield identifier, "".join(seq_parts)
                identifier = line[1:].strip().split()[0]
                seq_parts = []
            elif identifier is not None:
                seq_parts.append(line.strip())
        if identifier is not None:
            yield identifier, "".join(seq_parts)


def write_tsv(path, cols, rows):
    with open(path, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=cols, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def prepare(args):
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    proteomes = sorted(Path(args.ncbi_dir).rglob("protein.faa"))
    proteomes = [p for p in proteomes if ASSEMBLY_RE.match(p.parent.name)]
    if not proteomes:
        raise SystemExit("No accession/protein.faa files found in extracted NCBI Datasets directory")
    if args.max_assemblies:
        proteomes = proteomes[:args.max_assemblies]
        print("WARNING: limited assembly set. FOR SMOKE TEST ONLY; not for final conclusions.",file=sys.stderr)
    count_rows=[]
    seen_assemblies=set()
    combined=out / "combined_proteomes.faa"
    n_total=0
    with combined.open("w",encoding="utf-8") as f:
        for path in proteomes:
            accession = path.parent.name
            if accession in seen_assemblies:
                raise SystemExit(f"Duplicated assembly folder: {accession}")
            seen_assemblies.add(accession)
            seen_proteins=set()
            n_proteins=0
            n_empty=0
            for protein_id,sequence in read_fasta(path):
                if protein_id in seen_proteins:
                    raise SystemExit(f"Duplicate protein ID {protein_id} inside {accession}")
                seen_proteins.add(protein_id)
                sequence=sequence.upper()
                if not sequence:
                    n_empty += 1
                    continue
                bad=set(sequence)-AA_ALPHABET
                if bad:
                    raise SystemExit(f"Unexpected letters {bad} in {accession} {protein_id}")
                f.write(f">{accession}|{protein_id}\n")
                for start in range(0, len(sequence), 80):
                    f.write(sequence[start:start+80]+"\n")
                n_proteins += 1
            n_total += n_proteins
            count_rows.append({"assembly":accession,"n_proteins":n_proteins,"n_empty":n_empty,"source_file":str(path)})
    write_tsv(out/"assembly_counts.tsv",["assembly","n_proteins","n_empty","source_file"],count_rows)
    with (out/"prepare_stats.json").open("w",encoding="utf-8") as f:
        json.dump({"assemblies":len(count_rows),"proteins":n_total,"test_limited":bool(args.max_assemblies)},f,indent=2)
    print(f"Prepared {len(count_rows)} annotated assemblies, {n_total} proteins")
    print(f"FASTA: {combined}")


def summarize(args):
    out = Path(args.out)
    out.mkdir(parents=True,exist_ok=True)
    all_candidates=[k for k,v in read_fasta(args.candidates) if v]
    if len(all_candidates)!=len(set(all_candidates)):
        raise SystemExit("Duplicate candidate FASTA headers")
    if not all_candidates:
        raise SystemExit("Candidate FASTA empty")
    candidates=set(all_candidates)
    assemblies=[]
    with open(args.counts,newline="",encoding="utf-8") as f:
        for r in csv.DictReader(f,delimiter="\t"):
            if int(r['n_proteins'])>0:
                assemblies.append(r['assembly'])
    assemblies=set(assemblies)
    if not assemblies:
        raise SystemExit("No nonempty assemblies in counts file")
    best={}
    n_rows=0
    missing_ids=Counter()
    with open(args.hits,encoding="utf-8") as f:
        reader=csv.reader(f,delimiter="\t")
        for idx,row in enumerate(reader,1):
            if not row:
                continue
            if len(row)!=len(FIELDNAMES):
                raise SystemExit(f"Invalid DIAMOND format at line {idx}, expected {len(FIELDNAMES)} columns, got {len(row)}")
            d=dict(zip(FIELDNAMES,row))
            n_rows+=1
            qid=d['qseqid']
            if '|' not in qid:
                raise SystemExit(f"Expected prefix 'GCF_xxx|proteinID' in query ID: {qid}")
            asm,prot=qid.split('|',1)
            if asm not in assemblies:
                missing_ids[asm]+=1
                continue
            if d['sseqid'] not in candidates:
                raise SystemExit(f"Unknown candidate ID: {d['sseqid']}")
            try:
                v = {k:float(d[k]) for k in ('pident','qlen','slen','qstart','qend','sstart','send','evalue','bitscore')}
            except ValueError as e:
                raise SystemExit(f"Non-numeric DIAMOND output on line {idx}") from e
            qlen,slen=v['qlen'],v['slen']
            if qlen < 1 or slen < 1:
                raise SystemExit(f"Invalid sequence length on line {idx}")
            qcov=100*(abs(v['qend']-v['qstart'])+1)/qlen
            scov=100*(abs(v['send']-v['sstart'])+1)/slen
            key=(asm,d['sseqid'])
            entry={"assembly":asm,"candidate":d['sseqid'],"protein_id":prot,
                   "identity_pct":round(v['pident'],2),
                   "protein_coverage_pct":round(qcov,2),
                   "candidate_coverage_pct":round(scov,2),
                   "bitscore":v['bitscore'],"evalue":v['evalue']}
            rank=(v['bitscore'],min(qcov,scov),qcov+scov)
            if key not in best or rank>best[key][0]:
                best[key]=(rank,entry)
    if missing_ids:
        raise SystemExit(f"Hits from unknown assemblies: {dict(missing_ids)}")
    b_rows=[]
    tallies=Counter()
    for asm,cand in sorted(best):
        entry=best[(asm,cand)][1]
        entry['strict_hit']=int(entry['identity_pct']>=args.min_id and
                                entry['protein_coverage_pct']>=args.min_protein_cov and
                                entry['candidate_coverage_pct']>=args.min_candidate_cov)
        tallies[(cand,'any')]+=1
        tallies[(cand,'strict')]+=entry['strict_hit']
        b_rows.append(entry)
    bcols=['assembly','candidate','protein_id','identity_pct','protein_coverage_pct','candidate_coverage_pct','bitscore','evalue','strict_hit']
    write_tsv(out/'best_hits_by_assembly.tsv',bcols,b_rows)
    summary=[]
    matrix=[]
    for cand in all_candidates:
        strict=tallies[(cand,'strict')]
        summary.append({'candidate':cand,'n_assemblies_total':len(assemblies),
                        'n_assemblies_any_hit':tallies[(cand,'any')],
                        'n_assemblies_strict_hit':strict,
                        'strict_pct':round(100*strict/len(assemblies),2)})
    write_tsv(out/'candidate_presence_summary.tsv',
              ['candidate','n_assemblies_total','n_assemblies_any_hit','n_assemblies_strict_hit','strict_pct'],summary)
    lookup={(r['assembly'],r['candidate']):r['strict_hit'] for r in b_rows}
    with (out/'candidate_presence_matrix.tsv').open('w',encoding='utf-8',newline='') as f:
        w=csv.writer(f,delimiter='\t')
        w.writerow(['assembly']+all_candidates)
        for asm in sorted(assemblies):
            w.writerow([asm]+[lookup.get((asm,c),0) for c in all_candidates])
    with (out/'search_stats.json').open('w',encoding='utf-8') as f:
        json.dump({"assemblies":len(assemblies),"candidates":len(all_candidates),
                   "diamond_rows":n_rows,"candidate_assembly_hits":len(best),
                   "strict_thresholds":{"identity_pct":args.min_id,
                                        "protein_coverage_pct":args.min_protein_cov,
                                        "candidate_coverage_pct":args.min_candidate_cov}},f,indent=2)
    print(f"Processed {n_rows} DIAMOND rows across {len(assemblies)} assemblies")
    print('Saved: candidate_presence_summary.tsv, candidate_presence_matrix.tsv, best_hits_by_assembly.tsv')
    print('NOTE: zeros include missing/poor annotation, not necessarily biological absence.')


def main():
    ap=argparse.ArgumentParser(description=__doc__,formatter_class=argparse.RawDescriptionHelpFormatter)
    subs=ap.add_subparsers(dest='command',required=True)
    p=subs.add_parser('prepare',help='Create assembly-labelled protein FASTA from NCBI Datasets zip extract')
    p.add_argument('--ncbi-dir',required=True)
    p.add_argument('--out',required=True)
    p.add_argument('--max-assemblies',type=int,default=0,help='Only for technical testing')
    p.set_defaults(fn=prepare)
    s=subs.add_parser('summarize',help='Summarize DIAMOND blastp vs candidate proteins')
    s.add_argument('--hits',required=True)
    s.add_argument('--candidates',required=True)
    s.add_argument('--counts',required=True)
    s.add_argument('--out',required=True)
    s.add_argument('--min-id',type=float,default=40)
    s.add_argument('--min-protein-cov',type=float,default=70)
    s.add_argument('--min-candidate-cov',type=float,default=70)
    s.set_defaults(fn=summarize)
    args=ap.parse_args()
    args.fn(args)

if __name__=='__main__':
    main()
