#!/usr/bin/env python3
"""TEN C. auris: audit precomputed skani ANI+AF matrices against NCBI assembly manifest.

This analyzes similarity and putative near-clone components; it does NOT infer
published C. auris clades or a phylogenetic tree.

Dependencies: Python 3.9+, numpy, matplotlib.
Input:  skani triangle matrices (PHYLIP full or triangular), QC manifest TSV.
Output: numeric pair list, connected components, QC JSON and PNG heatmap.
"""
from __future__ import annotations

import argparse
import csv
import json
import math
import re
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np

ACC_RE = re.compile(r"GC[AF]_\d+\.\d+")


def read_manifest(path: Path):
    with path.open(encoding='utf-8', newline='') as fh:
        r = csv.DictReader(fh, delimiter='\t')
        if not r.fieldnames or not {'assembly', 'fasta'}.issubset(r.fieldnames):
            raise ValueError('QC manifest needs columns assembly and fasta')
        rows = list(r)
    out = {}
    basename_map = defaultdict(list)
    for row in rows:
        acc = row['assembly'].strip()
        name = row['fasta'].strip()
        if not ACC_RE.fullmatch(acc) or acc in out:
            raise ValueError(f'Invalid or duplicate assembly in manifest: {acc}')
        out[acc] = name
        basename_map[Path(name).name].append(acc)
    if not out:
        raise ValueError('Empty manifest')
    return out, basename_map


def parse_matrix(path: Path, directional: bool = False):
    """Read skani PHYLIP matrix. ANI is symmetric; AF is directional.

    A full AF matrix has A->B and B->A fractions relative to different
    genome lengths, so they need not be equal. Zero means missing/skipped.
    """
    with path.open(encoding='utf-8') as fh:
        lines = [line.strip() for line in fh if line.strip() and not line.lstrip().startswith('#')]
    if not lines:
        raise ValueError(f'Empty skani matrix: {path}')
    first = lines[0].split()
    if len(first) != 1 or not first[0].isdigit():
        raise ValueError(f'Expected PHYLIP count in first row of {path}, got {lines[0][:80]!r}')
    n = int(first[0])
    if n < 2 or len(lines) != n + 1:
        raise ValueError(f'Matrix {path}: declared {n} rows, found {len(lines)-1}')
    names = []
    raw = []
    for line in lines[1:]:
        # skani uses tabs, so a FASTA path may itself contain spaces.
        parts = line.split("\t") if "\t" in line else line.split()
        names.append(parts[0])
        try:
            raw.append([float(v) for v in parts[1:]])
        except ValueError as e:
            raise ValueError(f'Non-numeric value in {path}, {parts[0]}') from e
    if len(set(names)) != n:
        raise ValueError(f'Duplicate labels in {path}')
    lengths = [len(x) for x in raw]
    a = np.full((n, n), np.nan, float)
    if all(x == n for x in lengths):
        for i in range(n):
            a[i, :] = raw[i]
        # ANI can use a value from the opposite triangle if one is blank.
        # Do NOT do this to AF: one direction is not the other direction.
        if not directional:
            for i in range(n):
                for j in range(i+1, n):
                    x, y = a[i, j], a[j, i]
                    if (not np.isfinite(x) or x == 0) and np.isfinite(y) and y > 0:
                        a[i, j] = y
                    if (not np.isfinite(y) or y == 0) and np.isfinite(x) and x > 0:
                        a[j, i] = x
    elif directional:
        raise ValueError(f"AF requires skani --full-matrix to preserve both directional coverage fractions: {path}")
    elif lengths == [i+1 for i in range(n)]:
        for i in range(n):
            a[i, :i+1] = raw[i]
            for j in range(i):
                a[j, i] = raw[i][j]
    elif lengths == [i for i in range(n)]:
        for i in range(n):
            a[i, :i] = raw[i]
            for j in range(i):
                a[j, i] = raw[i][j]
    elif lengths == [n-i for i in range(n)]:
        for i in range(n):
            a[i, i:] = raw[i]
            for j in range(i+1,n):
                a[j, i] = raw[i][j-i]
    elif lengths == [n-i-1 for i in range(n)]:
        for i in range(n):
            a[i, i+1:] = raw[i]
            for j in range(i+1,n):
                a[j, i] = raw[i][j-i-1]
    else:
        raise ValueError(f'Unrecognized matrix layout in {path}; first widths {lengths[:8]}')
    if np.any(np.isfinite(a) & (a < 0)):
        raise ValueError(f'Negative matrix entry: {path}')
    for i in range(n):
        a[i, i] = np.nan
    # 0 may mean skani did not compute this ANI; treat it as missing.
    a[a == 0] = np.nan
    if directional:
        return names, a
    for i in range(n):
        for j in range(i+1,n):
            x, y = a[i,j], a[j,i]
            if np.isfinite(x) and np.isfinite(y):
                if abs(x-y) > 0.05:
                    raise ValueError(f'Non-symmetric matrix {path}, {names[i]} and {names[j]}: {x}, {y}')
                a[i,j] = a[j,i] = (x+y)/2
            elif np.isfinite(x):
                a[j,i] = x
            elif np.isfinite(y):
                a[i,j] = y
    return names, a


def index_labels(names, manifest, basename_map):
    mapped = []
    inverse_paths = {str(Path(p).resolve()): acc for acc, p in manifest.items()}
    for s in names:
        acc = None
        match = ACC_RE.search(s)
        if match and match.group(0) in manifest:
            acc = match.group(0)
        else:
            maybe = basename_map.get(Path(s).name, [])
            if len(maybe) == 1:
                acc = maybe[0]
            else:
                acc = inverse_paths.get(str(Path(s).resolve()))
        if acc is None:
            raise ValueError(f'Matrix label cannot be matched to manifest: {s}')
        mapped.append(acc)
    if len(set(mapped)) != len(mapped):
        raise ValueError('Multiple matrix labels resolve to the same assembly')
    missing_in_ani = sorted(set(manifest)-set(mapped))
    extra_in_ani = sorted(set(mapped)-set(manifest))
    if missing_in_ani or extra_in_ani:
        raise ValueError(f'Matrix/manifest mismatch. Missing in matrix: {missing_in_ani[:8]}; extra: {extra_in_ani[:8]}')
    return mapped


def reorder(matrix, labels, target_labels):
    index = {name: i for i,name in enumerate(labels)}
    order = [index[x] for x in target_labels]
    return matrix[np.ix_(order,order)]


class UnionFind:
    def __init__(self, n):
        self.parent = list(range(n))
        self.size = [1]*n
    def find(self, x):
        while self.parent[x] != x:
            self.parent[x] = self.parent[self.parent[x]]
            x = self.parent[x]
        return x
    def union(self, x, y):
        a, b = self.find(x), self.find(y)
        if a != b:
            if self.size[a] < self.size[b]: a,b=b,a
            self.parent[b] = a
            self.size[a] += self.size[b]


def components(ani, af, cutoff, afcut, accessions):
    uf = UnionFind(len(accessions))
    n = len(accessions)
    for i in range(n):
        for j in range(i+1,n):
            if (np.isfinite(ani[i,j]) and ani[i,j] >= cutoff and
                np.isfinite(af[i,j]) and af[i,j] >= afcut):
                uf.union(i,j)
    groups = defaultdict(list)
    for i, acc in enumerate(accessions):
        groups[uf.find(i)].append(i)
    groups_sorted = sorted(groups.values(), key=lambda members: (-len(members), min(accessions[i] for i in members)))
    output = {}
    for k, member in enumerate(groups_sorted,1):
        for i in member:
            output[i] = (f'G{k:03}',len(member))
    return output, groups_sorted


def emit_tsv(path, headers, rows):
    path.parent.mkdir(parents=True,exist_ok=True)
    with path.open('w', newline='', encoding='utf-8') as out:
        writer = csv.writer(out, delimiter='\t')
        writer.writerow(headers)
        writer.writerows(rows)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    ap.add_argument('--root', default=None, help='TEN root (normally inferred from scripts/ location)')
    ap.add_argument('--ani', default=None)
    ap.add_argument('--af', default=None)
    ap.add_argument('--manifest',default=None)
    ap.add_argument('--out',default=None)
    ap.add_argument('--clone-ani',type=float,default=99.9,help='Near-clone ANI cutoff; heuristic, not formal clade')
    ap.add_argument('--group-ani',type=float,default=99.0,help='Exploratory groups cutoff; not phylogenetic clades')
    ap.add_argument('--min-af',type=float,default=70.0,help='Minimum aligned fraction percent cutoff in BOTH directions')
    ap.add_argument('--no-plot',action='store_true')
    args=ap.parse_args()
    root = Path(args.root).resolve() if args.root else Path(__file__).resolve().parent.parent
    pop = root/'08_population'
    ani_path = Path(args.ani) if args.ani else pop/'skani_matrix.txt'
    af_path = Path(args.af) if args.af else pop/'skani_matrix.txt.af'
    manifest_path = Path(args.manifest) if args.manifest else pop/'qc'/'genome_manifest.tsv'
    out = Path(args.out) if args.out else pop/'results'/'ani_population'
    for p in (ani_path, af_path, manifest_path):
        if not p.is_file():
            raise SystemExit(f'ERROR missing required file: {p}')
    if not 0 <= args.min_af <= 100 or not 0 <= args.group_ani <= 100 or not 0 <= args.clone_ani <= 100:
        raise SystemExit('Cutoffs should be percentages between 0 and 100')
    manifest, by_base = read_manifest(manifest_path)
    ani_names, ani = parse_matrix(ani_path)
    af_names, af = parse_matrix(af_path, directional=True)
    acc = index_labels(ani_names, manifest, by_base)
    af_acc = index_labels(af_names, manifest, by_base)
    if set(acc) != set(af_acc):
        raise ValueError('ANI and AF assembly sets differ')
    af = reorder(af, af_acc, acc)
    # Avoid incorrectly interpreting a fractional AF matrix as 0..100.
    positive_af = af[np.isfinite(af) & (af > 0)]
    scale = 100.0 if len(positive_af) and float(np.nanmax(positive_af)) <= 1.0 else 1.0
    af = af*scale
    # Conservative pairwise coverage: both genomes must be sufficiently covered.
    # np.minimum yields NaN if either direction is missing.
    af_min = np.minimum(af, af.T)
    upper = np.triu_indices(len(acc),1)
    finite = np.isfinite(ani[upper])
    if len(acc) != len(manifest):
        raise ValueError('ANI matrix does not contain exactly QC genomes')
    for title, m in [('ANI',ani),('AF',af)]:
        bad=np.isfinite(m) & ((m < 0) | (m>100.000001))
        if np.any(bad):
            raise ValueError(f'{title} matrix has values outside 0-100')
    near, near_components = components(ani,af_min,args.clone_ani,args.min_af,acc)
    group, group_components = components(ani,af_min,args.group_ani,args.min_af,acc)
    out.mkdir(parents=True,exist_ok=True)
    near_rows=[]
    per_genome=[]
    for i,a in enumerate(acc):
        hits=[(ani[i,j],j,af_min[i,j]) for j in range(len(acc)) if j!=i and np.isfinite(ani[i,j]) and np.isfinite(af_min[i,j])]
        hits.sort(reverse=True)
        nearest = hits[0] if hits else None
        n_near = sum(1 for ai,j,afi in hits if ai>=args.clone_ani and afi>=args.min_af)
        per_genome.append([a, manifest[a],near[i][0],near[i][1],group[i][0],group[i][1],
                           acc[nearest[1]] if nearest else '',
                           f'{nearest[0]:.5f}' if nearest else '',
                           f'{nearest[2]:.3f}' if nearest else '',n_near])
        for j in range(i+1,len(acc)):
            if (np.isfinite(ani[i,j]) and np.isfinite(af_min[i,j]) and
                ani[i,j]>=args.clone_ani and af_min[i,j]>=args.min_af):
                near_rows.append([a,acc[j],f'{ani[i,j]:.5f}',f'{af_min[i,j]:.3f}',
                                  f'{af[i,j]:.3f}',f'{af[j,i]:.3f}'])
    emit_tsv(out/'genome_ani_groups.tsv',
             ['assembly','genome_fasta','near_clone_component','near_clone_component_size',
              'exploratory_component','exploratory_component_size','nearest_neighbor','nearest_ANI_pct',
              'nearest_AF_pct','near_clone_neighbors'],per_genome)
    emit_tsv(out/'near_clone_pairs.tsv', ['assembly_1','assembly_2','ANI_pct','AF_min_pct',
                                         'AF_1_to_2_pct','AF_2_to_1_pct'], near_rows)
    all_pairs=[]
    for i in range(len(acc)):
        for j in range(i+1,len(acc)):
            if np.isfinite(ani[i,j]):
                all_pairs.append([acc[i],acc[j],f'{ani[i,j]:.5f}',
                                  f'{af_min[i,j]:.3f}' if np.isfinite(af_min[i,j]) else '',
                                  f'{af[i,j]:.3f}' if np.isfinite(af[i,j]) else '',
                                  f'{af[j,i]:.3f}' if np.isfinite(af[j,i]) else ''])
    emit_tsv(out/'all_estimated_pairs.tsv',
             ['assembly_1','assembly_2','ANI_pct','AF_min_pct','AF_1_to_2_pct','AF_2_to_1_pct'],all_pairs)
    data=ani[upper]; valid=data[np.isfinite(data)]
    summary={
      'n_genomes':len(acc),
      'n_pairs_possible':len(data),
      'n_pairs_estimated':int(finite.sum()),
      'n_pairs_missing_or_skipped':int((~finite).sum()),
      'ani_min_pct': float(valid.min()) if valid.size else None,
      'ani_median_pct': float(np.median(valid)) if valid.size else None,
      'ani_max_pct': float(valid.max()) if valid.size else None,
      'near_clone_ani_cutoff_pct':args.clone_ani,
      'exploratory_ani_cutoff_pct':args.group_ani,
      'min_aligned_fraction_pct':args.min_af,
      'af_values_scaled_from_fraction':bool(scale==100.0),
      'af_pair_policy':'Both directions must meet min-AF; use min(AF_1_to_2,AF_2_to_1).',
      'n_pairs_with_both_af_directions':int(np.isfinite(af_min[upper]).sum()),
      'n_pairs_with_asymmetric_af_gt_0_05_pp':int(np.sum(
          np.isfinite(af[upper]) & np.isfinite(af.T[upper]) &
          (np.abs(af[upper]-af.T[upper])>0.05))),
      'n_near_clone_components':len(near_components),
      'n_exploratory_components':len(group_components),
      'n_near_clone_pairs':len(near_rows),
      'warning':'AF is directional, so it is not symmetrized. Connected components are threshold networks, not clades. ANI 0/missing is not a biological zero. Published clades need references/phylogenetics.'
    }
    (out/'ani_qc_summary.json').write_text(json.dumps(summary,indent=2,ensure_ascii=False)+'\n')
    if not args.no_plot:
        try:
            import matplotlib
            matplotlib.use('Agg')
            import matplotlib.pyplot as plt
            # Group near-clones together; missing comparisons are displayed in grey.
            order = sorted(range(len(acc)), key=lambda i: (group[i][0],near[i][0],acc[i]))
            view = np.ma.masked_invalid(ani[np.ix_(order,order)].copy())
            for i in range(len(order)):
                view[i,i]=100.0
            colors=plt.get_cmap('viridis').copy()
            colors.set_bad('#d0d0d0')
            fig,ax=plt.subplots(figsize=(9,8))
            im=ax.imshow(view,cmap=colors,interpolation='nearest',vmin=90,vmax=100)
            ax.set_title(f'skani ANI: {len(acc)} C. auris assemblies (cluster-grouped)')
            ax.set_xlabel('Assemblies, ordered by threshold groups')
            ax.set_ylabel('Assemblies, ordered by threshold groups')
            ax.set_xticks([]);ax.set_yticks([])
            cb=fig.colorbar(im,ax=ax,fraction=0.045,pad=0.03)
            cb.set_label('Estimated ANI (%)  |  grey = missing/skipped')
            fig.tight_layout()
            fig.savefig(out/'ani_heatmap.png',dpi=180)
            plt.close(fig)
        except ImportError:
            print('WARN: matplotlib is not installed; heatmap skipped.')
    print(json.dumps(summary,indent=2,ensure_ascii=False))
    print(f'Wrote results to: {out}')

if __name__ == '__main__':
    main()
