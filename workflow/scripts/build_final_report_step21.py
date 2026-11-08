#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""TEN / Lastovaz: step 21, evidence integration and transparent descriptive report.

Inputs are read from existing TEN stage directories; they are never modified.
Output: 12_final_report (tables, QC, figures, Russian Markdown report, ZIP).

This script does not infer vaccine efficacy/safety or design antigens/epitopes.
Requires Python >= 3.9; matplotlib/numpy for PNG figures (optional via --no-figures).

Usage:
    python scripts/build_final_report_step21.py
    python scripts/build_final_report_step21.py --check-only
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import re
import sys
import statistics
import zipfile
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path

VERSION = '1.0.0'
CLADES = ['I', 'II', 'III', 'IV', 'V', 'VI']

INPUTS = {
    'annotations': '06_blastp_on_human/auris_cluster_audit.tsv',
    'human': '06_blastp_on_human/auris_human_summary.tsv',
    'mapping': '08_population/results/mapping_audit/candidate_qc_summary.tsv',
    'mapping_by_clade': '08_population/results/mapping_audit/candidate_by_clade.tsv',
    'clade_assignments': '08_population/results/clade_typing/clade_assignments.tsv',
    'genome_manifest': '08_population/qc/genome_manifest.tsv',
    'sequence_inventory': '09_protein_sequences/sequence_inventory.tsv',
    'mafft_run': '10_mafft/mafft_summary.json',
    'conservation': '11_conservation_statistics/summary_by_candidate.tsv',
    'conservation_by_clade': '11_conservation_statistics/summary_by_clade.tsv',
    'conservation_run': '11_conservation_statistics/conservation_summary.json',
}
OPTIONAL = {
    'raw_input_provenance': '00_provenance/inputs_sha256.tsv',
    'workflow_configuration': '00_provenance/config.json',
    'uniprot_source_metadata': '01_surface_antigens/data/run_metadata.json',
    'initial_qc': '02_qc_catalog/qc_summary.json',
    'extraction_run': '09_protein_sequences/extraction_summary.json',
    'ani_qc': '08_population/results/ani_population/ani_qc_summary.json',
    'clade_run': '08_population/results/clade_typing/clade_summary.json',
    'mapping_audit_run': '08_population/results/mapping_audit/audit_summary.json',
}


def load_tsv(path: Path, required):
    with path.open('r', encoding='utf-8-sig', newline='') as f:
        dr = csv.DictReader(f, delimiter='\t')
        fields = dr.fieldnames or []
        absent = set(required) - set(fields)
        if absent:
            raise ValueError(f'{path}: нет столбцов {sorted(absent)}; имеются {fields}')
        rows = list(dr)
    if not rows:
        raise ValueError(f'{path}: пустая таблица')
    return rows


def load_json(path):
    return json.loads(path.read_text('utf-8'))


def index_unique(rows, key, label):
    found = {}
    for r in rows:
        k = r[key].strip()
        if not k or k in found:
            raise ValueError(f'{label}: пустой или повторяющийся {key}: {k!r}')
        found[k] = r
    return found


def number(value):
    if value is None or str(value).strip().lower() in ('', 'na', 'nan', 'none', 'null'):
        return None
    result = float(value)
    if not math.isfinite(result):
        raise ValueError(f'Некорректное числовое значение: {value!r}')
    return result


def integer(value):
    n = number(value)
    if n is None:
        return None
    if int(n) != n:
        raise ValueError(f'Ожидалось целое значение, получено {value!r}')
    return int(n)


def percent(value):
    n = number(value)
    if n is not None and not 0 <= n <= 100.00001:
        raise ValueError(f'Процент вне 0..100: {value!r}')
    return n


def yes(value):
    return str(value).strip().lower() in ('true', '1', 'yes', 'y')


def maybe(v, digits=2):
    return '' if v is None else round(v, digits)


def write_tsv(path, rows):
    path.parent.mkdir(parents=True, exist_ok=True)
    if not rows:
        raise ValueError(f'Пустые строки для {path}')
    cols = list(rows[0])
    with path.open('w', encoding='utf-8', newline='') as f:
        w = csv.DictWriter(f, fieldnames=cols, delimiter='\t', extrasaction='ignore')
        w.writeheader()
        w.writerows(rows)


def sha256(path):
    digest = hashlib.sha256()
    with path.open('rb') as f:
        for chunk in iter(lambda: f.read(1 << 20), b''):
            digest.update(chunk)
    return digest.hexdigest()


def as_md(value, decimals=1):
    return 'н/д' if value is None else f'{value:.{decimals}f}'


def chart_figures(out, candidates, byclade, clade_counts):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    import numpy as np
    plt.rcParams['font.family'] = 'DejaVu Sans'

    figures = out / 'figures'
    figures.mkdir(parents=True, exist_ok=True)
    created = []
    # Never silently depict no samples as 0% coverage.
    byclade_map = {(x['candidate'], x['clade']): x for x in byclade}
    subset = sorted(candidates, key=lambda r: (r['accession']))
    heat = np.full((len(subset), len(CLADES)), np.nan)
    for i, c in enumerate(subset):
        for j, cl in enumerate(CLADES):
            r = byclade_map[(c['accession'], cl)]
            if integer(r['n_genomes_in_clade']) > 0:
                heat[i,j] = percent(r['extraction_pct'])
    fig, ax = plt.subplots(figsize=(11, max(10, .28*len(subset)+2.5)))
    cmap = plt.get_cmap('viridis').copy()
    cmap.set_bad('#d8d8d8')
    im = ax.imshow(np.ma.masked_invalid(heat), vmin=0, vmax=100, aspect='auto', cmap=cmap)
    ax.set_yticks(range(len(subset)), [r['accession'] for r in subset], fontsize=7)
    ax.set_xticks(range(len(CLADES)), [f'{c}\nn={clade_counts[c]}' for c in CLADES])
    ax.set_xlabel('Предварительно назначенная клада (по ANI)')
    ax.set_title('Доля сборок с QC-извлечённым гомологом по кладам')
    fig.colorbar(im, ax=ax, fraction=.026, pad=.02, label='% сборок с извлечённым гомологом')
    fig.tight_layout()
    p = figures / '01_extraction_by_clade.png'
    fig.savefig(p, dpi=200, bbox_inches='tight'); plt.close(fig); created.append(p)

    xs, ys, names = [], [], []
    for c in candidates:
        x, y = c['extracted_pct'], c['mean_entropy_bits']
        if x is not None and y is not None:
            xs.append(x); ys.append(y); names.append(c['accession'])
    fig, ax = plt.subplots(figsize=(8.8, 6))
    ax.scatter(xs,ys,alpha=.72,s=42)
    ax.set(xlim=(-3, 103), xlabel='Сборки с QC-извлечённым гомологом (%)',
           ylabel='Средняя энтропия Шеннона (бит/позицию)',
           title='Представленность и вариабельность (описательная статистика)')
    for k in range(min(6,len(names))):
        ax.annotate(names[k],(xs[k],ys[k]),fontsize=7,xytext=(4,4),textcoords='offset points')
    ax.grid(alpha=.18)
    fig.tight_layout(); p = figures / '02_recovery_vs_entropy.png'
    fig.savefig(p, dpi=200); plt.close(fig); created.append(p)

    fig, ax = plt.subplots(figsize=(8.5, 4.8))
    count = Counter('Есть совпадение' if c['human_has_reported_hit'] else 'Нет сообщённых совпадений' for c in candidates)
    labels = ['Есть совпадение','Нет сообщённых совпадений']
    ax.bar(labels,[count[x] for x in labels])
    ax.set(ylabel='Количество белковых записей',title='DIAMOND: сравнение с протеомом человека')
    for i, k in enumerate(labels):
        ax.text(i,count[k]+.15,str(count[k]),ha='center')
    ax.text(.5,-.23,'Нет совпадения ≠ доказана иммунологическая безопасность',
            ha='center',transform=ax.transAxes,fontsize=9)
    fig.tight_layout(); p = figures / '03_human_reported_hits.png'
    fig.savefig(p,dpi=200,bbox_inches='tight'); plt.close(fig); created.append(p)

    fig, ax = plt.subplots(figsize=(9,5.2))
    selected=sorted(candidates,key=lambda r: (-r['n_mapping_high_severity'], -r['n_mapping_flags'],r['accession']))[:15]
    ax.barh([r['accession'] for r in selected][::-1], [r['n_mapping_high_severity'] for r in selected][::-1])
    ax.set(xlabel='Число пар геном–белок с критическими замечаниями',
           title='Контроль картирования: записи с наибольшим числом замечаний')
    fig.tight_layout(); p=figures/'04_mapping_qc_flags.png'
    fig.savefig(p,dpi=200);plt.close(fig);created.append(p)
    return created


def run(args):
    root = args.root.resolve()
    out = (args.out or root/'12_final_report').resolve()
    paths={k:root/v for k,v in INPUTS.items()}
    missing=[f'{k}: {p}' for k,p in paths.items() if not p.is_file()]
    if missing:
        raise FileNotFoundError('Обязательные данные отсутствуют:\n'+'\n'.join(missing))
    optional_paths={k:root/v for k,v in OPTIONAL.items() if (root/v).is_file()}
    annotation=load_tsv(paths['annotations'], ['accession','protein_name','reviewed','evidence_group'])
    human=load_tsv(paths['human'], ['accession','has_reported_hit','best_pident_pct','best_query_coverage_pct'])
    mapping=load_tsv(paths['mapping'],['candidate','n_genomes','n_pass','n_high_severity','n_any_flags'])
    mapping_clade=load_tsv(paths['mapping_by_clade'],['candidate','clade','n_genomes','n_pass'])
    conservation=load_tsv(paths['conservation'],['candidate','n_genomes_total','n_genomes_extracted','extraction_pct','n_unique_alleles','mean_entropy_bits','invariant_columns_pct','n_usable_columns'])
    byclade=load_tsv(paths['conservation_by_clade'], ['candidate','clade','n_genomes_in_clade','n_genomes_extracted','extraction_pct'])
    clade_assignment=load_tsv(paths['clade_assignments'], ['assembly','assigned_clade'])
    manifest=load_tsv(paths['genome_manifest'],['assembly'])
    inventory=load_tsv(paths['sequence_inventory'],['assembly','candidate'])
    mafft=load_json(paths['mafft_run'])
    conservation_meta=load_json(paths['conservation_run'])
    if mafft.get('complete') is not True:
        raise ValueError('MAFFT: complete != true; шаг 19 не завершён.')

    ann=index_unique(annotation,'accession','annotations')
    hum=index_unique(human,'accession','human')
    mp=index_unique(mapping,'candidate','mapping')
    cons=index_unique(conservation,'candidate','conservation')
    ids=set(ann)
    for name, obj in [('human',hum),('mapping',mp),('conservation',cons)]:
        if set(obj)!=ids:
            raise ValueError(f'Несовпадение кандидатов: {name}: missing {sorted(ids-set(obj))[:6]}, extra {sorted(set(obj)-ids)[:6]}')
    manifest_idx=index_unique(manifest,'assembly','manifest')
    clad_idx=index_unique(clade_assignment,'assembly','clades')
    if set(manifest_idx)!=set(clad_idx):
        raise ValueError('clade_assignments и genome_manifest содержат разный набор сборок')
    genome_n=len(manifest_idx)
    clade_n=Counter(r['assigned_clade'] if r['assigned_clade'] in CLADES else 'Unassigned' for r in clade_assignment)
    if sum(clade_n.values())!=genome_n:
        raise ValueError('Ошибка количества геномов')

    clade_group={}
    for r in byclade:
        key=(r['candidate'],r['clade'])
        if key in clade_group:
            raise ValueError(f'Повтор строки candidate/clade: {key}')
        clade_group[key]=r
    mapping_clade_idx={}
    for r in mapping_clade:
        key=(r['candidate'],r['clade'])
        if key in mapping_clade_idx:
            raise ValueError(f'Повтор в mapping_by_clade: {key}')
        mapping_clade_idx[key]=r
    expected_keys={(cid,cl) for cid in ids for cl in CLADES+['Unassigned']}
    if set(clade_group)!=expected_keys or set(mapping_clade_idx)!=expected_keys:
        raise ValueError(f'Неполные таблицы по кладам: conservation missing {len(expected_keys-set(clade_group))}, mapping missing {len(expected_keys-set(mapping_clade_idx))}')

    inv_counts=Counter((r['candidate'],r['assembly']) for r in inventory)
    if any(v!=1 for v in inv_counts.values()):
        raise ValueError('Повторные candidate×assembly в sequence_inventory')
    if any(c not in ids or a not in manifest_idx for c,a in inv_counts):
        raise ValueError('sequence_inventory содержит неизвестные белки или геномы')

    out_rows=[]; issues=[]; all_byclade=[]
    sample_clade_n={cl: clade_n[cl] for cl in CLADES+['Unassigned']}
    for cid in sorted(ids):
        a,h,m,c=ann[cid],hum[cid],mp[cid],cons[cid]
        total=integer(c['n_genomes_total']); extracted=integer(c['n_genomes_extracted'])
        mapped=integer(m['n_pass']); flags=integer(m['n_any_flags']); high=integer(m['n_high_severity'])
        if not all(x is not None for x in [total,extracted,mapped,flags,high]):
            raise ValueError(f'{cid}: пропущенные счётчики')
        if total!=genome_n or integer(m['n_genomes'])!=genome_n:
            raise ValueError(f'{cid}: число геномов не равно манифесту ({genome_n})')
        if not 0 <= extracted <= mapped <= total:
            raise ValueError(f'{cid}: n_extracted={extracted} > n_pass={mapped} или другие противоречия')
        if sum(1 for c2,a2 in inv_counts if c2==cid)!=extracted:
            raise ValueError(f'{cid}: n_extracted != фактической sequence_inventory')
        if number(c['extraction_pct']) is not None and abs(float(c['extraction_pct'])-100*extracted/genome_n) > .2:
            raise ValueError(f'{cid}: extraction_pct не соответствует счётчику')
        clade_values=[]
        for cl in CLADES+['Unassigned']:
            z=clade_group[(cid,cl)]; zmap=mapping_clade_idx[(cid,cl)]
            n=integer(z['n_genomes_in_clade']); k=integer(z['n_genomes_extracted'])
            mn=integer(zmap['n_genomes']); mk=integer(zmap['n_pass'])
            if n != clade_n[cl] or mn!=n or k is None or mk is None or not 0<=k<=mk<=n:
                raise ValueError(f'{cid}/{cl}: противоречие между результатами картирования, extraction и размером клады')
            all_byclade.append({
                'accession':cid,'clade':cl,'n_genomes':n,'n_mapped_pass':mk,
                'n_QC_extracted':k,'extraction_pct':maybe(percent(z['extraction_pct'])) if n else '',
                'wilson95_low_pct':z.get('wilson95_low_pct',''),
                'wilson95_high_pct':z.get('wilson95_high_pct',''),
                'n_unique_alleles':z.get('n_unique_alleles',''),
                'mean_entropy_bits':z.get('mean_entropy_bits',''),
                'caution':('no_sample' if n==0 else 'small_n_descriptive_only' if n<args.min_clade_n else z.get('caution','')),
            })
            if cl in CLADES and 0<n<args.min_clade_n:
                clade_values.append(cl)
        if sum(int(clade_group[(cid,cl)]['n_genomes_extracted']) for cl in CLADES+['Unassigned']) != extracted:
            raise ValueError(f'{cid}: сумма извлечённых последовательностей по кладам не совпадает')
        human_hit=yes(h['has_reported_hit'])
        hitid=h.get('best_human_hit','') if human_hit else ''
        # A location ECO:0000269 is a useful evidence flag, not a proof of surface exposure.
        loc_exp=yes(a.get('location_evidence_eco0000269',''))
        issues_c=[]
        if a.get('reviewed','').lower()!='reviewed': issues_c.append('unreviewed_uniprot_entry')
        if not loc_exp: issues_c.append('no_explicit_experimental_localization_evidence_in_catalog')
        if not a.get('cc_subcellular_location','').strip(): issues_c.append('no_subcellular_location_text')
        if high: issues_c.append('high_severity_mapping_flags')
        if flags: issues_c.append('mapping_review_flags')
        if mapped<genome_n: issues_c.append('mapping_not_passed_in_all_genomes')
        if extracted<mapped: issues_c.append('mapped_but_QC_excluded')
        if clade_values: issues_c.append('undersampled_clades_'+','.join(clade_values))
        if clade_n['Unassigned']>0: issues_c.append('unassigned_genomes')
        if human_hit: issues_c.append('reported_human_homology_review')
        if extracted<3: issues_c.append('insufficient_sequences_for_variability')
        for code in issues_c:
            issues.append({'accession':cid,'issue_code':code,'notes': ''})
        out_rows.append({
            'accession':cid,'protein_name':a['protein_name'],
            'uniprot_reviewed':a.get('reviewed',''),'location_annotation_group':a.get('evidence_group',''),
            'location_eco0000269_flag':loc_exp,
            'location_annotation':a.get('cc_subcellular_location',''),
            'go_cellular_component':a.get('go_c',''),
            'ft_signal':a.get('ft_signal',''), 'ft_lipid':a.get('ft_lipid',''),
            'homology_cluster_id':a.get('cluster_id',''),
            'homology_cluster_members':a.get('cluster_members',''),
            'homology_cluster_taxids':a.get('cluster_organism_taxids',''),
            'n_genomes_total':genome_n,'n_miniprot_pass':mapped,'miniprot_pass_pct':round(100*mapped/genome_n,2),
            'n_QC_extracted':extracted,'extracted_pct':round(100*extracted/genome_n,2),
            'wilson95_low_pct':c.get('wilson95_low_pct',''),
            'wilson95_high_pct':c.get('wilson95_high_pct',''),
            'n_unique_alleles':integer(c['n_unique_alleles']),
            'n_alignment_columns':integer(c.get('n_alignment_columns','')),
            'n_usable_alignment_columns':integer(c['n_usable_columns']),
            'mean_entropy_bits':number(c['mean_entropy_bits']),
            'median_entropy_bits':number(c.get('median_entropy_bits','')),
            'invariant_columns_pct':percent(c['invariant_columns_pct']),
            'mean_position_occupancy_pct':percent(c.get('mean_canonical_occupancy_pct','')),
            'human_has_reported_hit':human_hit,
            'human_best_hit':hitid,
            'human_best_identity_pct':percent(h['best_pident_pct']) if human_hit else None,
            'human_best_query_coverage_pct':percent(h['best_query_coverage_pct']) if human_hit else None,
            'human_best_evalue':h.get('best_evalue','') if human_hit else '',
            'n_mapping_high_severity':high,'n_mapping_flags':flags,
            'n_caveats':len(issues_c),'caveats':';'.join(issues_c),
            'analysis_interpretation':'descriptive_comparative_evidence_only',
        })

    overview={
        'generated_utc':datetime.now(timezone.utc).isoformat(),
        'software_version':VERSION,
        'n_candidates':len(ids),'n_genomes':genome_n,
        'n_candidate_genome_comparisons':len(ids)*genome_n,
        'clade_counts':sample_clade_n,
        'n_uniprot_reviewed':sum(r['uniprot_reviewed'].lower()=='reviewed' for r in out_rows),
        'n_with_explicit_experimental_location_flag':sum(r['location_eco0000269_flag'] for r in out_rows),
        'n_human_reported_hits':sum(r['human_has_reported_hit'] for r in out_rows),
        'n_no_human_reported_hits':sum(not r['human_has_reported_hit'] for r in out_rows),
        'n_with_mapping_high_severity_flags':sum(r['n_mapping_high_severity']>0 for r in out_rows),
        'n_without_QC_extracted_sequences':sum(r['n_QC_extracted']==0 for r in out_rows),
        'n_with_any_QC_extracted_sequence':sum(r['n_QC_extracted']>0 for r in out_rows),
        'median_QC_extraction_pct':statistics.median(r['extracted_pct'] for r in out_rows),
        'caveats':[
            'The UniProt annotation-based panel is incomplete and enriched by query design; counts are not antigen prevalence.',
            'ANI reference matching yields provisional clade assignments; phylogenetic validation remains desirable.',
            'Protein recovery and sequence conservation among recovered homologues do not establish expression, surface accessibility or protective immunity.',
            'DIAMOND absence of reported matches is not a safety assessment or exclusion of epitope-level cross-reactivity.',
            'No automated efficacy or safety ranking is computed.'
        ]
    }
    source_rows=[]
    for label, p in sorted({**paths,**optional_paths}.items()):
        source_rows.append({'source_key':label,'relative_path':str(p.relative_to(root)),
                            'bytes':p.stat().st_size,'sha256':sha256(p)})

    if args.check_only:
        print('[OK] Входные таблицы совместимы; запись результатов пропущена (--check-only)')
        print(json.dumps({k:overview[k] for k in ('n_candidates','n_genomes','clade_counts','n_human_reported_hits')},ensure_ascii=False,indent=2))
        return

    out.mkdir(parents=True,exist_ok=True)
    write_tsv(out/'integrated_candidates.tsv',out_rows)
    write_tsv(out/'candidate_by_clade.tsv',all_byclade)
    write_tsv(out/'evidence_issues.tsv',issues if issues else [{'accession':'','issue_code':'no_flags','notes':''}])
    write_tsv(out/'input_manifest_sha256.tsv',source_rows)

    figures=[]
    if not args.no_figures:
        try:
            figures=chart_figures(out,out_rows,byclade,clade_n)
        except ImportError as err:
            raise RuntimeError('Для графиков установи numpy matplotlib: python -m pip install numpy matplotlib') from err
    overview['figure_files']=[str(x.relative_to(out)) for x in figures]
    (out/'final_summary.json').write_text(json.dumps(overview,indent=2,ensure_ascii=False)+'\n',encoding='utf-8')

    represented=', '.join(f'{cl}={clade_n[cl]}' for cl in CLADES)
    human_examples=[r for r in out_rows if r['human_has_reported_hit']]
    examples='; '.join(f"{r['accession']} ({r['protein_name']})" for r in human_examples[:8]) or 'нет зарегистрированных совпадений'
    n_high=overview['n_with_mapping_high_severity_flags']
    n_unknown=clade_n['Unassigned']
    txt=f'''# Сводный биоинформатический отчёт: Candida auris (ТЕН «Ластоваз»)

Создан: {overview['generated_utc']} UTC. Данные объединены автоматически из шагов 01–20.

## Масштаб и основные наблюдения

- **Кандидатных белковых записей:** {len(ids)}.
- **Геномных сборок:** {genome_n}; в матрице **{genome_n*len(ids)} пар** кандидат × сборка.
- **Предварительное распределение по кладам:** {represented}. Неназначенных: {n_unknown}.
- **UniProtKB/Swiss-Prot reviewed среди кандидатов:** {overview['n_uniprot_reviewed']} из {len(ids)}.
- **Отдельный флаг экспериментальных данных о локализации (ECO:0000269):** {overview['n_with_explicit_experimental_location_flag']} записей. Это только один тип свидетельства.
- **Белков со значимыми совпадениями с человеком в представленном поиске DIAMOND:** {overview['n_human_reported_hits']} из {len(ids)}; остальные {overview['n_no_human_reported_hits']} не имеют сообщённых совпадений при прежних условиях поиска.
- **Записи со строгими QC-проблемами картирования хотя бы в одном геноме:** {n_high}.
- **Медиана доли геномов с QC-извлечённой последовательностью:** {overview['median_QC_extraction_pct']:.1f}%.
- **Примеры белков с зарегистрированным сходством с человеком:** {examples}.

## Что означают таблицы

1. `integrated_candidates.tsv` — объединение аннотаций, результатов MMseqs2, miniprot, MAFFT, DIAMOND. **Строки отсортированы по accession, а не по эффективности вакцины.**
2. `candidate_by_clade.tsv` — отдельные показатели по шести предварительным кладам и группе Unassigned. Если `n_genomes=0`, значение пустое, а не 0%.
3. `evidence_issues.tsv` — причины, по которым интерпретацию конкретного кандидата следует ограничить.
4. `input_manifest_sha256.tsv` — контрольные суммы всех использованных файлов для воспроизводимости.
5. `figures/` — иллюстрации расчётов с честными подписями.

## Главные методологические ограничения

- Сборки выбраны по качеству, а не вероятностной выборкой из популяции. По ним нельзя оценить эпидемиологическую распространённость антигенов во всём виде.
- Доля QC-извлечённых последовательностей характеризует успешный поиск/выделение гомологов, **не экспрессию белка и не его наличие на клеточной поверхности**.
- Энтропия относится только к наблюдаемым выравниваемым позициям и зависит от представительства геномов и качества выравнивания.
- Клады определены сравнением ANI с референсами и являются предварительными; требуется филогенетическая проверка.
- Белки с отсутствием значимого DIAMOND-совпадения **не объявляются безопасными**. Структурная, локальная и иммунологическая перекрёстная реактивность не исключены.
- Вычислительный анализ **не доказывает протективный иммунный ответ или эффективность вакцины**. Необходима критическая оценка литературы и последующая экспериментальная валидация.

## Что ещё не сделано автоматически

- Проверка публикаций для каждого белка, включая экспериментальную наружную локализацию и данные об иммуногенности.
- Независимая SNP/геномная филогенетическая верификация клад.
- Экспериментальная проверка защитного эффекта и побочных реакций.
- Выбор состава вакцины или дизайн эпитопов. Их нельзя обоснованно вычислить только по этим данным.

**Финальный научный вывод:** получен воспроизводимый сравнительный каталог свойств потенциальных поверхностных белков *C. auris* с количественными характеристиками, источниками и ограничениями. Он пригоден для обоснования гипотез, но не является доказательством эффективности вакцины.
'''
    (out/'RESULTS_RU.md').write_text(txt,encoding='utf-8')
    provenance = root / "00_provenance"
    exported_provenance = []
    for name in ("inputs_sha256.tsv", "config.json"):
        if (provenance / name).is_file():
            copy = out / ("raw_" + name)
            copy.write_bytes((provenance / name).read_bytes())
            exported_provenance.append(copy)
    reproducibility = """# Воспроизводимость

См. README репозитория для запуска Snakemake и описания входных данных.

- `input_manifest_sha256.tsv` — SHA-256 таблиц, использованных для сборки отчёта.
- `raw_inputs_sha256.tsv` и `raw_config.json` — входные снимки и параметры запуска, если отчёт собран через Snakemake.

Архив содержит таблицы, графики и отчёт без исходных геномов.
"""
    (out/'README_RU.md').write_text(reproducibility,encoding='utf-8')
    archive=out/'final_results_bundle.zip'
    with zipfile.ZipFile(archive,'w',compression=zipfile.ZIP_DEFLATED,compresslevel=6) as z:
        tables_and_text = [
            out / name for name in (
                "integrated_candidates.tsv", "candidate_by_clade.tsv", "evidence_issues.tsv",
                "input_manifest_sha256.tsv", "final_summary.json", "RESULTS_RU.md", "README_RU.md",
            )
        ]
        for path in sorted([*tables_and_text, *exported_provenance, *figures]):
            z.write(path, path.relative_to(out))
    print('[OK] Кандидатов:',len(ids),'Геномов:',genome_n)
    print('[OK] Клады:', ', '.join(f'{x}:{clade_n[x]}' for x in CLADES))
    print('[OK] Файлов графиков:',len(figures))
    print('[DONE] Отчёт:',out/'RESULTS_RU.md')
    print('[DONE] Таблица:',out/'integrated_candidates.tsv')
    print('[DONE] Архив:',archive)


def main():
    parser=argparse.ArgumentParser(description='TEN step21: integrate evidence, tables, QC and descriptive figures')
    parser.add_argument('--root',type=Path,default=Path(__file__).resolve().parent.parent)
    parser.add_argument('--out',type=Path)
    parser.add_argument('--min-clade-n',type=int,default=5)
    parser.add_argument('--check-only',action='store_true')
    parser.add_argument('--no-figures',action='store_true')
    args=parser.parse_args()
    if args.min_clade_n < 2:
        parser.error('--min-clade-n должен быть >=2')
    try:
        run(args)
    except (FileNotFoundError,ValueError,KeyError,RuntimeError,ZeroDivisionError) as e:
        print('[ERROR]',str(e),file=sys.stderr)
        return 1
    return 0


if __name__=='__main__':
    sys.exit(main())
