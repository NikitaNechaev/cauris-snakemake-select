#!/usr/bin/env python3
"""TEN / C. auris / step 20: descriptive amino-acid conservation statistics.

Read the exact FASTA/TSV format produced by steps 18 and 19. All outputs go
into 11_conservation_statistics. No vaccine efficacy, exposure or safety claims.

Inputs (relative to project root):
  08_population/qc/genome_manifest.tsv
  08_population/results/clade_typing/clade_assignments.tsv
  09_protein_sequences/sequence_inventory.tsv
  10_mafft/mafft_summary.json
  10_mafft/alignment_summary.tsv
  10_mafft/alignments_all/<candidate>.aln.faa
  09_protein_sequences/by_candidate/<candidate>.faa

Outputs:
  11_conservation_statistics/summary_by_candidate.tsv
  11_conservation_statistics/summary_by_clade.tsv
  11_conservation_statistics/clade_coverage_matrix.tsv
  11_conservation_statistics/position_profiles/*.tsv
  11_conservation_statistics/figures/*.png (if matplotlib installed)
  11_conservation_statistics/conservation_summary.json

Core definitions:
  - n_extracted: genomes with QC-approved translated miniprot alignment;
    NOT population prevalence or evidence of expression.
  - positional occupancy: canonical amino-acid residues / n_extracted.
    Gaps and ambiguous X/U/O do not enter entropy.
  - Shannon H = -sum(p_aa * log2(p_aa)) for A,C,D,...,Y.
  - Summary entropy/invariance: ONLY alignment columns with occupancy
    >= --min-occupancy. Invariance is within OBSERVED homologues.
  - By-clade counts are reported separately, with no balancing implied.

Python 3.9+, standard library; matplotlib optional for charts.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import os
import statistics
import sys
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path

VERSION = "1.0.0"
AA = set("ACDEFGHIKLMNPQRSTVWY")
OTHER = set("XUO")
GAP = set("-.")
CLADES = ["I", "II", "III", "IV", "V", "VI"]


def read_tsv(path: Path, required=()):
    if not path.is_file():
        raise FileNotFoundError(f"Нет исходного файла: {path}")
    with path.open(encoding="utf-8", newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        missing = set(required) - set(reader.fieldnames or ())
        if missing:
            raise ValueError(f"В {path} отсутствуют столбцы: {sorted(missing)}; имеются: {reader.fieldnames}")
        return list(reader)


def write_tsv(path: Path, fieldnames, rows):
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(path.name + ".partial")
    with tmp.open("w", encoding="utf-8", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)
    os.replace(tmp, path)


def write_json(path: Path, payload):
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(path.name + ".partial")
    tmp.write_text(json.dumps(payload, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    os.replace(tmp, path)


def file_sha(path: Path):
    digest = hashlib.sha256()
    with path.open("rb") as fh:
        for buf in iter(lambda: fh.read(1024 * 1024), b""):
            digest.update(buf)
    return digest.hexdigest()


def fasta(path: Path, *, alignment=False):
    if not path.is_file():
        raise FileNotFoundError(f"Отсутствует FASTA: {path}")
    records, ident, chunks = [], None, []
    seen = set()

    def flush():
        if ident is not None:
            sequence = "".join(chunks).upper()
            if not sequence:
                raise ValueError(f"Пустая FASTA-запись {ident} в {path}")
            allowed = AA | OTHER | (GAP if alignment else set())
            unexpected = set(sequence) - allowed
            if unexpected:
                raise ValueError(f"Недопустимые символы {unexpected} у {ident} в {path}")
            records.append((ident, sequence))

    with path.open(encoding="utf-8") as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                flush()
                ident = line[1:].split()[0]
                if not ident or ident in seen:
                    raise ValueError(f"Повтор/пустой ID FASTA: {ident!r} ({path})")
                seen.add(ident)
                chunks = []
            else:
                if ident is None:
                    raise ValueError(f"Последовательность без FASTA-заголовка ({path})")
                chunks.append(line)
    flush()
    return records


def wilson(k, n, z=1.959963984540054):
    if n == 0:
        return (None, None)
    p = k / n
    denom = 1 + z * z / n
    center = (p + z * z / (2 * n)) / denom
    width = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / denom
    return (100 * max(0, center - width), 100 * min(1, center + width))


def round_or_blank(v, digits=4):
    return "" if v is None else round(v, digits)


def position_stats(records, min_occ):
    """records are (FASTA_ID, aligned protein); same-length, nonempty."""
    if not records:
        return [], {}
    lengths = {len(seq) for _, seq in records}
    if len(lengths) != 1:
        raise ValueError(f"MAFFT alignment contains unequal lengths: {lengths}")
    cols = next(iter(lengths))
    n = len(records)
    positions = []
    for pos in range(cols):
        column = [seq[pos] for _, seq in records]
        canonical = [c for c in column if c in AA]
        counts = Counter(canonical)
        n_aa = len(canonical)
        gap_count = sum(c in GAP for c in column)
        other_count = n - n_aa - gap_count
        occupied = n_aa / n
        usable = occupied >= min_occ and n_aa > 0
        if counts:
            major, major_count = sorted(counts.items(), key=lambda x: (-x[1], x[0]))[0]
            probs = [cnt / n_aa for cnt in counts.values()]
            entropy = -sum(p * math.log2(p) for p in probs)
            major_frequency = 100 * major_count / n_aa
        else:
            major = ""
            entropy = None
            major_frequency = None
        positions.append({
            "position": pos + 1, "n_sequences": n,
            "n_canonical_aa": n_aa, "n_gap": gap_count, "n_ambiguous": other_count,
            "aa_occupancy_pct": 100 * occupied,
            "gap_pct": 100 * gap_count / n,
            "major_aa": major, "major_aa_frequency_pct": round_or_blank(major_frequency),
            "n_distinct_aa": len(counts),
            "entropy_bits": round_or_blank(entropy, 6),
            "passes_occupancy": int(usable),
            "invariant_among_observed": int(len(counts) == 1) if usable else "",
        })
    passing = [r for r in positions if r["passes_occupancy"]]
    summary = {
        "n_alignment_columns": cols,
        "n_usable_columns": len(passing),
        "mean_entropy_bits": statistics.mean([r["entropy_bits"] for r in passing]) if passing else None,
        "median_entropy_bits": statistics.median([r["entropy_bits"] for r in passing]) if passing else None,
        "invariant_columns_pct": 100 * sum(r["invariant_among_observed"] == 1 for r in passing) / len(passing) if passing else None,
        "mean_major_aa_frequency_pct": statistics.mean([r["major_aa_frequency_pct"] for r in passing]) if passing else None,
        "mean_canonical_occupancy_pct": statistics.mean([r["aa_occupancy_pct"] for r in positions]) if positions else None,
    }
    return positions, summary


def chart_figures(out, overall, byclade, clade_n, plot_top):
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        import numpy as np
    except ImportError as exc:
        print(f"[WARN] Графики пропущены: {exc}. Для графиков: python -m pip install matplotlib numpy", file=sys.stderr)
        return []
    figures = out / "figures"
    figures.mkdir(parents=True, exist_ok=True)
    files = []
    clade_map = {(r["candidate"], r["clade"]): r for r in byclade}
    ordered = sorted(overall, key=lambda r: (-int(r["n_genomes_extracted"]), r["candidate"]))
    matrix = []
    for r in ordered:
        vals = []
        for clade in CLADES:
            data = clade_map[(r["candidate"], clade)]
            vals.append(float(data["extraction_pct"]) if data["extraction_pct"] != "" else float("nan"))
        matrix.append(vals)
    values = np.ma.masked_invalid(np.array(matrix, dtype=float))
    fig, ax = plt.subplots(figsize=(9, max(9, 0.27 * len(ordered) + 2)))
    heat = ax.imshow(values, aspect="auto", vmin=0, vmax=100)
    ax.set_xticks(range(len(CLADES)), [f"{c} (n={clade_n[c]})" for c in CLADES])
    ax.set_yticks(range(len(ordered)), [r["candidate"] for r in ordered], fontsize=6)
    ax.set_xlabel("Предварительная клада (число сборок)")
    ax.set_ylabel("Белок UniProt accession")
    ax.set_title("Доля геномов с извлечённым QC-прошедшим гомологом, %")
    fig.colorbar(heat, ax=ax, label="% геномов в кладе")
    fig.tight_layout()
    path = figures / "recovery_by_clade_heatmap.png"
    fig.savefig(path, dpi=200)
    plt.close(fig)
    files.append(str(path))

    xs, ys, labels = [], [], []
    for r in overall:
        if r["mean_entropy_bits"] == "":
            continue
        xs.append(float(r["extraction_pct"]))
        ys.append(float(r["mean_entropy_bits"]))
        labels.append(r["candidate"])
    if xs:
        fig, ax = plt.subplots(figsize=(8, 5))
        ax.scatter(xs, ys)
        ax.set_xlabel("Сборки с QC-прошедшим гомологом, %")
        ax.set_ylabel("Средняя позиционная энтропия, бит")
        ax.set_title("Описательная вариабельность среди извлечённых гомологов")
        for i in range(min(plot_top, len(xs))):
            ax.annotate(labels[i], (xs[i], ys[i]), fontsize=7)
        fig.tight_layout()
        path = figures / "recovery_vs_entropy.png"
        fig.savefig(path, dpi=200)
        plt.close(fig)
        files.append(str(path))
    return files


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--root", type=Path, default=Path(__file__).resolve().parent.parent,
                    help="Путь к корню TEN (по умолчанию родитель каталога scripts)")
    ap.add_argument("--out", type=Path, help="По умолчанию <root>/11_conservation_statistics")
    ap.add_argument("--min-occupancy", type=float, default=0.8,
                    help="Минимальная доля канонических AA в колонке для расчёта сводных метрик (0–1)")
    ap.add_argument("--min-clade-n", type=int, default=5,
                    help="Малые клады (<n) отмечать как описательные, без уверенного вывода")
    ap.add_argument("--plot-top", type=int, default=5,
                    help="Сколько точек подписать на scatter-диаграмме")
    args = ap.parse_args()
    if not 0 < args.min_occupancy <= 1 or args.min_clade_n < 2 or args.plot_top < 0:
        ap.error("Проверь --min-occupancy (0,1], --min-clade-n (>=2), --plot-top (>=0)")
    root = args.root.expanduser().resolve()
    out = (args.out or root / "11_conservation_statistics").expanduser().resolve()
    pop = root / "08_population"
    sequences = root / "09_protein_sequences"
    mafft = root / "10_mafft"
    clade_file = pop / "results/clade_typing/clade_assignments.tsv"
    manifest_file = pop / "qc/genome_manifest.tsv"
    inventory_file = sequences / "sequence_inventory.tsv"
    alignment_file = mafft / "alignment_summary.tsv"
    run_file = mafft / "mafft_summary.json"
    if not run_file.is_file():
        raise FileNotFoundError(f"Не найден итог MAFFT: {run_file}")
    run_info = json.loads(run_file.read_text(encoding="utf-8"))
    if run_info.get("complete") is not True:
        raise RuntimeError("Шаг 19 MAFFT не завершён (mafft_summary.json: complete != true)")
    manifest = read_tsv(manifest_file, ("assembly",))
    clade_assignments = read_tsv(clade_file, ("assembly", "assigned_clade"))
    inventory = read_tsv(inventory_file, ("assembly", "candidate", "clade"))
    alignment_summary = read_tsv(alignment_file, ("candidate", "status", "n_genomes", "n_unique_alleles"))
    assemblies = [r["assembly"] for r in manifest]
    if len(set(assemblies)) != len(assemblies) or not assemblies:
        raise ValueError("Манифест содержит пустой список или повторные assembly")
    clades = {}
    for r in clade_assignments:
        acc = r["assembly"]
        if acc in clades:
            raise ValueError("Повтор в clade_assignments: " + acc)
        c = r["assigned_clade"] or "Unassigned"
        if c not in CLADES + ["Unassigned"]:
            raise ValueError(f"Неизвестная клада {c!r} для {acc}")
        clades[acc] = c
    if set(clades) != set(assemblies):
        raise ValueError("Набор assembly в clade_assignments и genome_manifest не совпадает")
    inventory_keys = {}
    for r in inventory:
        key = (r["candidate"], r["assembly"])
        if key in inventory_keys:
            raise ValueError("Повтор в sequence_inventory.tsv: " + str(key))
        if r["assembly"] not in clades or r["clade"] != clades[r["assembly"]]:
            raise ValueError(
                f"Устаревшая кладовая метка в sequence_inventory для {key}: "
                f"{r['clade']!r} != {clades.get(r['assembly'])!r}. "
                "Синхронизируй шаги 17–19 перед анализом."
            )
        inventory_keys[key] = r
    if len(alignment_summary) != len({r["candidate"] for r in alignment_summary}):
        raise ValueError("Дубликаты candidate в alignment_summary.tsv")
    if not alignment_summary:
        raise ValueError("Пустая alignment_summary.tsv")
    if run_info.get("n_candidates_reported") != len(alignment_summary):
        raise ValueError("Количество кандидатов в MAFFT summary не совпадает с таблицей")
    out.mkdir(parents=True, exist_ok=True)
    profile_dir = out / "position_profiles"
    profile_dir.mkdir(exist_ok=True)
    clade_n = {c: sum(clades[a] == c for a in assemblies) for c in CLADES + ["Unassigned"]}
    by_cand_rows, by_clade_rows, profile_counts = [], [], {}
    clade_matrix_rows = []
    for i, row in enumerate(alignment_summary, 1):
        cid = row["candidate"]
        expected_n = int(row["n_genomes"])
        data = {}
        for label, subset in (("all", assemblies),):
            aln = mafft / "alignments_all" / f"{cid}.aln.faa"
            raw = sequences / "by_candidate" / f"{cid}.faa"
            if expected_n == 0:
                if aln.is_file() and aln.stat().st_size > 0:
                    raise ValueError(f"{cid}: пустой кандидат имеет непустое выравнивание")
                if raw.is_file() and raw.stat().st_size > 0:
                    raise ValueError(f"{cid}: пустой кандидат имеет непустой исходный FASTA")
                data = {}
                ncols = 0
                stats = {}
                records = []
            else:
                records = fasta(aln, alignment=True)
                raw_records = fasta(raw)
                data = {}
                raw_by_id = dict(raw_records)
                if len(records) != expected_n or len(raw_records) != expected_n:
                    raise ValueError(f"{cid}: число FASTA-записей не равно n_genomes ({expected_n})")
                if set(n for n, _ in records) != set(raw_by_id):
                    raise ValueError(f"{cid}: изменились FASTA ID между 09 и 10 этапами")
                for name, seq in records:
                    prefix = cid + "|"
                    if not name.startswith(prefix):
                        raise ValueError(f"{cid}: неправильный заголовок FASTA {name}")
                    acc = name[len(prefix):]
                    if acc not in clades:
                        raise ValueError(f"{cid}: неизвестный геном {acc}")
                    if seq.replace("-", "").replace(".", "") != raw_by_id[name]:
                        raise ValueError(f"{cid}: MAFFT изменил исходную последовательность {name}")
                    if (cid, acc) not in inventory_keys:
                        raise ValueError(f"{cid}: отсутствует запись sequence_inventory для {acc}")
                    data[acc] = seq
                if {a for c, a in inventory_keys if c == cid} != set(data):
                    raise ValueError(f"{cid}: sequence_inventory содержит не те геномы, что FASTA")
                prof, stats = position_stats(records, args.min_occupancy)
                ncols = stats["n_alignment_columns"]
                write_tsv(profile_dir / f"{cid}.tsv", list(prof[0]), prof)
                profile_counts[cid] = len(prof)
        if expected_n != len(data):
            raise ValueError(f"{cid}: несоответствие количества записей")
        n_unique = len({seq for seq in data.values()})
        if n_unique != int(row["n_unique_alleles"]):
            raise ValueError(f"{cid}: аллели из FA по количеству не совпадают с MAFFT summary")
        lo, hi = wilson(expected_n, len(assemblies))
        candidate_row = {
            "candidate": cid, "n_genomes_total": len(assemblies),
            "n_genomes_extracted": expected_n,
            "extraction_pct": round_or_blank(100 * expected_n / len(assemblies), 2),
            "wilson95_low_pct": round_or_blank(lo, 2),
            "wilson95_high_pct": round_or_blank(hi, 2),
            "n_unique_alleles": n_unique,
            "n_alignment_columns": ncols,
            "n_usable_columns": stats.get("n_usable_columns", 0),
            "mean_entropy_bits": round_or_blank(stats.get("mean_entropy_bits")),
            "median_entropy_bits": round_or_blank(stats.get("median_entropy_bits")),
            "invariant_columns_pct": round_or_blank(stats.get("invariant_columns_pct"), 2),
            "mean_major_aa_frequency_pct": round_or_blank(stats.get("mean_major_aa_frequency_pct"), 2),
            "mean_canonical_occupancy_pct": round_or_blank(stats.get("mean_canonical_occupancy_pct"), 2),
            "interpretation_flag": ("no_QC_extracted_sequences" if expected_n == 0 else
                                    "few_genomes_descriptive_only" if expected_n < 3 else
                                    "single_allele_observed" if n_unique == 1 else "descriptive"),
        }
        by_cand_rows.append(candidate_row)
        by_clade = {}
        for clade in CLADES + ["Unassigned"]:
            group = [(cid + "|" + acc, seq) for acc, seq in data.items() if clades[acc] == clade]
            n = clade_n[clade]
            k = len(group)
            lo, hi = wilson(k, n)
            _, cstats = position_stats(group, args.min_occupancy)
            cr = {
                "candidate": cid, "clade": clade,
                "n_genomes_in_clade": n,
                "n_genomes_extracted": k,
                "extraction_pct": round_or_blank(100 * k / n, 2) if n else "",
                "wilson95_low_pct": round_or_blank(lo, 2),
                "wilson95_high_pct": round_or_blank(hi, 2),
                "n_unique_alleles": len({seq for _, seq in group}),
                "n_usable_columns": cstats.get("n_usable_columns", 0),
                "mean_entropy_bits": round_or_blank(cstats.get("mean_entropy_bits")),
                "invariant_columns_pct": round_or_blank(cstats.get("invariant_columns_pct"), 2),
                "caution": ("clade_absent" if n == 0 else
                            "small_clade_n" if n < args.min_clade_n else
                            "few_extracted_sequences" if k < 3 else "descriptive_only"),
            }
            by_clade_rows.append(cr)
            by_clade[clade] = cr
        matrix_row = {"candidate": cid}
        matrix_row.update({clade: by_clade[clade]["extraction_pct"] for clade in CLADES})
        clade_matrix_rows.append(matrix_row)
        print(f"[OK] {i}/{len(alignment_summary)} {cid}: "+
              f"n={expected_n}/{len(assemblies)}, alleles={n_unique}, columns={ncols}", flush=True)
    overall_fields = list(by_cand_rows[0])
    by_clade_fields = list(by_clade_rows[0])
    write_tsv(out / "summary_by_candidate.tsv", overall_fields, by_cand_rows)
    write_tsv(out / "summary_by_clade.tsv", by_clade_fields, by_clade_rows)
    write_tsv(out / "clade_coverage_matrix.tsv", ["candidate"] + CLADES, clade_matrix_rows)
    figures = chart_figures(out, by_cand_rows, by_clade_rows, clade_n, args.plot_top)
    summary = {
        "timestamp_utc": datetime.now(timezone.utc).isoformat(),
        "version": VERSION,
        "n_candidates": len(by_cand_rows),
        "n_genomes_manifest": len(assemblies),
        "n_QC_extracted_protein_sequences": sum(int(r["n_genomes_extracted"]) for r in by_cand_rows),
        "n_candidates_with_sequences": sum(int(r["n_genomes_extracted"]) > 0 for r in by_cand_rows),
        "n_position_profiles": len(profile_counts),
        "n_figures": len(figures),
        "clade_sizes": clade_n,
        "min_occupancy": args.min_occupancy,
        "min_clade_n_for_caution": args.min_clade_n,
        "source_sha256": {"genome_manifest": file_sha(manifest_file),
                          "clade_assignments": file_sha(clade_file),
                          "sequence_inventory": file_sha(inventory_file),
                          "alignment_summary": file_sha(alignment_file),
                          "mafft_summary": file_sha(run_file)},
        "interpretation_limits": [
            "n_genomes_extracted describes successful mapped protein sequence recovery, not biological prevalence or antigen expression.",
            "Positional entropy is conditional on the homologues passing previous QC; failures are excluded rather than scored as missing residues.",
            "Allele frequencies and entropy across the assembly panel are affected by sampling of related isolates and by unequal clade sizes.",
            "Wilson intervals are descriptive intervals for this selected assembly panel and do not resolve non-random sampling.",
            "These calculations cannot establish vaccine protection, antigen accessibility or absence of human cross-reactivity."
        ],
        "figure_files": figures,
    }
    write_json(out / "conservation_summary.json", summary)
    print("[DONE] " + json.dumps({k: summary[k] for k in (
        "n_candidates", "n_genomes_manifest", "n_QC_extracted_protein_sequences", "n_position_profiles", "n_figures")}, ensure_ascii=False))
    print("[OUTPUT] " + str(out))


if __name__ == "__main__":
    try:
        main()
    except (OSError, ValueError, RuntimeError, KeyError, TypeError) as exc:
        print("[ERROR] " + str(exc), file=sys.stderr)
        sys.exit(1)
