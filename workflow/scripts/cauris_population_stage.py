#!/usr/bin/env python3
"""Project TEN: QC, protein-to-genome mapping, and summaries for the C. auris panel.

Uses the staged project layout unless paths are supplied explicitly.
Python >= 3.9; only Python standard library. Requires miniprot for `map`.
The analysis reports sequence alignment; it does not establish orthology, antigen
expression, surface exposure, cross-reactivity, vaccine safety, or efficacy.
"""
from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import os
import re
import shutil
import subprocess
import sys
import tempfile
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path
from urllib.parse import unquote

ACCESSION_RE = re.compile(r"GC[AF]_\d+\.\d+")
FASTA_SUFFIXES = (".fna", ".fna.gz", ".fa", ".fa.gz", ".fasta", ".fasta.gz")


def default_paths(root: Path) -> dict[str, Path]:
    pop = root / "08_population"
    return {
        "root": root,
        "selected": pop / "cauris_selected_100_package" / "selected_accessions.txt",
        "genomes": pop / "ncbi_dataset" / "data",
        "candidates": root / "04_summarize_clusters" / "auris_priority.faa",
        "qc": pop / "qc",
        "manifest": pop / "qc" / "genome_manifest.tsv",
        "alignments": pop / "alignments",
        "results": pop / "results",
        "ani": pop / "skani_matrix.txt",
    }


def get_paths(args):
    root = Path(args.root).expanduser().resolve() if args.root else Path(__file__).resolve().parent.parent
    paths = default_paths(root)
    for key in ("selected", "genomes", "candidates", "qc", "manifest", "alignments", "results"):
        value = getattr(args, key, None)
        if value:
            paths[key] = Path(value).expanduser().resolve()
    return paths


def open_text(path):
    path = Path(path)
    if str(path).endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8")
    return path.open("r", encoding="utf-8")


def write_tsv(path, header, rows):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t")
        writer.writerow(header)
        writer.writerows(rows)


def write_json(path, obj):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(obj, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")


def digest_file(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as fh:
        for block in iter(lambda: fh.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def read_selected(path):
    path = Path(path)
    if not path.is_file():
        raise FileNotFoundError(f"Список accession не найден: {path}")
    with open_text(path) as fh:
        first = fh.readline().strip("\ufeff\n\r")
    with open_text(path) as fh:
        if first.startswith("Assembly Accession\t"):
            accessions = [r["Assembly Accession"].strip() for r in csv.DictReader(fh, delimiter="\t")]
        else:
            accessions = [line.strip().split()[0] for line in fh if line.strip() and not line.startswith("#")]
    if not accessions or any(not ACCESSION_RE.fullmatch(acc) for acc in accessions):
        raise ValueError(f"Некорректные accession в {path}")
    if len(set(accessions)) != len(accessions):
        raise ValueError(f"Повторяющиеся accession в {path}")
    return accessions


def read_manifest(path):
    path = Path(path)
    if not path.is_file():
        raise FileNotFoundError(f"Манифест не найден: {path}; сначала выполните команду qc")
    with path.open(encoding="utf-8", newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        if not reader.fieldnames or not {"assembly", "fasta"}.issubset(reader.fieldnames):
            raise ValueError("Ожидаются столбцы assembly и fasta в genome_manifest.tsv")
        manifest = list(reader)
    if not manifest:
        raise ValueError(f"Пустой манифест: {path}")
    assemblies = [row["assembly"] for row in manifest]
    if len(assemblies) != len(set(assemblies)):
        raise ValueError("Дубликаты accession в манифесте")
    return manifest


def fasta_lengths(path):
    result = {}
    current = None
    with open_text(path) as fh:
        for line in fh:
            if line.startswith(">"):
                current = line[1:].split()[0]
                if not current:
                    raise ValueError(f"Пустой FASTA ID: {path}")
                if current in result:
                    raise ValueError(f"Дублирующийся FASTA ID: {current}")
                result[current] = 0
            elif current is not None:
                result[current] += len(line.strip())
            elif line.strip():
                raise ValueError(f"Неверный FASTA (нет строки >): {path}")
    if not result or any(length <= 0 for length in result.values()):
        raise ValueError(f"Пустой FASTA или пустые последовательности: {path}")
    return result


def genome_stats(path):
    lengths, hashes = [], []
    seq_parts = []
    total = ambiguous = gc = unusual = 0

    def flush():
        nonlocal total, ambiguous, gc, unusual
        if not seq_parts:
            return
        seq = "".join(seq_parts).upper()
        total += len(seq)
        ambiguous += seq.count("N")
        gc += seq.count("G") + seq.count("C")
        unusual += sum(c not in "ACGTN" for c in seq)
        lengths.append(len(seq))
        hashes.append(hashlib.sha256(seq.encode("ascii", "replace")).hexdigest())
        seq_parts.clear()

    saw_header = False
    with open_text(path) as fh:
        for line in fh:
            if line.startswith(">"):
                saw_header = True
                flush()
            else:
                if line.strip() and not saw_header:
                    raise ValueError("Нет FASTA-заголовка")
                seq_parts.append(line.strip())
        flush()
    if total == 0 or not lengths:
        raise ValueError("Пустой или некорректный FASTA")
    lengths.sort(reverse=True)
    running = n50 = 0
    for size in lengths:
        running += size
        if running >= total / 2:
            n50 = size
            break
    combined_hash = hashlib.sha256("\n".join(sorted(hashes)).encode()).hexdigest()
    return (len(lengths), total, n50, lengths[0], round(ambiguous * 100 / total, 4),
            round(gc * 100 / total, 4), unusual, combined_hash)


def qc(args, paths):
    selected = read_selected(paths["selected"])
    folder = paths["genomes"]
    if not folder.is_dir():
        raise FileNotFoundError(f"Каталог геномов не найден: {folder}")
    out = paths["qc"]
    out.mkdir(parents=True, exist_ok=True)
    expected = set(selected)
    buckets = defaultdict(list)
    for folder_acc in sorted(folder.iterdir()):
        if not folder_acc.is_dir() or folder_acc.name not in expected:
            continue
        for fna in folder_acc.iterdir():
            if fna.is_file() and fna.name.lower().endswith(FASTA_SUFFIXES):
                buckets[folder_acc.name].append(fna.resolve())
    qc_records, manifest, issues, seen_hashes = [], [], [], {}
    for acc in selected:
        files = buckets.get(acc, [])
        if len(files) == 0:
            issues.append((acc, "ERROR", "Нет FASTA в каталоге accession"))
            continue
        if len(files) > 1:
            issues.append((acc, "ERROR", "Несколько FASTA для accession: " + ";".join(str(f) for f in files)))
            continue
        fna = files[0]
        try:
            nseq, total, n50, longest, n_pct, gc_pct, other, seq_hash = genome_stats(fna)
        except Exception as exc:
            issues.append((acc, "ERROR", f"FASTA QC: {exc}"))
            continue
        if not 11_800_000 <= total <= 13_500_000:
            issues.append((acc, "WARNING", f"Необычный размер генома: {total}"))
        if nseq > 50:
            issues.append((acc, "WARNING", f"Более 50 контингов: {nseq}"))
        if n_pct > 1:
            issues.append((acc, "WARNING", f"N > 1%: {n_pct}"))
        if other:
            issues.append((acc, "WARNING", f"Необычные символы, кроме ACGTN: {other}"))
        if seq_hash in seen_hashes:
            issues.append((acc, "WARNING", f"Точное совпадение последовательностей с {seen_hashes[seq_hash]}"))
        else:
            seen_hashes[seq_hash] = acc
        manifest.append((acc, str(fna)))
        qc_records.append((acc, str(fna), nseq, total, n50, longest, n_pct, gc_pct, other, seq_hash))
    write_tsv(out / "genome_manifest.tsv", ["assembly", "fasta"], manifest)
    write_tsv(out / "genome_qc.tsv", ["assembly", "fasta", "n_contigs", "length_bp", "N50_bp", "longest_contig_bp", "N_pct", "GC_pct", "non_ACGTN", "sequence_hash"], qc_records)
    write_tsv(out / "genome_issues.tsv", ["assembly", "severity", "message"], issues)
    (out / "genome_paths.txt").write_text("".join(path + "\n" for _, path in manifest), encoding="utf-8")
    mapped = {a for a, _ in manifest}
    (out / "genome_missing.txt").write_text("".join(a + "\n" for a in selected if a not in mapped), encoding="utf-8")
    levels = Counter(severity for _, severity, _ in issues)
    summary = {
        "expected_assemblies": len(selected),
        "valid_genome_paths": len(manifest),
        "n_errors": levels["ERROR"],
        "n_warnings": levels["WARNING"],
        "exact_duplicate_assemblies": sum("Точное совпадение" in txt for _, _, txt in issues),
        "selected_accessions": str(paths["selected"]),
        "genome_source": str(paths["genomes"]),
        "ignored_sibling_folder": str(paths["root"] / "08_population" / "data"),
        "qc_dir": str(out),
    }
    write_json(out / "qc_summary.json", summary)
    print(json.dumps(summary, indent=2, ensure_ascii=False))
    if levels["ERROR"]:
        raise RuntimeError("QC выявил ошибки. См. 08_population/qc/genome_issues.tsv")


def file_state(path):
    stat = Path(path).stat()
    return {"path": str(Path(path).resolve()), "size": stat.st_size, "mtime_ns": stat.st_mtime_ns}


def make_stamp(assembly, genome, candidates_hash, tool, code, threads):
    return {
        "assembly": assembly,
        "genome": file_state(genome),
        "candidates_sha256": candidates_hash,
        "tool": tool,
        "genetic_code": code,
        "threads": threads,
        "options": ["-T", str(code), "-I", "--gff"],
    }


def valid_stamp(path, expected):
    try:
        stamp = json.loads(Path(path).read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError):
        return False
    for k in ("assembly", "genome", "candidates_sha256", "tool", "genetic_code", "threads", "options"):
        if stamp.get(k) != expected.get(k):
            return False
    return True


def preflight_mapping(paths):
    manifest = read_manifest(paths["manifest"])
    selected = read_selected(paths["selected"])
    if set(row["assembly"] for row in manifest) != set(selected):
        raise ValueError(f"Манифест не содержит все {len(selected)} выбранных сборок. Повторите qc.")
    candidates = paths["candidates"]
    if not candidates.is_file():
        raise FileNotFoundError(f"Нет FASTA с кандидатами: {candidates}")
    cand_count = len(fasta_lengths(candidates))
    if cand_count != 51:
        print(f"[WARN] В FASTA {cand_count} последовательностей, ожидалось 51. Продолжаем.", file=sys.stderr)
    for row in manifest:
        fna = Path(row["fasta"])
        if not fna.is_file() or fna.stat().st_size == 0:
            raise FileNotFoundError(f"Нет генома для {row['assembly']}: {fna}")
    return manifest, digest_file(candidates)


def map_genomes(args, paths):
    manifest, candidates_hash = preflight_mapping(paths)
    executable = shutil.which(args.miniprot)
    if not executable:
        raise FileNotFoundError(f"miniprot не найден в PATH: {args.miniprot}")
    out = paths["alignments"]
    out.mkdir(parents=True, exist_ok=True)
    new = already = stale = 0
    if args.assembly:
        manifest = [row for row in manifest if row["assembly"] == args.assembly]
        if len(manifest) != 1:
            raise ValueError(f"Assembly {args.assembly} is missing from the QC manifest")
    for row in manifest:
        acc, fna = row["assembly"], Path(row["fasta"])
        result = out / f"{acc}.gff"
        done = out / f"{acc}.done.json"
        log = out / f"{acc}.log"
        expected = make_stamp(acc, fna, candidates_hash, str(Path(executable).resolve()), args.code, args.threads)
        if not args.force and result.is_file() and result.stat().st_size > 0 and valid_stamp(done, expected):
            already += 1
            continue
        if result.exists() or done.exists():
            stale += 1
        command = [executable, "-T", str(args.code), "-I", "-t", str(args.threads), "--gff", str(fna), str(paths["candidates"])]
        print(f"[MAP] {acc} ({new + 1} новое)  {fna.name}", flush=True)
        fd, tmp_name = tempfile.mkstemp(prefix=f".{acc}.", suffix=".gff.tmp", dir=out)
        try:
            with os.fdopen(fd, "w", encoding="utf-8") as tmp, log.open("w", encoding="utf-8") as err:
                err.write("COMMAND: " + " ".join(command) + "\n")
                err.flush()
                status = subprocess.run(command, stdout=tmp, stderr=err, check=False).returncode
            if status != 0 or os.stat(tmp_name).st_size == 0:
                raise RuntimeError(f"miniprot failed for {acc}, code={status}, log={log}")
            os.replace(tmp_name, result)
            expected["completed_utc"] = datetime.now(timezone.utc).isoformat()
            write_json(done, expected)
        finally:
            if os.path.exists(tmp_name):
                os.unlink(tmp_name)
        new += 1
        if args.limit and new >= args.limit:
            print(f"[STOP] Тестовый лимит: {args.limit}", flush=True)
            break
    print(f"Выполнено новых: {new}; пропущено завершённых: {already}; повторно запущено устаревших: {stale}")
    print(f"GFF: {out} ; для проверки: python scripts/cauris_population_stage.py status")


def parse_attrs(value):
    attrs = {}
    for token in value.split(";"):
        if "=" in token:
            key, item = token.split("=", 1)
            attrs[key] = unquote(item)
    return attrs


def parse_miniprot_gff(path, candidates):
    hits = []
    with open_text(path) as fh:
        for line_no, line in enumerate(fh, 1):
            if not line.strip() or line[0] == "#":
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 9 or fields[2] != "mRNA":
                continue
            attr = parse_attrs(fields[8])
            target = attr.get("Target", "").split()
            if len(target) < 3:
                raise ValueError(f"Нет Target в mRNA: {path}:{line_no}")
            protein = target[0]
            if protein not in candidates:
                continue
            qstart, qend = int(target[1]), int(target[2])
            if qstart < 1 or qend < qstart or qend > candidates[protein]:
                raise ValueError(f"Некорректные координаты белка: {path}:{line_no} Target={attr['Target']}")
            if "Identity" not in attr:
                raise ValueError(f"Нет Identity в mRNA: {path}:{line_no}")
            identity = float(attr["Identity"])
            identity_pct = identity * 100 if identity <= 1 else identity
            coverage_pct = 100 * (qend - qstart + 1) / candidates[protein]
            hits.append({
                "candidate": protein,
                "contig": fields[0],
                "genome_start": int(fields[3]),
                "genome_end": int(fields[4]),
                "strand": fields[6],
                "identity_pct": round(identity_pct, 3),
                "coverage_pct": round(coverage_pct, 3),
                "protein_start": qstart,
                "protein_end": qend,
                "frameshifts": int(attr.get("Frameshift", 0)),
                "stop_codons": int(attr.get("StopCodon", 0)),
            })
    return hits


def summarize(args, paths):
    manifest, candidates_hash = preflight_mapping(paths)
    candidates = fasta_lengths(paths["candidates"])
    out = paths["results"]
    out.mkdir(parents=True, exist_ok=True)
    long_rows, status_lookup = [], {}
    missing = []
    for row in manifest:
        acc = row["assembly"]
        gff = paths["alignments"] / f"{acc}.gff"
        stamp = paths["alignments"] / f"{acc}.done.json"
        if not gff.is_file() or gff.stat().st_size == 0 or not stamp.is_file():
            missing.append(acc)
            continue
        try:
            stamp_data = json.loads(stamp.read_text(encoding="utf-8"))
        except json.JSONDecodeError:
            missing.append(acc)
            continue
        if stamp_data.get("candidates_sha256") != candidates_hash or stamp_data.get("genome") != file_state(row["fasta"]):
            missing.append(acc)
            continue
        matches = defaultdict(list)
        for hit in parse_miniprot_gff(gff, candidates):
            matches[hit["candidate"]].append(hit)
        for cid in candidates:
            candidate_hits = matches[cid]
            strict = [hit for hit in candidate_hits if hit["identity_pct"] >= args.min_identity and hit["coverage_pct"] >= args.min_coverage]
            ranking = lambda hit: (hit["identity_pct"] * hit["coverage_pct"], hit["coverage_pct"], hit["identity_pct"])
            best = max(strict or candidate_hits, key=ranking) if candidate_hits else None
            status = "pass" if strict else ("partial_or_divergent" if candidate_hits else "no_alignment")
            status_lookup[cid, acc] = status
            long_rows.append([
                acc, cid, status, len(candidate_hits), len(strict),
                *( [best["identity_pct"], best["coverage_pct"], best["protein_start"], best["protein_end"],
                    best["contig"], best["genome_start"], best["genome_end"], best["strand"],
                    best["frameshifts"], best["stop_codons"]] if best else [""] * 10),
            ])
    if missing:
        raise RuntimeError(
            f"Не готовы или устарели выравнивания для {len(missing)} сборок, например: {', '.join(missing[:8])}. "
            "Завершите команду map и повторите summarize."
        )
    order = [row["assembly"] for row in manifest]
    columns = ["assembly", "candidate", "status", "n_alignments", "n_strict_alignments",
               "best_identity_pct", "best_coverage_pct", "protein_start", "protein_end", "contig",
               "genome_start", "genome_end", "strand", "frameshifts", "stop_codons"]
    write_tsv(out / "candidate_genome_details.tsv", columns, long_rows)
    matrix_binary, matrix_status, summary = [], [], []
    total = len(order)
    for cid in candidates:
        statuses = [status_lookup[cid, acc] for acc in order]
        npass = statuses.count("pass")
        matrix_binary.append([cid] + [1 if x == "pass" else 0 for x in statuses])
        matrix_status.append([cid] + statuses)
        summary.append([cid, total, npass, statuses.count("partial_or_divergent"), statuses.count("no_alignment"), round(100 * npass / total, 2)])
    summary.sort(key=lambda row: (-row[-1], row[0]))
    write_tsv(out / "candidate_presence_matrix.tsv", ["candidate"] + order, matrix_binary)
    write_tsv(out / "candidate_alignment_status_matrix.tsv", ["candidate"] + order, matrix_status)
    write_tsv(out / "candidate_presence_summary.tsv",
              ["candidate", "n_genomes", "n_pass", "n_partial_or_divergent", "n_no_alignment", "pass_pct"], summary)
    write_json(out / "summary.json", {
        "assemblies": total, "candidates": len(candidates),
        "mapping_status_counts": dict(Counter(x[2] for x in long_rows)),
        "thresholds": {"identity_pct": args.min_identity, "candidate_coverage_pct": args.min_coverage},
        "candidate_fasta": str(paths["candidates"]), "candidate_sha256": candidates_hash,
        "genome_manifest": str(paths["manifest"]),
        "interpretation": "Mapping of protein queries to genome assemblies, not proven orthology, expression, surface exposure or protection",
    })
    print(f"[OK] Обработано геномов: {total}; белков: {len(candidates)}")
    print(f"[OK] Результаты: {out}")
    print("[NOTE] Бинарная матрица: 1 = выравнивание выше заданных порогов, не доказанная ортология.")


def status(args, paths):
    manifest = read_manifest(paths["manifest"])
    out = paths["alignments"]
    current_hash = digest_file(paths["candidates"]) if paths["candidates"].is_file() else None
    count = Counter()
    for r in manifest:
        acc = r["assembly"]
        gff, stamp = out / f"{acc}.gff", out / f"{acc}.done.json"
        if gff.is_file() and gff.stat().st_size > 0 and stamp.is_file():
            try:
                meta = json.loads(stamp.read_text(encoding="utf-8"))
                up_to_date = (meta.get("candidates_sha256") == current_hash and
                              meta.get("genome") == file_state(r["fasta"]))
            except (OSError, json.JSONDecodeError):
                up_to_date = False
            count["done"] += 1 if up_to_date else 0
            count["stale"] += 0 if up_to_date else 1
        elif gff.exists() or (out / f"{acc}.log").exists():
            count["unfinished_or_legacy"] += 1
        else:
            count["pending"] += 1
    print(json.dumps({"expected": len(manifest), **count}, ensure_ascii=False, indent=2))
    print(f"ANI (готовый файл): {paths['ani']} = {'есть' if paths['ani'].is_file() else 'нет'}")


def show_paths(args, paths):
    for key, path in paths.items():
        marker = "[OK]" if path.exists() else "[--]"
        print(f"{marker:4} {key:11} {path}")


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    p.add_argument("--root", help="Корень TEN; по умолчанию каталог выше папки scripts/.")
    sub = p.add_subparsers(dest="command", required=True)
    sub.add_parser("paths", help="Показать, какие файлы скрипт найдёт автоматически")
    q = sub.add_parser("qc", help="Контроль FASTA и обновление манифеста")
    q.add_argument("--selected", help="Заменить путь к selected_accessions.txt")
    q.add_argument("--genomes", help="Заменить путь к ncbi_dataset/data/")
    q.add_argument("--qc", help="Заменить папку результатов QC")
    m = sub.add_parser("map", help="Запустить miniprot на выбранных геномах с автопродолжением")
    m.add_argument("--manifest", help="Путь к QC манифесту")
    m.add_argument("--selected", help="Путь к accession")
    m.add_argument("--candidates", help="FASTA белков-кандидатов")
    m.add_argument("--alignments", help="Каталог GFF")
    m.add_argument("--miniprot", default="miniprot", help="Путь к исполняемому miniprot")
    m.add_argument("--threads", type=int, default=14)
    m.add_argument("--code", type=int, default=12, help="Генетический код CUG-Ser")
    m.add_argument("--limit", type=int, default=0, help="Максимум НОВЫХ сборок: 0 = все")
    m.add_argument("--assembly", help="Map a single assembly from the validated manifest")
    m.add_argument("--force", action="store_true", help="Пересчитать уже готовые GFF")
    t = sub.add_parser("status", help="Посчитать завершённые и ожидающие карты")
    t.add_argument("--manifest", help="Путь к QC манифесту")
    t.add_argument("--alignments", help="Каталог GFF")
    s = sub.add_parser("summarize", help="Получить матрицу и таблицу покрытия")
    s.add_argument("--manifest", help="Путь к QC манифесту")
    s.add_argument("--selected", help="Путь к accession")
    s.add_argument("--candidates", help="FASTA белков-кандидатов")
    s.add_argument("--alignments", help="Каталог GFF")
    s.add_argument("--results", help="Папка результатов")
    s.add_argument("--min-identity", type=float, default=40.0)
    s.add_argument("--min-coverage", type=float, default=70.0)
    args = p.parse_args()
    if getattr(args, "threads", 1) < 1 or getattr(args, "limit", 0) < 0:
        p.error("--threads должен быть >=1, --limit должен быть >=0")
    paths = get_paths(args)
    try:
        if args.command == "paths": show_paths(args, paths)
        elif args.command == "qc": qc(args, paths)
        elif args.command == "map": map_genomes(args, paths)
        elif args.command == "status": status(args, paths)
        elif args.command == "summarize": summarize(args, paths)
    except (ValueError, RuntimeError, FileNotFoundError, OSError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        sys.exit(2)


if __name__ == "__main__":
    main()
