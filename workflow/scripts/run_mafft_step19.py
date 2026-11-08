#!/usr/bin/env python3
"""TEN / C. auris step 19: align unique protein alleles and expand by isolate.

Input:
  09_protein_sequences/by_candidate_unique/<candidate>.faa
  09_protein_sequences/by_candidate/<candidate>.faa
  09_protein_sequences/allele_counts.tsv
  09_protein_sequences/extraction_summary.json
Output:
  10_mafft/alignments_unique/*.aln.faa
  10_mafft/alignments_all/*.aln.faa
  10_mafft/logs/*.log
  10_mafft/alignment_summary.tsv
  10_mafft/mafft_summary.json

MAFFT is run only for >1 unique allele per candidate. One allele is copied
unchanged. Empty candidate files are reported but not aligned. The reconstruction
of isolate-specific alignments is based on exact ungapped sequences, not on
arbitrary output order; all counts and source IDs are validated.

Requires Python >=3.9 (standard library) and MAFFT v7.
"""
from __future__ import annotations
import argparse
import csv
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

VERSION = "1.0.0"
AA_LETTERS = set("ACDEFGHIKLMNPQRSTVWYXUO")
GAPPED_LETTERS = AA_LETTERS | {"-", "."}


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def read_fasta(path: Path, *, allow_gaps=False):
    if not path.is_file():
        raise FileNotFoundError(path)
    records, name, parts = [], None, []
    seen = set()
    def flush():
        if name is None:
            return
        seq = "".join(parts).upper()
        if not seq:
            raise ValueError(f"Empty sequence {name} in {path}")
        allowed = GAPPED_LETTERS if allow_gaps else AA_LETTERS
        invalid = set(seq) - allowed
        if invalid:
            raise ValueError(f"Invalid residue(s) {sorted(invalid)} for {name} in {path}")
        records.append((name, seq))
    with path.open(encoding="utf8") as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                flush()
                name = line[1:].split()[0]
                if not name or name in seen:
                    raise ValueError(f"Missing or duplicate FASTA ID {name!r} in {path}")
                seen.add(name)
                parts = []
            else:
                if name is None:
                    raise ValueError(f"Sequence before FASTA header in {path}")
                parts.append(line)
    flush()
    return records


def write_fasta_atomic(path: Path, records):
    path.parent.mkdir(parents=True, exist_ok=True)
    temp = path.with_name(path.name + ".part")
    with temp.open("w", encoding="utf8") as fh:
        for name, seq in records:
            fh.write(f">{name}\n")
            for i in range(0, len(seq), 80):
                fh.write(seq[i:i+80] + "\n")
    os.replace(temp, path)


def write_tsv_atomic(path: Path, fields, rows):
    path.parent.mkdir(parents=True, exist_ok=True)
    temp = path.with_name(path.name + ".part")
    with temp.open("w", encoding="utf8", newline="") as fh:
        writer = csv.DictWriter(fh, delimiter="\t", fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)
    os.replace(temp, path)


def write_json_atomic(path: Path, data):
    path.parent.mkdir(parents=True, exist_ok=True)
    temp = path.with_name(path.name + ".part")
    temp.write_text(json.dumps(data, indent=2, ensure_ascii=False) + "\n", encoding="utf8")
    os.replace(temp, path)


def read_alleles(path):
    if not path.is_file():
        raise FileNotFoundError(path)
    grouped = defaultdict(dict)
    with path.open(encoding="utf8", newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        expected = {"candidate", "allele", "count", "assemblies", "sequence_sha256"}
        if not reader.fieldnames or not expected <= set(reader.fieldnames):
            raise ValueError(f"Allele table missing {expected - set(reader.fieldnames or [])}")
        for row in reader:
            cand, allele = row["candidate"], row["allele"]
            if allele in grouped[cand]:
                raise ValueError(f"Duplicated allele {cand}:{allele}")
            count = int(row["count"])
            assemblies = [s for s in row["assemblies"].split(",") if s]
            if count < 1 or len(assemblies) != count or len(assemblies) != len(set(assemblies)):
                raise ValueError(f"Allele count/assembly list mismatch {cand}:{allele}")
            grouped[cand][allele] = {"count": count, "assemblies": assemblies, "hash": row["sequence_sha256"]}
    return grouped


def verify_inputs(cid, uniq, all_records, allele_table):
    by_seq = {seq: name for name, seq in uniq}
    if len(by_seq) != len(uniq):
        raise ValueError(f"{cid}: repeated sequence in by_candidate_unique, run step18 again")
    count_by_seq = Counter(seq for _, seq in all_records)
    if set(by_seq) != set(count_by_seq):
        raise ValueError(f"{cid}: unique variants and full sequences differ")
    if len(uniq) != len(allele_table):
        raise ValueError(f"{cid}: {len(uniq)} unique FASTA vs {len(allele_table)} allele table entries")
    if any(len(seq) < 20 for _, seq in uniq):
        raise ValueError(f"{cid}: unexpected protein under 20 amino acids")
    all_ids = set()
    for name, seq in all_records:
        prefix = f"{cid}|"
        if not name.startswith(prefix):
            raise ValueError(f"{cid}: full sequence ID {name!r} does not begin {prefix}")
        asm = name[len(prefix):]
        if asm in all_ids:
            raise ValueError(f"{cid}: duplicate assembly {asm}")
        all_ids.add(asm)
    expected_all_assemblies = set()
    for name, seq in uniq:
        prefix = f"{cid}|"
        if not name.startswith(prefix):
            raise ValueError(f"{cid}: unique ID {name!r} missing prefix")
        m = re.fullmatch(r"(allele_\d{3,})\|n=(\d+)", name[len(prefix):])
        if not m:
            raise ValueError(f"{cid}: unexpected allele ID {name!r}")
        allele, n = m.group(1), int(m.group(2))
        a = allele_table.get(allele)
        if a is None or a["count"] != n or count_by_seq[seq] != n:
            raise ValueError(f"{cid}: count mismatch for {allele}")
        if a["hash"] != hashlib.sha256(seq.encode()).hexdigest():
            raise ValueError(f"{cid}: sequence SHA mismatch for {allele}")
        expected_all_assemblies.update(a["assemblies"])
        for asm in a["assemblies"]:
            key = f"{cid}|{asm}"
            if not any(nm == key and sq == seq for nm, sq in all_records):
                raise ValueError(f"{cid}: {allele} assembly {asm} mismatches original FASTA")
    if expected_all_assemblies != all_ids:
        raise ValueError(f"{cid}: allele assemblies do not match full FASTA")


def validate_alignment(cid, aligned, uniq):
    input_by_id = dict(uniq)
    if len(aligned) != len(uniq) or {n for n,_ in aligned} != set(input_by_id):
        raise ValueError(f"{cid}: MAFFT changed number or names of sequences")
    lengths = {len(seq) for _, seq in aligned}
    if len(lengths) != 1:
        raise ValueError(f"{cid}: unequal alignment column counts: {lengths}")
    ncols = next(iter(lengths), 0)
    for name, seq in aligned:
        if seq.replace("-", "").replace(".", "") != input_by_id[name]:
            raise ValueError(f"{cid}: MAFFT altered ungapped sequence for {name}")
    return ncols


def get_mafft_version(executable):
    p = subprocess.run([executable, "--version"], stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    return (p.stdout or p.stderr).strip().splitlines()[:2]


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--root", type=Path, default=Path(__file__).resolve().parent.parent,
                    help="Project root (default parent of scripts/)")
    ap.add_argument("--input", type=Path, help="Default <root>/09_protein_sequences")
    ap.add_argument("--out", type=Path, help="Default <root>/10_mafft")
    ap.add_argument("--mafft", default="mafft", help="MAFFT executable or path")
    ap.add_argument("--threads", type=int, default=14)
    ap.add_argument("--method", choices=["auto", "linsi"], default="auto")
    ap.add_argument("--limit", type=int, default=0, help="Process at most N candidates this run; 0=all")
    ap.add_argument("--force", action="store_true", help="Recompute aligned candidates even if cached")
    args = ap.parse_args(argv)
    if args.threads < 1 or args.limit < 0:
        ap.error("--threads >= 1; --limit >= 0")
    root = args.root.expanduser().resolve()
    inp = (args.input or root / "09_protein_sequences").expanduser().resolve()
    out = (args.out or root / "10_mafft").expanduser().resolve()
    unique_dir, all_dir = inp / "by_candidate_unique", inp / "by_candidate"
    if not unique_dir.is_dir() or not all_dir.is_dir():
        raise FileNotFoundError(f"Expected step18 FASTA directories: {unique_dir}, {all_dir}")
    extraction_path = inp / "extraction_summary.json"
    if not extraction_path.is_file():
        raise FileNotFoundError(f"Step18 extraction_summary.json missing: {extraction_path}")
    extraction = json.loads(extraction_path.read_text(encoding="utf8"))
    alleles = read_alleles(inp / "allele_counts.tsv")
    uniq_files = sorted(unique_dir.glob("*.faa"))
    all_files = {p.stem: p for p in all_dir.glob("*.faa")}
    if not uniq_files or set(p.stem for p in uniq_files) != set(all_files):
        raise ValueError("Input FASTA sets differ or are empty")
    expected = int(extraction.get("n_candidates", 0))
    if expected and len(uniq_files) != expected:
        raise ValueError(f"Expected {expected} candidate FASTAs from step18, found {len(uniq_files)}")
    mafft = shutil.which(args.mafft)
    if mafft is None:
        raise FileNotFoundError(f"MAFFT missing from PATH: {args.mafft}")
    mafft = str(Path(mafft).resolve())
    version = get_mafft_version(mafft)
    for p in (out / "alignments_unique", out / "alignments_all", out / "logs", out / "cache"):
        p.mkdir(parents=True, exist_ok=True)
    summary_rows = []
    processed_now = cached_now = 0
    total_records, total_uniq = 0, 0
    for j, uniq_file in enumerate(uniq_files, 1):
        cid = uniq_file.stem
        full_file = all_files[cid]
        uniq = read_fasta(uniq_file)
        originals = read_fasta(full_file)
        verify_inputs(cid, uniq, originals, alleles.get(cid, {}))
        total_records += len(originals)
        total_uniq += len(uniq)
        uout = out / "alignments_unique" / f"{cid}.aln.faa"
        aout = out / "alignments_all" / f"{cid}.aln.faa"
        cache_path = out / "cache" / f"{cid}.json"
        log_path = out / "logs" / f"{cid}.log"
        stamp = {"pipeline_version": VERSION, "unique_fasta_sha256": sha256(uniq_file),
                 "full_fasta_sha256": sha256(full_file), "allele_counts_sha256": sha256(inp / "allele_counts.tsv"),
                 "mafft": mafft, "mafft_version": version, "method": args.method, "threads": args.threads}
        if len(uniq) == 0:
            status, ncols, from_cache = "no_extracted_sequences", 0, False
            # Do not emit misleading empty alignments; delete results from a previous run.
            for stale in (uout, aout, cache_path):
                stale.unlink(missing_ok=True)
        else:
            from_cache = False
            if not args.force and cache_path.is_file() and uout.is_file() and aout.is_file():
                try:
                    cache = json.loads(cache_path.read_text(encoding="utf8"))
                    if cache.get("stamp") == stamp and cache.get("sha256_unique_alignment") == sha256(uout) and cache.get("sha256_all_alignment") == sha256(aout):
                        from_cache = True
                        status = cache["status"]
                        ncols = cache["n_columns"]
                except (ValueError, OSError, KeyError):
                    pass
            if not from_cache:
                if args.limit and processed_now >= args.limit:
                    print(f"[LIMIT] Processed {processed_now} candidates; resume without --limit", flush=True)
                    break
                if len(uniq) == 1:
                    aligned = uniq
                    status = "single_unique_allele"
                    log_path.write_text("Only one unique sequence; no MAFFT alignment necessary.\n", encoding="utf8")
                else:
                    cmd = [mafft, "--amino", "--anysymbol", "--thread", str(args.threads)]
                    if args.method == "auto":
                        cmd += ["--auto"]
                    else:
                        cmd += ["--localpair", "--maxiterate", "1000"]
                    cmd += [str(uniq_file)]
                    with tempfile.TemporaryDirectory(prefix="mafft_", dir=out) as tmp:
                        tmp_out = Path(tmp) / "aligned.faa"
                        with tmp_out.open("w", encoding="utf8") as aln_fh, log_path.open("w", encoding="utf8") as log_fh:
                            log_fh.write("COMMAND: " + " ".join(cmd) + "\n")
                            log_fh.flush()
                            p = subprocess.run(cmd, stdout=aln_fh, stderr=log_fh, text=True, encoding="utf8")
                        if p.returncode:
                            raise RuntimeError(f"MAFFT failed for {cid}: exit {p.returncode}; see {log_path}")
                        aligned = read_fasta(tmp_out, allow_gaps=True)
                    status = "aligned"
                ncols = validate_alignment(cid, aligned, uniq)
                seq_to_alignment = {dict(uniq)[name]: aln for name, aln in aligned}
                expanded = [(name, seq_to_alignment[seq]) for name, seq in originals]
                if len(expanded) != len(originals):
                    raise ValueError(f"{cid}: allele expansion mismatch")
                write_fasta_atomic(uout, aligned)
                write_fasta_atomic(aout, expanded)
                write_json_atomic(cache_path, {
                    "stamp": stamp, "status": status, "n_columns": ncols,
                    "n_unique": len(uniq), "n_assemblies": len(originals),
                    "sha256_unique_alignment": sha256(uout),
                    "sha256_all_alignment": sha256(aout)
                })
                processed_now += 1
                print(f"[RUN] {j}/{len(uniq_files)} {cid}: {status}, variants={len(uniq)}, genomes={len(originals)}, columns={ncols}", flush=True)
            else:
                cached_now += 1
                print(f"[CACHED] {j}/{len(uniq_files)} {cid}", flush=True)
        summary_rows.append({"candidate": cid, "status": status, "n_genomes": len(originals),
                             "n_unique_alleles": len(uniq), "alignment_columns": ncols,
                             "max_length_aa": max((len(seq) for _, seq in originals), default=0),
                             "min_length_aa": min((len(seq) for _, seq in originals), default=0),
                             "source_unique_sha256": sha256(uniq_file),
                             "source_all_sha256": sha256(full_file)})
    fields = ["candidate", "status", "n_genomes", "n_unique_alleles", "alignment_columns",
              "max_length_aa", "min_length_aa", "source_unique_sha256", "source_all_sha256"]
    write_tsv_atomic(out / "alignment_summary.tsv", fields, summary_rows)
    run_summary = {
        "timestamp_utc": datetime.now(timezone.utc).isoformat(), "version": VERSION,
        "root": str(root), "input": str(inp), "out": str(out),
        "n_candidate_fasta_files": len(uniq_files), "n_candidates_reported": len(summary_rows),
        "n_candidates_aligned": sum(r["status"] == "aligned" for r in summary_rows),
        "n_candidates_single_unique": sum(r["status"] == "single_unique_allele" for r in summary_rows),
        "n_candidates_no_extracted_sequences": sum(r["status"] == "no_extracted_sequences" for r in summary_rows),
        "n_isolate_sequences_in_processed_candidates": sum(r["n_genomes"] for r in summary_rows),
        "n_unique_alleles_in_processed_candidates": sum(r["n_unique_alleles"] for r in summary_rows),
        "processed_this_run": processed_now, "reused_cache_this_run": cached_now,
        "complete": len(summary_rows) == len(uniq_files), "mafft_version": version,
        "method": args.method, "threads": args.threads,
        "note": "Alignments only; no inference of conservation or phylogeny. Alleles expanded exactly to original per-genome FASTA."
    }
    write_json_atomic(out / "mafft_summary.json", run_summary)
    print("[SUMMARY] " + json.dumps({k:run_summary[k] for k in ("n_candidates_reported", "n_candidates_aligned", "n_candidates_single_unique", "n_candidates_no_extracted_sequences", "complete")}, ensure_ascii=False), flush=True)
    return 0


if __name__ == "__main__":
    try:
        sys.exit(main())
    except (OSError, ValueError, RuntimeError, subprocess.SubprocessError) as exc:
        print("[ERROR] " + str(exc), file=sys.stderr)
        sys.exit(1)
