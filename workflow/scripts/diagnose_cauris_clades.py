#!/usr/bin/env python3
"""Audit nearest-reference ANI clade decisions without re-running skani or changing results.

Usage:
    python scripts/diagnose_cauris_clades.py
"""
from __future__ import annotations
import argparse
import csv
import json
from collections import Counter, defaultdict
from pathlib import Path
from statistics import median


def tsv(path):
    with path.open(encoding="utf-8", newline="") as f:
        return list(csv.DictReader(f, delimiter="\t"))


def value(s):
    try:
        return float(s)
    except (ValueError, TypeError):
        return None


def fmt(v, nd=4):
    return f"{v:.{nd}f}" if v is not None else "NA"


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--folder", type=Path,
                   default=Path("08_population/results/clade_typing"))
    p.add_argument("--min-af", type=float, default=70.0)
    args = p.parse_args()
    folder = args.folder
    comparisons = tsv(folder / "all_reference_comparisons.tsv")
    assignments = tsv(folder / "clade_assignments.tsv")
    assignment_by_id = {r["assembly"]: r for r in assignments}
    grouped = defaultdict(list)
    for r in comparisons:
        ani = value(r.get("ANI_pct"))
        afr = value(r.get("AF_ref_pct"))
        afq = value(r.get("AF_query_pct"))
        if ani is None or afr is None or afq is None:
            continue
        grouped[r["assembly"]].append({
            "clade": r["ref_clade"], "reference": r["reference_accession"],
            "ANI": ani, "AFref": afr, "AFquery": afq,
            "valid": min(afr, afq) >= args.min_af,
        })

    audits = []
    for assembly in sorted(assignment_by_id):
        row = assignment_by_id[assembly]
        hits = sorted((h for h in grouped.get(assembly, []) if h["valid"]),
                      key=lambda h: h["ANI"], reverse=True)
        top = hits[0] if hits else None
        nxt = hits[1] if len(hits) > 1 else None
        gap = top["ANI"] - nxt["ANI"] if nxt else None
        expected = row.get("known_control_clade", "").strip()
        audits.append({
            "assembly": assembly, "control_expected": expected,
            "nearest": top["clade"] if top else "",
            "nearest_reference": top["reference"] if top else "",
            "nearest_ANI": fmt(top["ANI"] if top else None),
            "nearest_min_AF": fmt(min(top["AFref"],top["AFquery"]) if top else None,2),
            "second": nxt["clade"] if nxt else "",
            "second_ANI": fmt(nxt["ANI"] if nxt else None),
            "margin_pp": fmt(gap,5),
            "nearest_matches_control": ("YES" if expected and top and top["clade"]==expected else
                                        "NO" if expected else ""),
            "old_status": row.get("assignment_status", ""),
        })

    out_path = folder / "clade_diagnostics.tsv"
    with out_path.open("w", encoding="utf-8", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(audits[0]) if audits else ["assembly"],delimiter="\t")
        writer.writeheader()
        writer.writerows(audits)

    controls = [r for r in audits if r["control_expected"]]
    margin_vals = [value(r["margin_pp"]) for r in audits]
    margin_vals = [v for v in margin_vals if v is not None]
    nearest_counts = Counter(r["nearest"] or "No_valid_AF" for r in audits)
    ctrl_ok = sum(r["nearest_matches_control"] == "YES" for r in controls)
    print("=== Сводка первичной диагностики ===")
    print(f"Геномов в отчёте: {len(audits)}; сравнений в исходном TSV: {len(comparisons)}")
    print(f"Ближайшие референсы (без порога отрыва): {dict(sorted(nearest_counts.items()))}")
    print(f"Известные контроли: {ctrl_ok}/{len(controls)} имеют правильный ближайший референс")
    if margin_vals:
        print(f"Отрыв ANI, п.п.: min {min(margin_vals):.5f}; median {median(margin_vals):.5f}; max {max(margin_vals):.5f}")
        for cutoff in (0, .02, .05, .1, .2, .5):
            print(f"  Отрыв ≥{cutoff:.2f}: {sum(x>=cutoff for x in margin_vals)}/{len(margin_vals)}")
    print("\n=== Контрольные геномы ===")
    print("assembly\texpected\tnearest\tbest_ANI\tsecond\tsecond_ANI\tmargin_pp\tmatch")
    for r in controls:
        print("\t".join(r[k] for k in (
            "assembly","control_expected","nearest","nearest_ANI",
            "second","second_ANI","margin_pp","nearest_matches_control")))
    print(f"\nПолный диагностический отчёт: {out_path}")
    if controls and ctrl_ok != len(controls):
        print("WARNING: Есть контроли с неправильным ближайшим референсом. Пороги менять НЕЛЬЗЯ до проверки референсов/файлов.")
    elif controls:
        print("Все известные контроли имеют правильный ближайший референс. Неоднозначность связана с выбранным порогом и требует оценки величины отрывов.")

if __name__ == "__main__":
    main()
