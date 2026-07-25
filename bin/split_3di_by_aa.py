#!/usr/bin/env python3
"""Filter a cohort-wide predictions dir to one chunk's protein IDs.

Used inside the PHOLD_COMPARE scatter: SEQKIT_SPLIT2 splits the cohort's
input AA FASTA into M chunks (e.g. 100 k proteins each). Each chunk needs
its own predictions_dir (phold_aa.fasta + phold_3di.fasta) limited to
that chunk's protein IDs, so the per-chunk `phold proteins-compare` call
sees a self-contained input.

Inputs:
    --chunk-aa            Per-chunk AA FASTA (the chunk SEQKIT_SPLIT2 produced).
                          May be .gz. Defines the protein ID set for the chunk.
    --merged-predict-dir  Cohort-wide MERGE_PHOLD_PREDICTIONS output containing
                          at least phold_aa.fasta and phold_3di.fasta.
    --out-dir             Output dir to write the filtered predictions_dir.

The output dir gets:
    phold_aa.fasta                           — AA entries matching the chunk's IDs
    phold_3di.fasta                          — 3Di entries matching the chunk's IDs
    phold_prostT5_3di_mean_probabilities.csv — per-protein prostT5 confidence,
                                               filtered to chunk's IDs. PHOLD compare
                                               reads this to attach `prostt5_confidence`
                                               to the per-CDS TSV — if missing, compare
                                               raises FileNotFoundError.
    phold_prostT5_3di_all_probabilities.json — per-residue probabilities dict,
                                               filtered to chunk's IDs. May be empty
                                               ({}) if PHOLD_PREDICT was run without
                                               --save_per_residue_embeddings.
"""

from __future__ import annotations

import argparse
import gzip
import io
import json
import sys
from pathlib import Path
from typing import IO, Iterator


def _open(path: Path) -> IO[str]:
    if str(path).endswith(".gz"):
        return io.TextIOWrapper(gzip.open(path, "rb"))
    return open(path, "r")


def _iter_records(fh: IO[str]) -> Iterator[tuple[str, list[str]]]:
    """Yield (header_id, [header_line, *seq_lines]) for each FASTA record.

    Memory-light: streams the file; never holds more than one record at a time.
    """
    cur_id: str | None = None
    cur: list[str] = []
    for line in fh:
        if line.startswith(">"):
            if cur_id is not None:
                yield cur_id, cur
            cur = [line]
            # PHOLD's protein IDs contain `::` and `-` but no whitespace; take
            # the first whitespace-separated token after '>' and strip newline.
            cur_id = line[1:].split(None, 1)[0].rstrip()
        else:
            cur.append(line)
    if cur_id is not None:
        yield cur_id, cur


def _collect_ids(aa_path: Path) -> set[str]:
    ids: set[str] = set()
    with _open(aa_path) as fh:
        for line in fh:
            if line.startswith(">"):
                ids.add(line[1:].split(None, 1)[0].rstrip())
    return ids


def _filter_to(src: Path, dst: Path, want: set[str]) -> tuple[int, int]:
    """Copy FASTA records from src to dst whose ID is in `want`. Returns (kept, total)."""
    kept = total = 0
    with _open(src) as fh, open(dst, "w") as out:
        for rid, lines in _iter_records(fh):
            total += 1
            if rid in want:
                out.writelines(lines)
                kept += 1
    return kept, total


def _filter_csv(src: Path, dst: Path, want: set[str]) -> tuple[int, int]:
    """Filter the mean_probabilities CSV to rows whose first column is in `want`.

    PHOLD's CSV has NO header — every line is `<protein_id>,<confidence>`. Keep
    rows whose first comma-separated field is in `want`.
    """
    kept = total = 0
    with open(src, "r") as fh, open(dst, "w") as out:
        for line in fh:
            if not line.strip():
                continue
            total += 1
            rid = line.split(",", 1)[0]
            if rid in want:
                out.write(line)
                kept += 1
    return kept, total


def _filter_json(src: Path, dst: Path, want: set[str]) -> tuple[int, int]:
    """Filter the all_probabilities JSON dict to keys in `want`.

    Loads the whole dict into memory — fine for our scale (mouse: 0-150 MB
    typical when populated; 500-sample: <2 GB; bigger cohorts would warrant
    a streaming parser but aren't on the table yet).
    """
    with open(src, "r") as fh:
        data = json.load(fh)
    if not isinstance(data, dict):
        # PHOLD always emits a dict; bail loudly if not.
        print(f"WARN: {src} is not a JSON dict; writing empty dict", file=sys.stderr)
        filtered: dict = {}
    else:
        filtered = {k: v for k, v in data.items() if k in want}
    with open(dst, "w") as out:
        json.dump(filtered, out)
    return len(filtered), len(data) if isinstance(data, dict) else 0


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--chunk-aa", required=True, type=Path,
                    help="Per-chunk AA FASTA (defines protein ID set for the chunk)")
    ap.add_argument("--merged-predict-dir", required=True, type=Path,
                    help="MERGE_PHOLD_PREDICTIONS output dir (contains phold_aa.fasta + phold_3di.fasta)")
    ap.add_argument("--out-dir", required=True, type=Path,
                    help="Output dir for the filtered predictions")
    args = ap.parse_args()

    args.out_dir.mkdir(parents=True, exist_ok=True)

    aa_src   = args.merged_predict_dir / "phold_aa.fasta"
    di_src   = args.merged_predict_dir / "phold_3di.fasta"
    csv_src  = args.merged_predict_dir / "phold_prostT5_3di_mean_probabilities.csv"
    json_src = args.merged_predict_dir / "phold_prostT5_3di_all_probabilities.json"
    for p in (aa_src, di_src, csv_src, json_src):
        if not p.is_file():
            print(f"ERROR: missing {p}", file=sys.stderr)
            return 1

    want = _collect_ids(args.chunk_aa)
    if not want:
        print(f"ERROR: no protein IDs found in {args.chunk_aa}", file=sys.stderr)
        return 1

    kept_aa,   total_aa   = _filter_to (aa_src,   args.out_dir / "phold_aa.fasta",                            want)
    kept_di,   total_di   = _filter_to (di_src,   args.out_dir / "phold_3di.fasta",                           want)
    kept_csv,  total_csv  = _filter_csv(csv_src,  args.out_dir / "phold_prostT5_3di_mean_probabilities.csv",  want)
    kept_json, total_json = _filter_json(json_src, args.out_dir / "phold_prostT5_3di_all_probabilities.json", want)

    print(f"[split_3di_by_aa] chunk IDs: {len(want)}", file=sys.stderr)
    print(f"[split_3di_by_aa] phold_aa.fasta:                            kept {kept_aa}/{total_aa}", file=sys.stderr)
    print(f"[split_3di_by_aa] phold_3di.fasta:                           kept {kept_di}/{total_di}", file=sys.stderr)
    print(f"[split_3di_by_aa] phold_prostT5_3di_mean_probabilities.csv:  kept {kept_csv}/{total_csv}", file=sys.stderr)
    print(f"[split_3di_by_aa] phold_prostT5_3di_all_probabilities.json:  kept {kept_json}/{total_json}", file=sys.stderr)

    # Sanity: every chunk ID should be findable in AA + 3Di + CSV. JSON may be
    # empty ({}) if PHOLD_PREDICT was run without per-residue embeddings, so
    # only warn about CSV deficits — not JSON.
    missing_aa  = len(want) - kept_aa
    missing_di  = len(want) - kept_di
    missing_csv = len(want) - kept_csv
    if missing_aa or missing_di or missing_csv:
        print(
            f"WARN: chunk has {missing_aa} missing AA, {missing_di} missing 3Di, "
            f"{missing_csv} missing prostT5-confidence rows — PHOLD_COMPARE may skip "
            "those proteins or surface NaN confidence values downstream",
            file=sys.stderr,
        )

    return 0


if __name__ == "__main__":
    sys.exit(main())
