#!/usr/bin/env python3
"""Write a samplesheet pointing at the reads QC produced.

The input sheet is carried through rather than rebuilt: every column it had survives, including
`biome`, which the classification and MAG workflows read and which nothing else in the run knows.
Only `fastq_1` and `fastq_2` are replaced.

A sample that did not come out of KneadData is dropped and named on stderr. Writing a sheet that
points at a file which is not there would only move the failure to whatever runs next.
"""

import argparse
import csv
import sys


def main() -> int:
    parsed = argparse.ArgumentParser(description=__doc__)
    parsed.add_argument("--samplesheet", required=True, help="the validated input sheet")
    parsed.add_argument("--pairs", required=True, help="TSV of sample, mate 1 name, mate 2 name")
    parsed.add_argument("--reads-dir", required=True, help="where the clean reads were published")
    parsed.add_argument("--output", required=True)
    args = parsed.parse_args()

    clean: dict[str, tuple[str, str]] = {}
    with open(args.pairs, encoding="utf-8") as handle:
        for line in handle:
            if not line.strip():
                continue
            sample, mate1, mate2 = line.rstrip("\n").split("\t")
            clean[sample] = (mate1, mate2)

    root = args.reads_dir.rstrip("/")
    with open(args.samplesheet, newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    if not rows:
        print("samplesheet is empty", file=sys.stderr)
        return 1

    kept, missing = [], []
    for row in rows:
        pair = clean.get(row["sample"])
        if pair is None:
            missing.append(row["sample"])
            continue
        row["fastq_1"] = f"{root}/{pair[0]}"
        row["fastq_2"] = f"{root}/{pair[1]}"
        kept.append(row)

    if missing:
        print(f"no clean reads for: {', '.join(sorted(missing))}", file=sys.stderr)
    if not kept:
        print("no sample produced clean reads", file=sys.stderr)
        return 1

    with open(args.output, "w", newline="", encoding="utf-8") as handle:
        # LF, not the csv module's default CRLF: this file is read back by the pipeline
        # and by whatever the user points at it.
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()), lineterminator="\n")
        writer.writeheader()
        writer.writerows(kept)
    return 0


if __name__ == "__main__":
    sys.exit(main())
