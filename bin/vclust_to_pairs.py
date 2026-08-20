#!/usr/bin/env python
"""Convert Vclust cluster assignments into the pipeline's clustering interface.

Vclust writes `<sequence>\\t<cluster_id>` with a numeric cluster id. Everything
downstream of CLUSTER_SEQUENCES consumes MMseqs2 easy-cluster's shape instead: a
`representative<TAB>member` table plus a FASTA of representatives. This script
converts one into the other so the criterion can change without touching any
consumer.

The representative of a cluster is its longest member, which is what CheckV's
aniclust does; the sequence name breaks ties so the choice is deterministic
across runs and across machines.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import sys


def read_lengths(fasta: str) -> tuple[dict[str, int], list[str]]:
    lengths: dict[str, int] = {}
    order: list[str] = []
    name = None
    total = 0
    with open(fasta) as handle:
        for line in handle:
            if line.startswith(">"):
                if name is not None:
                    lengths[name] = total
                name = line[1:].split()[0]
                order.append(name)
                total = 0
            else:
                total += len(line.strip())
    if name is not None:
        lengths[name] = total
    return lengths, order


def read_clusters(path: str) -> dict[str, list[str]]:
    members: dict[str, list[str]] = {}
    with open(path) as handle:
        reader = csv.reader(handle, delimiter="\t")
        first = next(reader, None)
        if first is not None and first[:1] != ["object"]:
            # No header: the first row is data.
            handle.seek(0)
            reader = csv.reader(handle, delimiter="\t")
        for row in reader:
            if len(row) < 2:
                continue
            members.setdefault(row[1], []).append(row[0])
    return members


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--fasta", required=True, help="the clustered input FASTA (uncompressed)")
    parser.add_argument("--clusters", required=True, help="vclust cluster output")
    parser.add_argument("--out-pairs", required=True, help="representative<TAB>member TSV to write")
    parser.add_argument("--out-representatives", required=True, help="gzipped representative FASTA to write")
    args = parser.parse_args(argv)

    lengths, order = read_lengths(args.fasta)
    members = read_clusters(args.clusters)

    assigned = sum(len(v) for v in members.values())
    if assigned != len(lengths):
        sys.exit(
            f"vclust assigned {assigned} sequences but the input has {len(lengths)}; "
            "refusing to write a partial catalog"
        )

    representatives: dict[str, list[str]] = {}
    for ids in members.values():
        ids.sort(key=lambda i: (-lengths[i], i))
        representatives[ids[0]] = ids
    if len(representatives) != len(members):
        sys.exit("two clusters selected the same representative; input names are not unique")

    # newline="\n" matters: csv.writer defaults to \r\n, and the clusters table is
    # read downstream by awk, cut and pandas on Unix. A trailing \r silently
    # becomes part of the last field and breaks every join against it.
    with open(args.out_pairs, "w", newline="\n") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        # Emit in input order so the table is stable between runs.
        for name in order:
            if name in representatives:
                for member in representatives[name]:
                    writer.writerow([name, member])

    keep = set(representatives)
    written = 0
    with open(args.fasta) as src, gzip.open(args.out_representatives, "wt") as dst:
        emit = False
        for line in src:
            if line.startswith(">"):
                emit = line[1:].split()[0] in keep
                written += emit
            if emit:
                dst.write(line)
    if written != len(keep):
        sys.exit(f"wrote {written} of {len(keep)} representatives")

    print(f"clusters: {len(members)}  representatives: {len(keep)}  sequences: {assigned}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
