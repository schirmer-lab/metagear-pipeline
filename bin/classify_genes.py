#!/usr/bin/env python
"""Classify each gene-catalog representative by the contig classes of its members.

Reads per-sample per-contig TSVs from classification (contig_id →
primary_class), streams the gene cluster TSV (rep<TAB>member), and emits one
row per representative with the set of primary_classes seen across members,
per-class counts, and the representative's own contig class.

Output schema (TSV):

    rep_id        cluster representative gene id
    classes       comma-separated set of *non-unknown* primary_classes seen
                  across members, sorted alphabetically. `unknown` is dropped
                  from this column when at least one real class is present —
                  members from unclassified contigs reflect incomplete
                  coverage, not a biological feature of the cluster. When
                  ALL members are unknown the column shows `unknown` so the
                  all-unclassified case stays visible.
    num_members   total members in the cluster
    class_counts  e.g. `bacteria=10,unknown=5,virus=2` — sorted alphabetically.
                  ALWAYS includes the unknown count when present, so coverage
                  can be reconstructed independently of `classes`.
    rep_class     the representative's own contig class (can be `unknown`)
    multi_class   `yes` if ≥2 distinct *non-unknown* classes are seen across
                  members. Use this to surface shared/mixed clusters with
                  real cross-class biology (HGT, AMGs, mis-split proviruses)
                  without bacteria+unknown coverage-noise dominating the
                  signal. Filter on `class_counts !~ /unknown=/` to find
                  clusters where every member has a real classification.

This is a clean re-derivation from per-contig signal; it does NOT merge with
virus's `.draft` aggregated TSV. The per-contig TSV's `primary_class`
is treated as authoritative.

Gene IDs are assumed to be `<contig>::<idx>::<start>::<stop>::<strand>` so the
contig is everything before the first `::`. This matches Prodigal output as
emitted by the metagear gene-call subworkflow.
"""

from __future__ import annotations

import argparse
import csv
import logging
import sys
from collections import defaultdict
from pathlib import Path
from typing import Dict, Iterable, List

LOG = logging.getLogger("classify_genes")

OUTPUT_HEADER = [
    "rep_id",
    "classes",
    "num_members",
    "class_counts",
    "rep_class",
    "multi_class",
]

UNKNOWN = "unknown"


def _load_per_contig_tsvs(paths: Iterable[Path]) -> Dict[str, str]:
    """Load all per-sample per-contig TSVs into a single contig→class map.

    Skips zero-byte or missing files. Reads `contig_id` and `primary_class`
    by column name so additions to the per-contig TSV schema don't break us.
    """
    out: Dict[str, str] = {}
    for p in paths:
        if not p.exists() or p.stat().st_size == 0:
            LOG.warning("skipping empty/missing per-contig TSV: %s", p)
            continue
        with p.open() as fh:
            reader = csv.DictReader(fh, delimiter="\t")
            fields = reader.fieldnames or []
            if "contig_id" not in fields or "primary_class" not in fields:
                LOG.warning("skipping %s: missing contig_id or primary_class column", p)
                continue
            n = 0
            for row in reader:
                cid = row["contig_id"]
                if cid:
                    out[cid] = row["primary_class"] or UNKNOWN
                    n += 1
            LOG.info("loaded %d contig→class rows from %s", n, p.name)
    LOG.info("total distinct contigs in per-contig map: %d", len(out))
    return out


def _contig_from_gene_id(gene_id: str) -> str:
    """Gene ID format: `<contig>::<idx>::<start>::<stop>::<strand>`.

    Returns the contig prefix (everything before the first `::`).
    """
    sep = gene_id.find("::")
    return gene_id[:sep] if sep > 0 else gene_id


def classify(contig_class: Dict[str, str], clusters_tsv: Path, out_path: Path) -> None:
    members: Dict[str, int] = defaultdict(int)
    class_counts: Dict[str, Dict[str, int]] = defaultdict(lambda: defaultdict(int))
    rep_class: Dict[str, str] = {}
    rep_order: List[str] = []
    seen_reps: set[str] = set()

    n_rows = 0
    n_members_without_contig_in_map = 0

    with clusters_tsv.open() as fh:
        for raw in fh:
            line = raw.rstrip("\n")
            if not line:
                continue
            parts = line.split("\t")
            if len(parts) < 2:
                continue
            rep, member = parts[0], parts[1]

            if rep not in seen_reps:
                seen_reps.add(rep)
                rep_order.append(rep)
                rep_contig = _contig_from_gene_id(rep)
                rep_class[rep] = contig_class.get(rep_contig, UNKNOWN)

            member_contig = _contig_from_gene_id(member)
            if member_contig in contig_class:
                cls = contig_class[member_contig]
            else:
                cls = UNKNOWN
                n_members_without_contig_in_map += 1
            class_counts[rep][cls] += 1
            members[rep] += 1
            n_rows += 1

    LOG.info("streamed %d cluster rows; %d members had no entry in the per-contig map (→ unknown)",
             n_rows, n_members_without_contig_in_map)
    LOG.info("distinct representatives: %d", len(rep_order))

    multi_class_count = 0
    with out_path.open("wt", newline="") as out_fh:
        writer = csv.writer(out_fh, delimiter="\t", lineterminator="\n")
        writer.writerow(OUTPUT_HEADER)
        for rep in rep_order:
            counts = class_counts[rep]
            # Drop `unknown` from the classes column when other classes are
            # present — "bacteria,unknown" reflects incomplete coverage, not
            # biologically interesting cross-class variation. Keep the
            # all-unknown case as `unknown` so it stays visible.
            non_unknown = {k for k in counts.keys() if k != UNKNOWN}
            if non_unknown:
                classes_str = ",".join(sorted(non_unknown))
            else:
                classes_str = UNKNOWN
            counts_str = ",".join(f"{k}={v}" for k, v in sorted(counts.items()))
            multi = "yes" if len(non_unknown) > 1 else "no"
            if multi == "yes":
                multi_class_count += 1
            writer.writerow([
                rep,
                classes_str,
                members[rep],
                counts_str,
                rep_class[rep],
                multi,
            ])

    LOG.info("wrote %d rows to %s (%d multi_class)", len(rep_order), out_path, multi_class_count)


def _parse_args(argv: Iterable[str] | None = None) -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__.split("\n", 1)[0])
    p.add_argument("--per-contig-tsvs", required=True, nargs="+", type=Path,
                   help="One or more per-sample per-contig TSVs from classification.")
    p.add_argument("--clusters-tsv", required=True, type=Path,
                   help="Full gene cluster TSV (rep<TAB>member).")
    p.add_argument("--output", required=True, type=Path)
    p.add_argument("--log-level", default="INFO", choices=["DEBUG", "INFO", "WARNING", "ERROR"])
    return p.parse_args(list(argv) if argv is not None else None)


def main(argv: Iterable[str] | None = None) -> int:
    args = _parse_args(argv)
    logging.basicConfig(level=args.log_level, format="%(asctime)s | %(levelname)s | %(message)s")
    contig_class = _load_per_contig_tsvs(args.per_contig_tsvs)
    if not contig_class:
        LOG.error("no contig→class entries loaded; refusing to write empty classification")
        return 1
    classify(contig_class, args.clusters_tsv, args.output)
    return 0


if __name__ == "__main__":
    sys.exit(main())
