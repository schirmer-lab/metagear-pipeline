#!/usr/bin/env python
"""Enrich the per-sample contigs.tsv from integrated_classification with
cohort-level MAG taxonomy and cluster IDs from cohort_dereplication.

For each contig that was bin-attributable in v1 (classifier == 'binette'),
we look up:

    1. which Binette bin owns the contig    (from final_contig_to_bin.tsv)
    2. which dRep secondary cluster that bin belongs to    (from Cdb.csv)
    3. which genome is the cluster winner    (from Wdb.csv)
    4. which GTDB lineage the winner received    (from gtdbtk.*.summary.tsv)

The resulting columns are appended:

    cluster_id       e.g. "1_2"     (blank for non-bin contigs)
    taxonomy_source  'gtdb-tk' | 'mmseqs2' | ''       (provenance of lineage)
    lineage          existing column; now also filled for bin contigs

Non-bin contigs (virus / plasmid / unbinned-but-mmseqs-classified) pass
through unchanged except for the new columns sitting blank (or with
taxonomy_source='mmseqs2' when the v1 row already had a mmseqs lineage).
"""

from __future__ import annotations

import argparse
import csv
import logging
import sys
from pathlib import Path
from typing import Dict

LOG = logging.getLogger("enrich_per_contig_tsv")


def _load_cdb(cdb_csv: Path) -> Dict[str, str]:
    """genome (filename) -> secondary_cluster ('1_2'). Cdb.csv is dRep's
    cluster table. The `genome` column carries the full FASTA filename."""
    out: Dict[str, str] = {}
    with cdb_csv.open() as fh:
        reader = csv.DictReader(fh)
        for row in reader:
            out[row["genome"]] = row["secondary_cluster"]
    LOG.info("loaded %d genome→cluster mappings from %s", len(out), cdb_csv.name)
    return out


def _load_wdb(wdb_csv: Path) -> Dict[str, str]:
    """secondary_cluster -> winning genome filename."""
    out: Dict[str, str] = {}
    with wdb_csv.open() as fh:
        reader = csv.DictReader(fh)
        for row in reader:
            out[row["cluster"]] = row["genome"]
    LOG.info("loaded %d cluster→winner mappings from %s", len(out), wdb_csv.name)
    return out


def _load_gtdb(gtdb_dir: Path) -> Dict[str, str]:
    """Combine all gtdbtk.*.summary.tsv files into one map. GTDB-Tk emits
    per-domain summaries (gtdbtk.bac120.summary.tsv, gtdbtk.ar53.summary.tsv);
    both have columns `user_genome` and `classification`. user_genome is the
    genome name WITHOUT the .fa suffix."""
    out: Dict[str, str] = {}
    if not gtdb_dir.exists():
        LOG.warning("GTDB-Tk dir %s missing; lineage will be unfilled", gtdb_dir)
        return out
    for tsv in sorted(gtdb_dir.glob("*summary.tsv")):
        with tsv.open() as fh:
            reader = csv.DictReader(fh, delimiter="\t")
            if "user_genome" not in (reader.fieldnames or []):
                LOG.warning("skipping %s: no user_genome column", tsv.name)
                continue
            for row in reader:
                out[row["user_genome"]] = row.get("classification", "")
    LOG.info("loaded %d GTDB-Tk classifications", len(out))
    return out


def _load_contig_to_bin(c2b_tsv: Path, sample_id: str) -> Dict[str, str]:
    """Binette's final_contig_to_bin.tsv: contig_name<TAB>bin_name.
    Returns contig_name -> renamed_bin (`<sample>.<binette_bin_name>.fa`)
    so the value can be joined directly against Cdb's `genome` column."""
    out: Dict[str, str] = {}
    if not c2b_tsv.exists():
        LOG.warning("contig→bin %s missing; no bin attributions for sample", c2b_tsv)
        return out
    with c2b_tsv.open() as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) < 2:
                continue
            contig, bin_name = parts[0], parts[1]
            # PREPARE_DREP_INPUTS renames each bin to `<sample>.<orig>.fa`.
            out[contig] = f"{sample_id}.{bin_name}.fa"
    LOG.info("loaded %d contig→bin mappings for sample %s", len(out), sample_id)
    return out


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--contigs", required=True, type=Path,
                        help="v1 per-sample contigs.tsv")
    parser.add_argument("--contig-to-bin", required=True, type=Path,
                        help="Binette final_contig_to_bin.tsv for the same sample")
    parser.add_argument("--cdb", required=True, type=Path,
                        help="dRep Cdb.csv (cohort cluster table)")
    parser.add_argument("--wdb", required=True, type=Path,
                        help="dRep Wdb.csv (cluster winners)")
    parser.add_argument("--gtdb-dir", required=True, type=Path,
                        help="Directory with gtdbtk.*.summary.tsv files")
    parser.add_argument("--sample", required=True,
                        help="Sample ID (must match the prefix used in PREPARE_DREP_INPUTS)")
    parser.add_argument("--out", required=True, type=Path)
    parser.add_argument("--log-level", default="INFO")
    args = parser.parse_args()

    logging.basicConfig(level=args.log_level,
                        format="%(asctime)s %(levelname)s %(message)s")

    cdb = _load_cdb(args.cdb)
    wdb = _load_wdb(args.wdb)
    gtdb = _load_gtdb(args.gtdb_dir)
    contig2bin = _load_contig_to_bin(args.contig_to_bin, args.sample)

    n_rows = 0
    n_enriched = 0
    n_lineage_filled = 0

    with args.contigs.open() as in_fh, args.out.open("w") as out_fh:
        reader = csv.DictReader(in_fh, delimiter="\t")
        if reader.fieldnames is None:
            LOG.error("input %s has no header", args.contigs)
            return 1
        fieldnames = list(reader.fieldnames)
        for col in ("cluster_id", "taxonomy_source"):
            if col not in fieldnames:
                fieldnames.append(col)
        if "lineage" not in fieldnames:
            fieldnames.append("lineage")

        writer = csv.DictWriter(out_fh, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()

        for row in reader:
            n_rows += 1
            contig = row["contig_id"]
            bin_name = contig2bin.get(contig)
            cluster = cdb.get(bin_name) if bin_name else None
            winner = wdb.get(cluster) if cluster else None
            # winner is "<sample>.binette_binN.fa"; GTDB-Tk's user_genome
            # is the same name minus the .fa suffix.
            winner_key = winner[:-3] if winner and winner.endswith(".fa") else winner
            gtdb_lineage = gtdb.get(winner_key) if winner_key else None

            row["cluster_id"] = cluster or ""
            if gtdb_lineage:
                row["lineage"] = gtdb_lineage
                row["taxonomy_source"] = "gtdb-tk"
                n_lineage_filled += 1
            elif row.get("lineage"):
                # Preserve the v1 mmseqs lineage if present and provenance was
                # not yet set. v1 didn't tag taxonomy_source so we infer it.
                row.setdefault("taxonomy_source", "mmseqs2")
            else:
                row.setdefault("taxonomy_source", "")

            if cluster:
                n_enriched += 1
            writer.writerow(row)

    LOG.info("processed %d rows; %d gained cluster_id; %d gained gtdb-tk lineage",
             n_rows, n_enriched, n_lineage_filled)
    return 0


if __name__ == "__main__":
    sys.exit(main())
