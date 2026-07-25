#!/usr/bin/env python
"""Build a cohort-level mag_catalog.csv summarizing every dRep cluster
together with its winning genome's QC + GTDB-Tk taxonomy.

Columns (in order):

    cluster_id          dRep secondary_cluster ("1_2", "3_1", …)
    winning_genome      Wdb winner, e.g. SAMPLE-0.binette_bin3.fa
    n_members           number of genomes assigned to the cluster
    n_member_samples    number of distinct samples those genomes came from
    member_genomes      semicolon-separated full list (Cdb membership)
    completeness        winner's completeness from genomeInfo.csv
    contamination       winner's contamination from genomeInfo.csv
    gtdb_lineage        winner's lineage from gtdbtk.*.summary.tsv, or blank
"""

from __future__ import annotations

import argparse
import csv
import logging
import re
import sys
from collections import defaultdict
from pathlib import Path
from typing import Dict, List, Tuple

LOG = logging.getLogger("summarize_mag_catalog")

SAMPLE_PREFIX_RE = re.compile(r"^(?P<sample>[^.]+)\.")


def _load_clusters(cdb_csv: Path) -> Tuple[Dict[str, List[str]], Dict[str, str]]:
    """Returns:
        cluster -> [member_genome, …]   (preserves Cdb row order)
        genome -> cluster
    """
    members: Dict[str, List[str]] = defaultdict(list)
    g2c: Dict[str, str] = {}
    with cdb_csv.open() as fh:
        reader = csv.DictReader(fh)
        for row in reader:
            cluster = row["secondary_cluster"]
            genome = row["genome"]
            members[cluster].append(genome)
            g2c[genome] = cluster
    LOG.info("loaded %d clusters spanning %d genomes from %s",
             len(members), len(g2c), cdb_csv.name)
    return members, g2c


def _load_winners(wdb_csv: Path) -> Dict[str, str]:
    out: Dict[str, str] = {}
    with wdb_csv.open() as fh:
        reader = csv.DictReader(fh)
        for row in reader:
            out[row["cluster"]] = row["genome"]
    LOG.info("loaded %d cluster→winner mappings", len(out))
    return out


def _load_genome_info(csv_path: Path) -> Dict[str, Tuple[str, str]]:
    """genome -> (completeness, contamination)"""
    out: Dict[str, Tuple[str, str]] = {}
    with csv_path.open() as fh:
        reader = csv.DictReader(fh)
        for row in reader:
            out[row["genome"]] = (
                row.get("completeness", ""),
                row.get("contamination", ""),
            )
    LOG.info("loaded %d genome QC entries", len(out))
    return out


def _load_gtdb(gtdb_dir: Path) -> Dict[str, str]:
    out: Dict[str, str] = {}
    if not gtdb_dir.exists():
        LOG.warning("GTDB-Tk dir %s missing; lineages will be blank", gtdb_dir)
        return out
    for tsv in sorted(gtdb_dir.glob("*summary.tsv")):
        with tsv.open() as fh:
            reader = csv.DictReader(fh, delimiter="\t")
            if "user_genome" not in (reader.fieldnames or []):
                continue
            for row in reader:
                out[row["user_genome"]] = row.get("classification", "")
    LOG.info("loaded %d GTDB-Tk classifications", len(out))
    return out


def _sample_of(genome_filename: str) -> str:
    """Extract sample id from SAMPLE-N.binette_binM.fa naming used upstream."""
    m = SAMPLE_PREFIX_RE.match(genome_filename)
    return m.group("sample") if m else ""


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cdb", required=True, type=Path)
    parser.add_argument("--wdb", required=True, type=Path)
    parser.add_argument("--genome-info", required=True, type=Path,
                        help="The same genomeInfo.csv passed to dRep")
    parser.add_argument("--gtdb-dir", required=True, type=Path)
    parser.add_argument("--out", required=True, type=Path)
    parser.add_argument("--log-level", default="INFO")
    args = parser.parse_args()

    logging.basicConfig(level=args.log_level,
                        format="%(asctime)s %(levelname)s %(message)s")

    members, _ = _load_clusters(args.cdb)
    winners = _load_winners(args.wdb)
    ginfo = _load_genome_info(args.genome_info)
    gtdb = _load_gtdb(args.gtdb_dir)

    fieldnames = [
        "cluster_id", "winning_genome", "n_members", "n_member_samples",
        "member_genomes", "completeness", "contamination", "gtdb_lineage",
    ]

    n_clusters = 0
    n_with_lineage = 0
    with args.out.open("w") as out_fh:
        writer = csv.DictWriter(out_fh, fieldnames=fieldnames)
        writer.writeheader()
        for cluster in sorted(members.keys()):
            n_clusters += 1
            mems = members[cluster]
            winner = winners.get(cluster, "")
            winner_key = winner[:-3] if winner.endswith(".fa") else winner
            comp, cont = ginfo.get(winner, ("", ""))
            lineage = gtdb.get(winner_key, "")
            if lineage:
                n_with_lineage += 1
            samples = sorted({_sample_of(m) for m in mems if _sample_of(m)})
            writer.writerow({
                "cluster_id": cluster,
                "winning_genome": winner,
                "n_members": len(mems),
                "n_member_samples": len(samples),
                "member_genomes": ";".join(mems),
                "completeness": comp,
                "contamination": cont,
                "gtdb_lineage": lineage,
            })

    LOG.info("wrote %d clusters; %d have a GTDB-Tk lineage", n_clusters, n_with_lineage)
    return 0


if __name__ == "__main__":
    sys.exit(main())
