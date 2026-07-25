#!/usr/bin/env python
"""Join PHOLD per-CDS predictions onto the viral catalog via the cluster table.

Reads:
    --phold-tsv      phold_per_cds_predictions.tsv from PHOLD_COMPARE.
                     Keyed by `cds_id` column (the protein/gene ID we fed
                     to PHOLD).
    --join-table     viral_join_table.tsv from FILTER_STRUCTURES_INPUTS.
                     Columns: viral_rep_id, phold_rep_id, relation.

Emits:
    --out-tsv        virus.proteins.phold.tsv. One row per viral cluster
                     rep, with all PHOLD columns copied from the row keyed
                     by phold_rep_id when relation ∈ {direct, propagated}
                     (the `relation` column is preserved verbatim). When
                     relation == 'none' the PHOLD columns are blank — the
                     row is still emitted so downstream analyses can
                     enumerate all viral reps and see which lack annotation.

Logic:
    direct       phold_rep_id == viral_rep_id; copy that row from PHOLD output.
    propagated   phold_rep_id != viral_rep_id; copy the row keyed by
                 phold_rep_id (the protein-cluster rep) but the row's
                 own cds_id is rewritten to viral_rep_id so downstream
                 consumers can join on the viral catalog directly.
    none         no PHOLD row exists for this viral rep at all (no
                 annotation possible); emit a row with viral_rep_id and
                 empty PHOLD columns.
"""

from __future__ import annotations

import argparse
import csv
import logging
import sys
from pathlib import Path
from typing import Dict, Iterable, List, Optional

LOG = logging.getLogger("merge_viral_phold")

# PHOLD's per-CDS output keys the CDS by the first column. Across PHOLD
# versions the first-column header has been: 'contig_id', 'cds_id',
# 'protein_id', 'gene_name'. We detect dynamically and rename to a stable
# label.
PHOLD_CDS_KEYS = ("cds_id", "protein_id", "contig_id", "gene_name", "ID")


def _detect_cds_column(header: List[str]) -> str:
    """Return the column name in `header` that holds the CDS/protein ID."""
    for c in PHOLD_CDS_KEYS:
        if c in header:
            return c
    raise ValueError(
        f"PHOLD output header has no recognised CDS key (looked for {PHOLD_CDS_KEYS}). "
        f"Got header: {header!r}"
    )


def _load_phold_table(path: Path) -> tuple[List[str], Dict[str, Dict[str, str]]]:
    """Returns (header, {cds_id: row_as_dict})."""
    with path.open(newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        header = reader.fieldnames or []
        if not header:
            raise ValueError(f"empty PHOLD output: {path}")
        cds_col = _detect_cds_column(header)
        out: Dict[str, Dict[str, str]] = {}
        for row in reader:
            cid = row.get(cds_col, "")
            if cid:
                out[cid] = row
    LOG.info("loaded %d PHOLD rows from %s (cds column: %s)", len(out), path, cds_col)
    return header, out


def _load_join_table(path: Path) -> List[tuple[str, str, str]]:
    """Returns list of (viral_rep_id, phold_rep_id, relation) tuples."""
    rows: List[tuple[str, str, str]] = []
    with path.open(newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for r in reader:
            rows.append((
                r.get("viral_rep_id", ""),
                r.get("phold_rep_id", ""),
                r.get("relation", ""),
            ))
    LOG.info("loaded %d viral join rows from %s", len(rows), path)
    return rows


def _emit(
    phold_header: List[str],
    phold_rows: Dict[str, Dict[str, str]],
    join_rows: List[tuple[str, str, str]],
    out_path: Path,
) -> None:
    cds_col = _detect_cds_column(phold_header)
    # Output schema: replace PHOLD's CDS-key column with `viral_rep_id`,
    # keep all other PHOLD columns, prepend `relation` (so analysts see
    # at-a-glance whether a row is direct/propagated/none).
    annotation_cols = [c for c in phold_header if c != cds_col]
    out_header = ["viral_rep_id", "relation", *annotation_cols]

    n_direct = n_propagated = n_none = 0
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with out_path.open("w", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t", lineterminator="\n")
        writer.writerow(out_header)
        for viral_id, phold_id, relation in join_rows:
            if relation == "direct":
                src = phold_rows.get(phold_id)
                if src is None:
                    # Filter promised this row but PHOLD didn't emit it
                    # (e.g. the protein failed prediction). Surface as 'none'
                    # so the consumer can see the gap.
                    writer.writerow([viral_id, "none", *[""] * len(annotation_cols)])
                    n_none += 1
                else:
                    writer.writerow([
                        viral_id, "direct",
                        *(src.get(c, "") for c in annotation_cols),
                    ])
                    n_direct += 1
            elif relation == "propagated":
                src = phold_rows.get(phold_id)
                if src is None:
                    writer.writerow([viral_id, "none", *[""] * len(annotation_cols)])
                    n_none += 1
                else:
                    writer.writerow([
                        viral_id, "propagated",
                        *(src.get(c, "") for c in annotation_cols),
                    ])
                    n_propagated += 1
            else:
                writer.writerow([viral_id, "none", *[""] * len(annotation_cols)])
                n_none += 1
    LOG.info("wrote %d rows to %s (%d direct, %d propagated, %d none)",
             n_direct + n_propagated + n_none, out_path, n_direct, n_propagated, n_none)


def main(argv: Optional[Iterable[str]] = None) -> int:
    p = argparse.ArgumentParser(description=__doc__.split("\n\n", 1)[0])
    p.add_argument("--phold-tsv", required=True, type=Path,
                   help="phold_per_cds_predictions.tsv from PHOLD_COMPARE")
    p.add_argument("--join-table", required=True, type=Path,
                   help="viral_join_table.tsv from FILTER_STRUCTURES_INPUTS")
    p.add_argument("--out-tsv", required=True, type=Path,
                   help="Output viral PHOLD annotations TSV")
    p.add_argument("--log-level", default="INFO",
                   choices=["DEBUG", "INFO", "WARNING", "ERROR"])
    args = p.parse_args(list(argv) if argv is not None else None)
    logging.basicConfig(level=args.log_level, format="%(asctime)s | %(levelname)s | %(message)s")

    header, phold_rows = _load_phold_table(args.phold_tsv)
    join_rows = _load_join_table(args.join_table)
    _emit(header, phold_rows, join_rows, args.out_tsv)
    return 0


if __name__ == "__main__":
    sys.exit(main())
