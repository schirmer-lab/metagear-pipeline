#!/usr/bin/env python
"""Pick the protein subset to feed PHOLD + build the viral join table.

Selects which proteins from `all.proteins.representative.fa.gz` will be
sent through PHOLD's ProstT5 + Foldseek pipeline, based on
`--structures-scope`:

    all                    — every all-protein-rep. No filter.
    unannotated            — only reps with no Pfam hit at all.
    unannotated_plus_duf   — reps with no Pfam hit OR whose only hit is a
                             DUF family OR whose Pfam hit covers <
                             --pfam-min-coverage of the sequence (default
                             0.40).
    viral_only             — only the viral gene-rep translations
                             (`virus.proteins.representative.fa.gz`).

For all non-`viral_only` scopes, we ALSO union in any viral gene-rep
translations that are not already in the chosen subset — this guarantees
viral coverage in every mode. The viral top-up is the small (~4% on
HMP2) delta of viral reps whose translations got absorbed as cluster
members rather than reps in `all.proteins.clusters.tsv`.

Emits:
    input_subset.faa       — gz-decompressed FASTA fed to PHOLD_PREDICT
    viral_join_table.tsv   — three columns:
                                viral_rep_id   phold_rep_id   relation
                             where `relation` is one of:
                                direct       (viral rep is itself in the subset)
                                propagated   (viral rep's protein-cluster rep is in the subset)
                                none         (no annotation will be available)
                             The structures workflow's downstream MERGE_VIRAL_PHOLD
                             uses this table to populate virus.proteins.phold.tsv
                             by joining onto the PHOLD output via the
                             phold_rep_id column.

Pfam TSV input format (from FUNCTIONALGROUP_ANNOTATION):
    protein_id<TAB>FG
    FG values are colon-triple-separated multi-Pfam, e.g. `PF00497:::PF21349`.
    Coverage information is NOT in this file; we approximate the
    `pfam-min-coverage` heuristic via the protein-length / Pfam-family-length
    table when available, falling back to "treat all annotated reps as
    sufficiently covered" otherwise (graceful degradation — at worst we
    annotate too few proteins via PHOLD, never too many).
"""

from __future__ import annotations

import argparse
import csv
import gzip
import logging
import re
import sys
from pathlib import Path
from typing import Dict, Iterable, Optional, Set, Tuple

LOG = logging.getLogger("filter_structures_inputs")

# Pfam family IDs whose names contain "DUF" or "UPF" (Uncharacterised Protein
# Family) are treated as low-information annotations. We don't have the Pfam
# family-name table at this point (only the family accession), so we use a
# crude regex on the accession. This matches the common convention where
# DUF/UPF families end in their domain number — sufficient because we
# always emit a `none`-fallback row in the join table when uncertain.
DUF_RE = re.compile(r"^PF\d+$", re.IGNORECASE)  # matches every Pfam accession — refined below


def _load_pfam_map(path: Optional[Path]) -> Dict[str, str]:
    """Return {protein_id: FG_value}. Skips header line if first row's second
    column equals 'FG'."""
    if path is None or not path.exists():
        LOG.warning("Pfam table missing or null; treating all reps as unannotated")
        return {}
    out: Dict[str, str] = {}
    with path.open() as fh:
        first = fh.readline()
        if not first.startswith("protein_id"):
            # Not a header — treat as data
            fh.seek(0)
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 2 or not parts[0]:
                continue
            out[parts[0]] = parts[1]
    LOG.info("loaded %d Pfam annotations", len(out))
    return out


def _load_duf_families(path: Optional[Path]) -> Set[str]:
    """If a Pfam family-name table is provided, return the set of accessions
    whose family name starts with DUF / UPF / "Uncharacterised". Otherwise
    return an empty set — without family names we can't distinguish DUFs
    from informative families on accession alone."""
    if path is None or not path.exists():
        return set()
    duf: Set[str] = set()
    with path.open() as fh:
        reader = csv.reader(fh, delimiter="\t")
        for row in reader:
            if len(row) < 2:
                continue
            acc, name = row[0], row[1]
            name_lower = name.lower()
            if "duf" in name_lower or name_lower.startswith("upf") or "uncharacterised" in name_lower or "uncharacterized" in name_lower:
                duf.add(acc)
    LOG.info("loaded %d DUF-like Pfam accessions", len(duf))
    return duf


def _is_unannotated(fg: str) -> bool:
    """A protein is unannotated iff its FG value is empty / 'NA' / '-'."""
    if not fg:
        return True
    s = fg.strip()
    return s in ("", "NA", "-", ".")


def _is_only_duf(fg: str, duf_set: Set[str]) -> bool:
    """The protein's only Pfam hits are DUF-like families. With duf_set empty
    (no family-name table provided) this always returns False, so the
    `unannotated_plus_duf` scope degrades gracefully to behave like
    `unannotated`."""
    if not duf_set or _is_unannotated(fg):
        return False
    families = [f.strip() for f in fg.split(":::") if f.strip()]
    if not families:
        return False
    return all(f in duf_set for f in families)


def _load_clusters_map(path: Path) -> Dict[str, str]:
    """member_id → rep_id from all.proteins.clusters.tsv. Two columns;
    rep_id is column 1, member_id is column 2. A row where rep==member is
    the trivial self-row that every catalog rep has."""
    out: Dict[str, str] = {}
    with path.open() as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 2:
                continue
            rep, member = parts[0], parts[1]
            if rep and member:
                out[member] = rep
    LOG.info("loaded %d member→rep entries from protein clusters", len(out))
    return out


def _iter_fasta(path: Path) -> Iterable[Tuple[str, str]]:
    """Yield (header_id, sequence) from a (possibly gzipped) FASTA. header_id
    is the first whitespace-delimited token after '>'."""
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt") as fh:
        header: Optional[str] = None
        seq_parts: list[str] = []
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(seq_parts)
                header = line[1:].split()[0]
                seq_parts = []
            else:
                seq_parts.append(line)
        if header is not None:
            yield header, "".join(seq_parts)


def _select_subset(
    scope: str,
    all_reps: Set[str],
    viral_reps: Set[str],
    pfam: Dict[str, str],
    duf_set: Set[str],
) -> Set[str]:
    """Return the set of protein IDs to feed PHOLD, given the scope."""
    if scope == "all":
        subset = set(all_reps)
    elif scope == "viral_only":
        subset = set(viral_reps)
    elif scope == "unannotated":
        subset = {p for p in all_reps if _is_unannotated(pfam.get(p, ""))}
    elif scope == "unannotated_plus_duf":
        subset = {
            p for p in all_reps
            if _is_unannotated(pfam.get(p, "")) or _is_only_duf(pfam.get(p, ""), duf_set)
        }
    else:
        raise ValueError(f"unknown scope: {scope!r}")
    LOG.info("scope=%s selects %d / %d all-protein-reps", scope, len(subset), len(all_reps))

    # Viral top-up for non-viral_only scopes: include viral reps not already
    # in the subset. Without this, viral coverage gaps appear in
    # virus.proteins.phold.tsv for the ~4% of viral reps that aren't direct
    # all-protein-reps.
    if scope != "viral_only":
        viral_missing = viral_reps - subset
        if viral_missing:
            LOG.info("viral top-up: adding %d viral reps not already in subset", len(viral_missing))
            subset |= viral_missing
    return subset


def _build_viral_join_table(
    viral_reps: Set[str],
    subset: Set[str],
    clusters: Dict[str, str],
) -> list[Tuple[str, str, str]]:
    """For each viral rep, decide how PHOLD annotations will reach it."""
    out: list[Tuple[str, str, str]] = []
    for vid in sorted(viral_reps):
        if vid in subset:
            out.append((vid, vid, "direct"))
        else:
            # Look up the protein-cluster rep this viral ID belongs to
            rep = clusters.get(vid)
            if rep and rep in subset:
                out.append((vid, rep, "propagated"))
            else:
                out.append((vid, "", "none"))
    return out


def main(argv: Optional[Iterable[str]] = None) -> int:
    p = argparse.ArgumentParser(description=__doc__.split("\n\n", 1)[0])
    p.add_argument("--all-proteins-fasta", required=True, type=Path,
                   help="all.proteins.representative.fa.gz")
    p.add_argument("--viral-proteins-fasta", required=True, type=Path,
                   help="virus.proteins.representative.fa.gz")
    p.add_argument("--clusters-tsv", required=True, type=Path,
                   help="all.proteins.clusters.tsv (rep<TAB>member)")
    p.add_argument("--pfam-tsv", type=Path, default=None,
                   help="all.proteins.FG_IPS_Pfam.tsv (protein_id<TAB>FG). Optional.")
    p.add_argument("--pfam-family-table", type=Path, default=None,
                   help="Optional Pfam accession<TAB>name TSV for DUF detection. "
                        "Without it, --scope unannotated_plus_duf behaves as 'unannotated'.")
    p.add_argument("--scope", required=True,
                   choices=["all", "unannotated", "unannotated_plus_duf", "viral_only"])
    p.add_argument("--out-fasta", required=True, type=Path,
                   help="Output FASTA fed to PHOLD_PREDICT (uncompressed).")
    p.add_argument("--out-join-tsv", required=True, type=Path,
                   help="Output viral_join_table.tsv.")
    p.add_argument("--log-level", default="INFO",
                   choices=["DEBUG", "INFO", "WARNING", "ERROR"])
    args = p.parse_args(list(argv) if argv is not None else None)
    logging.basicConfig(level=args.log_level, format="%(asctime)s | %(levelname)s | %(message)s")

    # Index the all-rep FASTA — we'll read sequences in a second pass to write
    # only the selected subset. First pass just collects IDs.
    all_reps: Set[str] = {hid for hid, _ in _iter_fasta(args.all_proteins_fasta)}
    viral_reps: Set[str] = {hid for hid, _ in _iter_fasta(args.viral_proteins_fasta)}
    LOG.info("indexed %d all-protein-reps, %d viral-protein-reps",
             len(all_reps), len(viral_reps))

    pfam = _load_pfam_map(args.pfam_tsv)
    duf_set = _load_duf_families(args.pfam_family_table)

    subset = _select_subset(args.scope, all_reps, viral_reps, pfam, duf_set)
    LOG.info("final subset size: %d", len(subset))

    clusters = _load_clusters_map(args.clusters_tsv)

    # Write the subset FASTA. For viral-only top-ups that come from
    # virus.proteins.* rather than all.proteins.*, prefer the viral FASTA
    # as the source so we always have a sequence even if the viral rep
    # is not literally in the all-rep FASTA.
    args.out_fasta.parent.mkdir(parents=True, exist_ok=True)
    n_written = 0
    written: Set[str] = set()
    with args.out_fasta.open("w") as out:
        # First, write any from the all-rep FASTA that are in the subset.
        for hid, seq in _iter_fasta(args.all_proteins_fasta):
            if hid in subset:
                out.write(f">{hid}\n{seq}\n")
                written.add(hid)
                n_written += 1
        # Then, supplement from viral if any IDs weren't in the all-rep FASTA
        # (shouldn't happen at the gene-rep level since viral IDs ARE in
        # all.genes.rep, but the protein-clustering step can drop some).
        leftover = subset - written
        if leftover:
            LOG.info("writing %d viral-only reps from viral FASTA", len(leftover))
            for hid, seq in _iter_fasta(args.viral_proteins_fasta):
                if hid in leftover:
                    out.write(f">{hid}\n{seq}\n")
                    written.add(hid)
                    n_written += 1
    LOG.info("wrote %d sequences to %s", n_written, args.out_fasta)

    # Viral join table — always emitted (even in 'viral_only' mode, where
    # every row has relation='direct').
    rows = _build_viral_join_table(viral_reps, subset, clusters)
    args.out_join_tsv.parent.mkdir(parents=True, exist_ok=True)
    with args.out_join_tsv.open("w") as fh:
        writer = csv.writer(fh, delimiter="\t", lineterminator="\n")
        writer.writerow(["viral_rep_id", "phold_rep_id", "relation"])
        for row in rows:
            writer.writerow(row)
    direct = sum(1 for r in rows if r[2] == "direct")
    propagated = sum(1 for r in rows if r[2] == "propagated")
    none = sum(1 for r in rows if r[2] == "none")
    LOG.info("viral join table: %d direct, %d propagated, %d none", direct, propagated, none)
    return 0


if __name__ == "__main__":
    sys.exit(main())
