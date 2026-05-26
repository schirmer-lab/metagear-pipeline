#!/usr/bin/env python
"""Merge per-contig classification evidence into a single TSV per sample.

Consumes evidence from the integrated_classification chain:

    * Assembly contigs FASTA — defines the universe of contig_ids + lengths
    * Viral IDs (one ID per line, optional ``seq_name`` header) — from MERGE_TABLES
    * Plasmid IDs — from MERGE_TABLES
    * Binette final bins directory — one ``.fa`` per MAG, headers = bin members
    * MMseqs2 easy-taxonomy LCA TSV (optional) — lineage for unbinned contigs
    * Tiara classifications (.txt or .txt.gz) — per-contig domain calls

Priority for ``primary_class`` (highest to lowest, first match wins):

    1. ``virus``     — contig in viral_ids
    2. ``plasmid``   — contig in plasmid_ids
    3. ``bacteria``  — contig is a member of a Binette bin (assume bacteria;
                       v2 cohort_dereplication + GTDB-Tk will refine into
                       bacteria vs archaea)
    4. ``eukaryote`` — Tiara called it eukarya
    5. ``archaea``   — Tiara called it archaea AND no bin/virus/plasmid above
    6. lineage from MMseqs2 LCA's superkingdom — when available, set
                       primary_class to bacteria / archaea / eukaryote / virus
                       based on the parsed superkingdom prefix
    7. ``unknown``   — none of the above hit

The ``classifier`` column records which step assigned ``primary_class``;
``lineage`` is populated when available (mmseqs2 LCA's full lineage string).
"""

from __future__ import annotations

import argparse
import csv
import gzip
import logging
import os
import sys
from collections import defaultdict
from pathlib import Path
from typing import IO, Dict, Iterable, Iterator, Set, Tuple

LOG = logging.getLogger("merge_contig_classification")
TSV_HEADER = [
    "contig_id",
    "sample",
    "length",
    "primary_class",
    "classifier",
    "lineage",
    "confidence",
    "bin_id",
]


def _open_text(path: Path) -> IO[str]:
    """Open a (possibly gzipped) text file for reading."""
    if str(path).endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path, "rt")


def _iter_fasta_records(path: Path) -> Iterator[Tuple[str, int]]:
    """Yield (contig_id, length) for each record in a FASTA (gzipped or not).

    Length is the number of non-whitespace residues in the sequence. Contig_id
    is the first token of the header line (whitespace-split, ``>`` stripped) —
    matching the convention every other tool in this pipeline uses.
    """
    contig_id: str | None = None
    length = 0
    with _open_text(path) as fh:
        for line in fh:
            if not line:
                continue
            if line.startswith(">"):
                if contig_id is not None:
                    yield contig_id, length
                # First whitespace-delimited token after '>'
                contig_id = line[1:].split()[0] if len(line) > 1 else ""
                length = 0
            else:
                length += len(line.strip())
        if contig_id is not None:
            yield contig_id, length


def _load_id_list(path: Path | None) -> Set[str]:
    """Load a list of contig IDs (one per line). Tolerates the seq_name header.

    Returns the empty set when ``path`` is None, missing, or zero-byte.
    """
    if path is None or not path.exists() or path.stat().st_size == 0:
        return set()
    ids: Set[str] = set()
    with _open_text(path) as fh:
        for raw in fh:
            line = raw.strip()
            if not line or line.startswith("#"):
                continue
            if line == "seq_name":  # header line from MERGE_TABLES
                continue
            # Strip trailing per-row metadata; we only care about the ID.
            ids.add(line.split()[0])
    return ids


def _load_bin_membership(bin_dir: Path | None) -> Dict[str, str]:
    """Return ``{contig_id: bin_id}`` for every member contig of every bin.

    Bin files are ``${bin_dir}/*.fa`` (Binette default). ``bin_id`` is the
    filename stem (e.g. ``bin_3.fa`` → ``bin_3``).
    """
    membership: Dict[str, str] = {}
    if bin_dir is None or not bin_dir.exists():
        return membership
    if not bin_dir.is_dir():
        LOG.warning("bin_dir %s is not a directory; skipping", bin_dir)
        return membership
    for bin_file in sorted(bin_dir.iterdir()):
        if not bin_file.is_file():
            continue
        if bin_file.suffix not in (".fa", ".fasta") and not bin_file.name.endswith((".fa.gz", ".fasta.gz")):
            continue
        bin_id = bin_file.name
        for suffix in (".fa.gz", ".fasta.gz", ".fa", ".fasta"):
            if bin_id.endswith(suffix):
                bin_id = bin_id[: -len(suffix)]
                break
        with _open_text(bin_file) as fh:
            for line in fh:
                if line.startswith(">"):
                    contig_id = line[1:].split()[0] if len(line) > 1 else ""
                    if contig_id:
                        membership[contig_id] = bin_id
    return membership


def _load_mmseqs_lca(path: Path | None) -> Dict[str, Tuple[str, str]]:
    """Parse the ``${prefix}_lca.tsv`` output of ``mmseqs easy-taxonomy``.

    Columns (verified against mmseqs2 17.b804f output):
        1. query contig id
        2. taxonomy id
        3. rank
        4. taxonomy name
        5. number of fragments retained
        6. number of fragments assigned to the LCA
        7. agreement / consensus count
        8. confidence (fraction at the LCA, 0..1)
        9. lineage (semicolon-separated; e.g. d_Bacteria;p_Bacteroidota;…)
           — present when easy-taxonomy was invoked with ``--tax-lineage 1``,
           which the pipeline always sets.

    Earlier versions of mmseqs2 emitted only 8 columns (no agreement count);
    we tolerate both by indexing from the end for the lineage and confidence,
    so the parser keeps working if mmseqs2 changes the upstream columns again.

    Returns ``{contig_id: (lineage_str, confidence_str)}`` keyed by contig id.
    Empty when ``path`` is None or zero-byte.
    """
    out: Dict[str, Tuple[str, str]] = {}
    if path is None or not path.exists() or path.stat().st_size == 0:
        return out
    with _open_text(path) as fh:
        for raw in fh:
            line = raw.rstrip("\n")
            if not line or line.startswith("#"):
                continue
            cols = line.split("\t")
            if len(cols) < 8:
                continue
            contig_id = cols[0]
            # Lineage is always the last column when --tax-lineage 1 is set;
            # confidence is the column immediately before it. This is robust to
            # the 8-col (older mmseqs) vs 9-col (current) layout.
            lineage = cols[-1]
            confidence = cols[-2]
            out[contig_id] = (lineage, confidence)
    return out


def _tiara_kingdom_to_class(label: str) -> str | None:
    """Map a Tiara label to our primary_class enum, or None if unmappable.

    Tiara labels (per tiara docs): archaea, bacteria, eukarya, mitochondrion,
    plast, unknown, organelle. We collapse mitochondrion/plast/organelle to
    eukaryote (they're eukaryotic organelles) but plasmid/virus calls win
    upstream of this anyway.
    """
    label = label.strip().lower()
    if label in ("archaea",):
        return "archaea"
    if label in ("bacteria",):
        return "bacteria"
    if label in ("eukarya", "mitochondrion", "plast", "organelle"):
        return "eukaryote"
    return None


def _load_tiara(path: Path | None) -> Dict[str, str]:
    """Parse Tiara's classifications .txt.

    File layout (Tiara default):

        sequence_id    class_fst_stage    class_snd_stage

    The 2nd stage is more refined (e.g. resolves "eukarya" → "mitochondrion")
    but is only set when first stage was eukaryotic. We prefer 2nd-stage when
    populated, else 1st-stage.
    """
    out: Dict[str, str] = {}
    if path is None or not path.exists() or path.stat().st_size == 0:
        return out
    with _open_text(path) as fh:
        header_skipped = False
        for raw in fh:
            line = raw.rstrip("\n")
            if not line or line.startswith("#"):
                continue
            if not header_skipped:
                header_skipped = True
                # Tiara always emits a header line
                if "sequence_id" in line.lower() or "class_fst" in line.lower():
                    continue
            cols = line.split("\t")
            if len(cols) < 2:
                continue
            contig_id = cols[0]
            label = cols[2].strip() if len(cols) >= 3 and cols[2].strip() and cols[2].strip() != "n/a" else cols[1]
            mapped = _tiara_kingdom_to_class(label)
            if mapped:
                out[contig_id] = mapped
    return out


def _lineage_to_class(lineage: str) -> str | None:
    """Extract primary_class from a GTDB-style lineage prefix.

    Examples:
        ``d__Bacteria;p__Bacteroidota;...`` → bacteria
        ``d__Archaea;p__...``               → archaea
        ``d__Eukaryota;p__...``             → eukaryote
        ``d__Viruses;...``                  → virus

    Falls back to None if no recognizable domain prefix is found.
    """
    if not lineage:
        return None
    head = lineage.split(";", 1)[0].lower()
    if "bacteria" in head:
        return "bacteria"
    if "archaea" in head:
        return "archaea"
    if "eukaryota" in head or "eukarya" in head:
        return "eukaryote"
    if "virus" in head or "viruses" in head:
        return "virus"
    return None


def merge(
    contigs_fa: Path,
    sample: str,
    viral_ids_path: Path | None,
    plasmid_ids_path: Path | None,
    bin_dir: Path | None,
    mmseqs_lca_path: Path | None,
    tiara_path: Path | None,
    out_tsv: Path,
) -> None:
    LOG.info("Loading evidence channels…")
    viral_ids = _load_id_list(viral_ids_path)
    plasmid_ids = _load_id_list(plasmid_ids_path)
    bin_membership = _load_bin_membership(bin_dir)
    mmseqs_lca = _load_mmseqs_lca(mmseqs_lca_path)
    tiara_labels = _load_tiara(tiara_path)

    LOG.info(
        "Loaded %d viral, %d plasmid, %d binned, %d mmseqs LCA, %d tiara labels",
        len(viral_ids),
        len(plasmid_ids),
        len(bin_membership),
        len(mmseqs_lca),
        len(tiara_labels),
    )

    written = 0
    class_counts: Dict[str, int] = defaultdict(int)
    classifier_counts: Dict[str, int] = defaultdict(int)

    with open(out_tsv, "wt", newline="") as out_fh:
        writer = csv.writer(out_fh, delimiter="\t", lineterminator="\n")
        writer.writerow(TSV_HEADER)

        for contig_id, length in _iter_fasta_records(contigs_fa):
            primary_class = "unknown"
            classifier = "unknown"
            lineage = ""
            confidence = ""
            bin_id = ""

            if contig_id in viral_ids:
                primary_class = "virus"
                classifier = "genomad"
            elif contig_id in plasmid_ids:
                primary_class = "plasmid"
                classifier = "genomad"
            elif contig_id in bin_membership:
                # v1: every Binette bin is assumed bacterial (we filter at MIMAG
                # MQ+ which is a bacterial-MAG threshold). v2 cohort_dereplication
                # + GTDB-Tk will refine bacteria vs archaea at the species level.
                primary_class = "bacteria"
                classifier = "binette"
                bin_id = bin_membership[contig_id]
            elif contig_id in mmseqs_lca:
                lineage, confidence = mmseqs_lca[contig_id]
                class_from_lineage = _lineage_to_class(lineage)
                if class_from_lineage:
                    primary_class = class_from_lineage
                    classifier = "mmseqs2"
            if primary_class == "unknown" and contig_id in tiara_labels:
                primary_class = tiara_labels[contig_id]
                classifier = "tiara"

            writer.writerow(
                [contig_id, sample, length, primary_class, classifier, lineage, confidence, bin_id]
            )
            written += 1
            class_counts[primary_class] += 1
            classifier_counts[classifier] += 1

    LOG.info("Wrote %d rows to %s", written, out_tsv)
    LOG.info(
        "primary_class counts: %s",
        ", ".join(f"{k}={v}" for k, v in sorted(class_counts.items())),
    )
    LOG.info(
        "classifier counts:    %s",
        ", ".join(f"{k}={v}" for k, v in sorted(classifier_counts.items())),
    )


def _parse_args(argv: Iterable[str] | None = None) -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__.split("\n", 1)[0])
    p.add_argument("--contigs", required=True, type=Path, help="Per-sample assembly contigs FASTA (gzipped ok).")
    p.add_argument("--sample", required=True, type=str, help="Sample id (filled into the `sample` column).")
    p.add_argument("--viral-ids", type=Path, default=None, help="Viral contig IDs (one per line; optional `seq_name` header tolerated).")
    p.add_argument("--plasmid-ids", type=Path, default=None, help="Plasmid contig IDs.")
    p.add_argument("--bin-dir", type=Path, default=None, help="Directory of Binette final_bins/*.fa")
    p.add_argument("--mmseqs-lca", type=Path, default=None, help="mmseqs2 easy-taxonomy <prefix>_lca.tsv")
    p.add_argument("--tiara", type=Path, default=None, help="Tiara classifications .txt (gzipped ok).")
    p.add_argument("--output", required=True, type=Path, help="Output per-contig classification TSV.")
    p.add_argument("--log-level", default="INFO", choices=["DEBUG", "INFO", "WARNING", "ERROR"])
    return p.parse_args(argv)


def main(argv: Iterable[str] | None = None) -> int:
    args = _parse_args(argv)
    logging.basicConfig(level=args.log_level, format="%(asctime)s | %(levelname)s | %(message)s")
    merge(
        contigs_fa=args.contigs,
        sample=args.sample,
        viral_ids_path=args.viral_ids,
        plasmid_ids_path=args.plasmid_ids,
        bin_dir=args.bin_dir,
        mmseqs_lca_path=args.mmseqs_lca,
        tiara_path=args.tiara,
        out_tsv=args.output,
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
