#!/usr/bin/env python3
"""
VirSorter2-style affi builder for DRAM-v using ONLY:
  - Pharokka per-protein TSV (ID, phrog, annot)
  - Protein FASTA with headers that encode coordinates: {contig}::{idx}::{start}::{end}::{strand}
Optional:
  - Contigs FASTA (to set exact contig lengths)
  - Circular contig list

Why this version?
- No Prodigal needed. We enumerate ORFs from your existing protein FASTA headers.
- Joins Pharokka by (contig,start,end,strand) with ±5 nt tolerance (configurable).
- Emits VirSorter2-like gene IDs (contig__1..N) and 12-column, pipe-delimited lines.
- Guarantees every contig in the contigs FASTA appears in the affi (N=0 allowed).

Category codes:
  0 = viral hallmark (structural/packaging), 1 = viral (non-hallmark), 2 = unknown

Usage:
  python build_affi_from_pharokka.py \
    --pharokka pharokka.tsv \
    --proteins-faa proteins_all.faa \
    --contigs-fasta viral_contigs.fa \
    --out-affi viral-affi-contigs-for-dramv.tab
"""

import argparse, csv, gzip, io, re, sys
from collections import defaultdict, namedtuple
from typing import Dict, List, Optional, Tuple

# ----------------------------
# IO helpers
# ----------------------------


def open_maybe_gzip(path: Optional[str]):
    if path is None:
        return None
    if path.endswith(".gz"):
        return io.TextIOWrapper(gzip.open(path, "rb"), encoding="utf-8", newline="")
    return open(path, "r", encoding="utf-8", newline="")


def fasta_lengths(path: Optional[str]) -> Dict[str, int]:
    if not path:
        return {}
    lens: Dict[str, int] = {}
    with open_maybe_gzip(path) as f:
        name, L = None, 0
        for line in f:
            if line.startswith(">"):
                if name is not None:
                    lens[name] = L
                name = line[1:].strip().split()[0]
                L = 0
            else:
                L += len(line.strip())
        if name is not None:
            lens[name] = L
    return lens


def read_id_set(path: Optional[str]) -> set:
    s = set()
    if not path:
        return s
    with open_maybe_gzip(path) as f:
        for line in f:
            t = line.strip()
            if t and not t.startswith("#"):
                s.add(t)
    return s


# ----------------------------
# Hallmark classification (conservative, tuned to PHROGs vocab)
# ----------------------------

EXCLUDE_LYSIS = re.compile(r"\b(holin|endolysin|lysozyme|spanin)\b", re.I)
HALLMARK_PATTERN = re.compile(
    r"""(?xi)
    # Packaging & head
    \bportal\b|
    \bterminase\b|
    \bmajor\s+head\b|\bminor\s+head\b|
    \bhead\s+(?:assembly|closure|decoration|fiber|fibers|scaffolding|maturation|morphogenesis|protease|protein)\b|
    \bhead\s+and\s+packaging\b|
    # Connector / neck / collar
    \bhead[-\s]?tail\s+(?:connector|adaptor)\b|\bconnector\b|
    \b(?:pre-?)?neck(?:\s*\d+)?\b|\b(?:upper|lower)\s+collar\b|\bcollar\b|
    # Baseplate
    \bbaseplate(?:\s+(?:assembly|hub|wedge|spike|protein))?\b|
    # Tail morphogenesis
    \bdistal\s+tail\b|\bcentral\s+tail\b|\bdit\b|
    \btail\s+(?:protein|tube|sheath|spike|terminator|completion|assembly|chaperone|
              collar(?:\s+protein|\s+fiber)?|
              fiber(?:\s+(?:assembly|protein|and\s+host\s+specificity))?|
              length\s+tape\s+measure)\b|
    \breceptor\s+binding\s+tail\s+protein\b|\bpre-?neck\s+appendage\b
"""
)


def is_hallmark(annot: str) -> bool:
    if not annot:
        return False
    a = annot.lower()
    if EXCLUDE_LYSIS.search(a):
        # Only treat as hallmark if clearly embedded in tail structure context
        if "tail" in a and any(
            k in a
            for k in (
                "fiber",
                "tube",
                "sheath",
                "spike",
                "collar",
                "tape measure",
                "terminator",
                "assembly",
                "chaperone",
            )
        ):
            return True
        return False
    return bool(HALLMARK_PATTERN.search(a))


# ----------------------------
# Pharokka TSV parsing
# ----------------------------


def detect_pharokka_header(header: List[str]) -> Dict[str, int]:
    lc = [h.strip().lower() for h in header]
    need = {}
    for k in ("id", "phrog", "annot"):
        if k not in lc:
            raise KeyError(
                f"Required column '{k}' not found in Pharokka TSV header: {header}"
            )
        need[k] = lc.index(k)
    if "length" in lc:
        need["length"] = lc.index("length")
    if "category" in lc:
        need["category"] = lc.index("category")
    return need


def parse_pharokka_id(gid: str) -> Tuple[str, int, int, int, str]:
    # {contig}::{gene_index}::{start}::{end}::{strand}
    parts = gid.split("::")
    if len(parts) < 5:
        raise ValueError(f"ID not in expected pattern: {gid}")
    contig = parts[0]
    gi = int(parts[1])
    start, end = int(parts[2]), int(parts[3])
    strand = parts[4] if parts[4] in ("+", "-") else "+"
    return contig, gi, start, end, strand


# ----------------------------
# Protein FASTA parsing (authoritative ORF list)
# ----------------------------

ORF = namedtuple("ORF", "contig gi start end strand")


def parse_protein_header(h: str) -> Optional[ORF]:
    """
    Expect: >contig::idx::start::end::strand
    Returns ORF(contig, idx, start, end, strand) or None if unparsable.
    """
    h = h.strip()
    if h.startswith(">"):
        h = h[1:]
    # Split off trailing description if present
    h = h.split()[0]
    parts = h.split("::")
    if len(parts) < 5:
        return None
    try:
        contig = parts[0]
        gi = int(parts[1])
        start = int(parts[2])
        end = int(parts[3])
        strand = parts[4] if parts[4] in ("+", "-") else "+"
        return ORF(contig, gi, start, end, strand)
    except Exception:
        return None


def load_orfs_from_proteins(faa_path: str) -> Dict[str, List[ORF]]:
    """
    Build contig -> list of ORFs from protein FASTA headers.
    For duplicated coords/indices, keep the first seen.
    """
    per_contig: Dict[str, List[ORF]] = defaultdict(list)
    seen = set()
    with open_maybe_gzip(faa_path) as f:
        for line in f:
            if line.startswith(">"):
                orf = parse_protein_header(line.strip())
                if orf is None:
                    continue
                key = (orf.contig, orf.gi, orf.start, orf.end, orf.strand)
                if key in seen:
                    continue
                seen.add(key)
                per_contig[orf.contig].append(orf)
    # Sort by provided index, then by start for stability; enforce 1..N numbering later
    for c in per_contig:
        per_contig[c].sort(key=lambda r: (r.gi, r.start, r.end))
    return per_contig


# ----------------------------
# Core builder
# ----------------------------


def build_affi(
    pharokka_tsv: str,
    proteins_faa: str,
    out_affi: str,
    contigs_fasta: Optional[str],
    circular_list: Optional[str],
    labels_out: Optional[str],
    coord_tolerance: int,
) -> None:

    contig_len = fasta_lengths(contigs_fasta)
    circular_ids = read_id_set(circular_list)
    orfs_by_contig = load_orfs_from_proteins(proteins_faa)

    # Map Pharokka by (contig,start,end,strand)
    pk: Dict[Tuple[str, int, int, str], Tuple[str, str]] = {}
    with open_maybe_gzip(pharokka_tsv) as pf:
        reader = csv.reader(pf, delimiter="\t")
        header = next(reader)
        idx = detect_pharokka_header(header)
        for row in reader:
            if not row or len(row) < len(header):
                continue
            gid = row[idx["id"]].strip()
            phrog = row[idx["phrog"]].strip()
            annot = (row[idx["annot"]] or "").strip()
            try:
                contig, gi, start, end, strand = parse_pharokka_id(gid)
            except Exception:
                continue
            pk[(contig, start, end, strand)] = (phrog, annot)

    # tolerant lookup
    def pharokka_hit(contig: str, s: int, e: int, strand: str):
        hit = pk.get((contig, s, e, strand))
        if hit is not None or coord_tolerance <= 0:
            return hit
        for ds in range(-coord_tolerance, coord_tolerance + 1):
            for de in range(-coord_tolerance, coord_tolerance + 1):
                h = pk.get((contig, s + ds, e + de, strand))
                if h is not None:
                    return h
        return None

    # If contigs FASTA was provided, collect its contig order to also emit headers for contigs with 0 ORFs
    fasta_contigs = []
    if contigs_fasta:
        with open_maybe_gzip(contigs_fasta) as f:
            for line in f:
                if line.startswith(">"):
                    fasta_contigs.append(line[1:].strip().split()[0])

    # Emit
    total_orfs = 0
    matched = 0
    with open(out_affi, "w", encoding="utf-8") as out:
        # First, contigs with ORFs from proteins_faa
        for contig in sorted(orfs_by_contig):
            orfs = orfs_by_contig[contig]
            # Enforce continuous numbering 1..N in VirSorter2 style, preserving order
            orfs_sorted = sorted(orfs, key=lambda r: (r.gi, r.start, r.end))
            N = len(orfs_sorted)
            topo = "c" if contig in circular_ids else "l"
            clen = contig_len.get(contig, max(r.end for r in orfs_sorted))
            out.write(f">{contig}|{N}|{topo}\n")

            for i, r in enumerate(orfs_sorted, start=1):
                total_orfs += 1
                hit = pharokka_hit(contig, r.start, r.end, r.strand)
                if hit:
                    phrog, annot = hit
                    if phrog and phrog != "No_PHROG":
                        top_hit = f"phrog_{phrog}"
                        cat = "0" if is_hallmark(annot) else "1"
                    else:
                        top_hit = "-"
                        cat = "2"
                    matched += 1 if phrog and phrog != "No_PHROG" else 0
                else:
                    top_hit = "-"
                    cat = "2"

                # VirSorter2-style 12 columns
                row = [
                    f"{contig}__{i}",
                    str(r.start),
                    str(r.end),
                    str(clen),
                    r.strand,
                    top_hit,
                    "nan",  # score
                    "-",  # evalue
                    cat,  # category
                    "-",  # pfam_id
                    "nan",  # pfam_score
                    "-",  # trailing
                ]
                out.write("|".join(row) + "\n")

        # Then, contigs present in FASTA but with 0 ORFs seen in proteins_faa -> still emit header
        if fasta_contigs:
            for contig in fasta_contigs:
                if contig in orfs_by_contig:
                    continue
                topo = "c" if contig in circular_ids else "l"
                clen = contig_len.get(contig, 0)
                out.write(f">{contig}|0|{topo}\n")  # no per-ORF lines

    # Summary
    print(
        f"[INFO] Contigs with ORFs (proteins_faa): {len(orfs_by_contig)}",
        file=sys.stderr,
    )
    if fasta_contigs:
        missing = [c for c in fasta_contigs if c not in orfs_by_contig]
        if missing:
            print(
                f"[WARN] {len(missing)} contig(s) in FASTA had no ORFs in proteins_faa. Headers with N=0 were emitted.",
                file=sys.stderr,
            )
    print(f"[INFO] Total ORFs written            : {total_orfs}", file=sys.stderr)
    print(f"[INFO] ORFs with PHROG annotations  : {matched}", file=sys.stderr)
    if matched == 0:
        print(
            "[WARN] No PHROG matches were joined; check coordinate patterns and tolerance.",
            file=sys.stderr,
        )


# ----------------------------
# CLI
# ----------------------------


def main():
    ap = argparse.ArgumentParser(
        description="VirSorter2-style affi for DRAM-v from Pharokka + protein FASTA (no Prodigal needed)."
    )
    ap.add_argument(
        "--pharokka",
        required=True,
        help="Pharokka TSV (columns: ID, phrog, annot). IDs encode coords.",
    )
    ap.add_argument(
        "--proteins-faa",
        required=True,
        help="Protein FASTA with headers: contig::idx::start::end::strand",
    )
    ap.add_argument(
        "--out-affi", required=True, help="Output affi file (pipe-delimited)"
    )
    ap.add_argument(
        "--contigs-fasta",
        default=None,
        help="Contigs FASTA (sets contig lengths; ensures headers for all contigs)",
    )
    ap.add_argument(
        "--circular-list",
        default=None,
        help="File with circular contig IDs (one per line)",
    )
    ap.add_argument(
        "--labels-out",
        default=None,
        help="Optional per-ORF label table (not emitted here to keep code compact)",
    )
    ap.add_argument(
        "--coord-tolerance",
        type=int,
        default=5,
        help="±N nt allowed when matching Pharokka to proteins (default 5)",
    )
    args = ap.parse_args()

    build_affi(
        pharokka_tsv=args.pharokka,
        proteins_faa=args.proteins_faa,
        out_affi=args.out_affi,
        contigs_fasta=args.contigs_fasta,
        circular_list=args.circular_list,
        labels_out=args.labels_out,
        coord_tolerance=args.coord_tolerance,
    )


if __name__ == "__main__":
    main()
