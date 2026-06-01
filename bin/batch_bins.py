#!/usr/bin/env python
"""Partition cohort bins into dRep batches via connected components.

Reads a sparse skani triangle TSV (Ref_file, Query_file, ANI, ...) and the
cohort genomeInfo CSV (genome, completeness, contamination), and emits one
directory per batch under the working directory:

    batch_<NNN>/
        bins/<bin>.fa               # symlinks to the per-batch bin FASTAs
        drep_work_seed/             # mirrors STAGE_DREP_WORK's layout so the
            genomeInfo.csv          #   downstream dRep call can reuse the
                                    #   same `--genomeInfo` ext.args path.

A pair (A, B) shares an edge iff their skani ANI is >= --ani-threshold
(default 0.90). Connected components of the resulting graph define the
multi-bin batches. Size-1 components (bins with no neighbours above the
threshold) are merged into a single "singletons" batch — dRep runs on it
just like any other batch, so the -comp/-con quality filter is applied
uniformly across the cohort. Because every bin in the singletons batch is
<--ani-threshold ANI from every other singleton, dRep's primary clustering
puts each in its own primary cluster (no secondary clustering happens) and
each emerges as its own winner. This is semantically equivalent to running
dRep once per singleton, but bundled so dRep doesn't error on N=1 input
(which it does: `linkage` can't operate on an empty distance matrix).

Edge case: if there is exactly ONE size-1 component, we attach it to the
smallest multi-bin batch instead (a singletons batch of N=1 would re-
introduce the same dRep N=1 failure). The attached singleton is <ANI
threshold from every bin in the host batch, so it still emerges as its own
cluster — the result is unchanged, only the batch boundary moves.

Threshold safety:
  dRep's secondary clustering threshold defaults to 0.95 ANI. Any pair at
  ANI >= 0.95 will, by single-linkage at 0.90, share a component (since
  0.95 >= 0.90 means the edge exists). So within-batch dRep at 0.95 yields
  the same cohort-wide representatives as a single all-vs-all dRep on the
  full cohort. The 0.05 ANI buffer (0.90 vs 0.95) is conservatism, not a
  correctness requirement — same tool (skani) on both sides means there is
  no measurement-noise gap to budget for.

skani triangle output uses bin FASTA paths in its first two columns. We
key the graph on basenames (`<sample>_binN.fa`) so matching against
genomeInfo.csv's `genome` column (also `<sample>_binN.fa`) is direct.
"""

from __future__ import annotations

import argparse
import csv
import logging
import os
import sys
from collections import defaultdict
from pathlib import Path
from typing import Dict, Iterable, List

LOG = logging.getLogger("batch_bins")


def _load_genome_info(path: Path) -> Dict[str, List[str]]:
    """Load genomeInfo.csv as {genome_basename: [genome, completeness, contam]}.

    Preserves the original row text so we can write the sliced CSVs back
    out unchanged (and don't have to worry about float formatting).
    """
    rows: Dict[str, List[str]] = {}
    with path.open(newline="") as fh:
        reader = csv.reader(fh)
        header = next(reader)
        if [c.strip().lower() for c in header[:3]] != ["genome", "completeness", "contamination"]:
            raise ValueError(
                f"genomeInfo.csv header must start with genome,completeness,contamination — got {header!r}"
            )
        for row in reader:
            if not row or not row[0]:
                continue
            rows[row[0]] = row
    LOG.info("loaded %d genomes from %s", len(rows), path)
    return rows


def _index_bins(bin_paths: Iterable[Path]) -> Dict[str, Path]:
    """Map bin basename -> Path. The basename matches the genomeInfo `genome`
    column and the skani triangle Ref_file/Query_file basename."""
    idx: Dict[str, Path] = {}
    for p in bin_paths:
        idx[p.name] = p
    LOG.info("indexed %d bin FASTAs", len(idx))
    return idx


def _build_components(
    ani_tsv: Path, threshold: float, all_bins: Iterable[str]
) -> List[List[str]]:
    """Union-find over the sparse ANI TSV at the given threshold.

    Nodes are bin basenames. Every bin in `all_bins` is seeded so that
    singletons (no edges) end up as their own component.
    """
    parent: Dict[str, str] = {b: b for b in all_bins}

    def find(x: str) -> str:
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    def union(a: str, b: str) -> None:
        ra, rb = find(a), find(b)
        if ra != rb:
            parent[ra] = rb

    n_edges = 0
    n_skipped_unknown = 0
    with ani_tsv.open() as fh:
        first = fh.readline()
        if not first.startswith("Ref_file"):
            fh.seek(0)
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 3:
                continue
            ref, qry, ani_str = parts[0], parts[1], parts[2]
            try:
                ani = float(ani_str)
            except ValueError:
                continue
            # skani reports ANI as a percentage (e.g. 95.3); convert to fraction.
            if ani > 1.0:
                ani = ani / 100.0
            if ani < threshold:
                continue
            ref_b, qry_b = os.path.basename(ref), os.path.basename(qry)
            if ref_b not in parent or qry_b not in parent:
                n_skipped_unknown += 1
                continue
            union(ref_b, qry_b)
            n_edges += 1
    LOG.info("kept %d edges at ANI >= %.3f (skipped %d unknown-bin edges)",
             n_edges, threshold, n_skipped_unknown)

    components: Dict[str, List[str]] = defaultdict(list)
    for b in parent:
        components[find(b)].append(b)
    out = list(components.values())
    out.sort(key=lambda c: (-len(c), c[0]))
    LOG.info("built %d components (largest=%d, singletons=%d)",
             len(out),
             len(out[0]) if out else 0,
             sum(1 for c in out if len(c) == 1))
    return out


def _bundle_singletons(components: List[List[str]]) -> List[List[str]]:
    """Collapse all size-1 components into a single 'singletons' batch.

    Rationale: dRep crashes on N=1 input (scipy linkage can't run on an
    empty distance matrix). By bundling all singletons into one batch, dRep
    sees N>=2 bins that are all <ani_threshold from each other; its primary
    clustering puts each in its own primary cluster, no secondary clustering
    happens, the -comp/-con filter still applies, and each bin emits as
    its own winner — semantically equivalent to per-singleton dRep.

    Edge case: if there is exactly ONE singleton across the cohort, the
    singletons batch would itself be N=1. We attach it to the smallest
    multi-bin batch instead, which is safe because the attached singleton
    is <ani_threshold from every host bin (so dRep keeps it as its own
    cluster).
    """
    multi = [c for c in components if len(c) > 1]
    singletons = [c[0] for c in components if len(c) == 1]
    if not singletons:
        return multi
    if len(singletons) == 1:
        if not multi:
            # No multi-batch to attach to. Cohort is degenerate (effectively
            # one bin). Surface this loudly — dRep would fail downstream.
            raise ValueError(
                "Cannot batch a cohort of exactly one bin: no multi-bin "
                "batches to attach the lone singleton to."
            )
        # Append to the smallest multi-batch so we perturb the run-time
        # distribution as little as possible.
        host = min(multi, key=len)
        host.append(singletons[0])
        LOG.info("attached lone singleton %r to smallest multi-batch (size now %d)",
                 singletons[0], len(host))
        return multi
    LOG.info("bundled %d singletons into one batch", len(singletons))
    return multi + [singletons]


def _link_bin(src: Path, dst: Path) -> None:
    """Hardlink src → dst, falling back to symlink if the FS doesn't allow
    hardlinks (e.g. across mounts). Resolves src first so symlink fallbacks
    point at the real file instead of a chain of symlinks."""
    src_real = src.resolve()
    try:
        os.link(src_real, dst)
    except OSError:
        os.symlink(src_real, dst)


def _write_batches(
    components: List[List[str]],
    genome_info: Dict[str, List[str]],
    bins_by_name: Dict[str, Path],
    out_root: Path,
    header: List[str],
) -> int:
    out_root.mkdir(parents=True, exist_ok=True)
    width = max(3, len(str(len(components))))
    for i, comp in enumerate(components, start=1):
        batch_dir = out_root / f"batch_{i:0{width}d}"
        batch_bins = batch_dir / "bins"
        batch_seed = batch_dir / "drep_work_seed"
        batch_bins.mkdir(parents=True, exist_ok=True)
        batch_seed.mkdir(parents=True, exist_ok=True)
        # Per-batch genomeInfo (dRep's --genomeInfo target).
        with (batch_seed / "genomeInfo.csv").open("w", newline="") as fh:
            writer = csv.writer(fh)
            writer.writerow(header)
            for b in sorted(comp):
                row = genome_info.get(b)
                if row is None:
                    LOG.warning("no genomeInfo row for %s — emitting placeholder", b)
                    writer.writerow([b, "0", "0"])
                else:
                    writer.writerow(row)
        # Stage the per-batch bin FASTAs (hardlink, fallback symlink).
        for b in sorted(comp):
            src = bins_by_name.get(b)
            if src is None:
                LOG.warning("no FASTA found for %s — batch will be missing this bin", b)
                continue
            _link_bin(src, batch_bins / b)
    return len(components)


def _parse_args(argv: Iterable[str] | None = None) -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__.split("\n", 1)[0])
    p.add_argument("--ani-tsv", required=True, type=Path,
                   help="Sparse skani triangle output (Ref_file, Query_file, ANI, ...)")
    p.add_argument("--genome-info", required=True, type=Path,
                   help="Cohort genomeInfo CSV (genome, completeness, contamination)")
    p.add_argument("--bins", required=True, nargs="+", type=Path,
                   help="All cohort bin FASTAs. Basename must match the genome "
                        "column of --genome-info and the basenames in --ani-tsv.")
    p.add_argument("--ani-threshold", type=float, default=0.90,
                   help="Single-linkage edge threshold as fraction (default 0.90). Must be "
                        "strictly less than dRep's secondary clustering threshold (default 0.95).")
    p.add_argument("--out-root", type=Path, default=Path("batches"),
                   help="Directory to write batch_<NNN>/ subdirs into.")
    p.add_argument("--log-level", default="INFO", choices=["DEBUG", "INFO", "WARNING", "ERROR"])
    return p.parse_args(list(argv) if argv is not None else None)


def main(argv: Iterable[str] | None = None) -> int:
    args = _parse_args(argv)
    logging.basicConfig(level=args.log_level, format="%(asctime)s | %(levelname)s | %(message)s")

    if not (0 < args.ani_threshold < 0.95):
        LOG.error("--ani-threshold must be in (0, 0.95); got %.3f", args.ani_threshold)
        return 1

    genome_info = _load_genome_info(args.genome_info)
    if not genome_info:
        LOG.error("no rows loaded from genomeInfo CSV")
        return 1

    bins_by_name = _index_bins(args.bins)
    missing = [b for b in genome_info if b not in bins_by_name]
    if missing:
        LOG.warning("%d genomeInfo entries have no matching FASTA (e.g. %s)",
                    len(missing), missing[:3])

    header = ["genome", "completeness", "contamination"]
    components = _build_components(args.ani_tsv, args.ani_threshold, genome_info.keys())
    try:
        components = _bundle_singletons(components)
    except ValueError as e:
        LOG.error(str(e))
        return 1
    n_batches = _write_batches(components, genome_info, bins_by_name, args.out_root, header)
    LOG.info("wrote %d batches to %s", n_batches, args.out_root)
    return 0


if __name__ == "__main__":
    sys.exit(main())
