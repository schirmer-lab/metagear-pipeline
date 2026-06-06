#!/usr/bin/env python
"""Assemble the offsite-predict bundle for metagear `structures`.

Stages the protein shards + viral_join_table + a manifest + a templated
README and runner script into a single directory the user can rsync to a
GPU server. The GPU operator runs `bash run_phold_predict.sh --phold-db
<path>` there; the resulting predict outputs come back via rsync and
metagear resumes the workflow on CPU.

Inputs:
    --shards-dir          Directory containing the per-shard FASTAs
                          (e.g. *.faa.gz from SEQKIT_SPLIT2).
    --viral-join-table    viral_join_table.tsv to carry through to the
                          source machine's resume phase. Copied verbatim.
    --readme-template     Path to README.md.in (templated).
    --script-template     Path to run_phold_predict.sh.in (templated).
    --out-dir             Output bundle directory to assemble.

Template placeholders supported:
    __NUM_SHARDS__        e.g. 28
    __BUNDLE_SIZE_MB__    e.g. 6
    __SHARD_SIZE__        per-shard sequence count (taken from the largest
                          shard; the last shard may be smaller)
    __PHOLD_VERSION__     hardcoded from the pipeline (passed by caller via
                          --phold-version)
    __PIPELINE_COMMIT__   pipeline git commit (passed by caller, optional)
    __STRUCTURES_SCOPE__  the scope param used to build the bundle
    __SUGGESTED_WALLTIME__ e.g. "08:00:00"
    __ETA_H100__, __ETA_A100__, __ETA_V100__, __ETA_T4__, __ETA_CPU__
                          estimated total runtimes per platform
    __CREATED_AT__        ISO 8601 timestamp
"""

from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import os
import shutil
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, Iterable


def _count_sequences(path: Path) -> int:
    """Count '>' header lines in a (possibly gzipped) FASTA."""
    opener = gzip.open if str(path).endswith(".gz") else open
    n = 0
    with opener(path, "rt") as fh:
        for line in fh:
            if line.startswith(">"):
                n += 1
    return n


def _sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def _eta_hint(n_shards: int, per_shard_min: float, suffix: str = "") -> str:
    """Return a human-readable total-time string for the given per-shard rate."""
    total_min = n_shards * per_shard_min
    if total_min < 60:
        return f"{total_min:.0f} min{suffix}"
    elif total_min < 60 * 24:
        return f"{total_min/60:.1f} h{suffix}"
    else:
        return f"{total_min/60/24:.1f} d{suffix}"


def _fill(text: str, subs: Dict[str, str]) -> str:
    for key, value in subs.items():
        text = text.replace(key, value)
    return text


def main(argv: Iterable[str] | None = None) -> int:
    p = argparse.ArgumentParser(description=__doc__.split("\n\n", 1)[0])
    p.add_argument("--shards-dir", required=True, type=Path,
                   help="Directory containing the per-shard FASTAs (*.faa.gz)")
    p.add_argument("--viral-join-table", required=True, type=Path)
    p.add_argument("--readme-template", required=True, type=Path)
    p.add_argument("--script-template", required=True, type=Path)
    p.add_argument("--out-dir", required=True, type=Path)
    p.add_argument("--structures-scope", default="unknown")
    p.add_argument("--phold-version", default="1.2.5")
    p.add_argument("--pipeline-commit", default="unknown")
    args = p.parse_args(list(argv) if argv is not None else None)

    out_root: Path = args.out_dir
    out_shards = out_root / "shards"
    out_shards.mkdir(parents=True, exist_ok=True)

    # Stage shards as symlinks (rsync -L will dereference; saves space in
    # the work dir vs. copying gigabytes around inside Nextflow).
    shards = sorted(args.shards_dir.glob("*.faa.gz"))
    if not shards:
        print(f"ERROR: no *.faa.gz shards in {args.shards_dir}", file=sys.stderr)
        return 1

    manifest_shards = []
    total_bytes = 0
    max_seq_count = 0
    for shard in shards:
        dst = out_shards / shard.name
        if dst.exists() or dst.is_symlink():
            dst.unlink()
        # Hardlink if same FS, fall back to copy
        try:
            os.link(shard.resolve(), dst)
        except OSError:
            shutil.copy2(shard, dst)
        size = shard.stat().st_size
        nseq = _count_sequences(shard)
        max_seq_count = max(max_seq_count, nseq)
        total_bytes += size
        manifest_shards.append({
            "name": shard.name,
            "size_bytes": size,
            "n_sequences": nseq,
            "sha256": _sha256(shard),
            "expected_output_subdir": f"outputs/predict_{shard.stem.replace('.faa', '')}",
        })

    # Pass through the viral join table verbatim
    join_dst = out_root / "viral_join_table.tsv"
    shutil.copy2(args.viral_join_table, join_dst)

    n_shards = len(shards)
    bundle_size_mb = round(total_bytes / 1024 / 1024, 1)

    # ETA estimates per platform (per-shard minutes from prior PHOLD benchmarks)
    eta_h100 = _eta_hint(n_shards, 2.0)
    eta_a100 = _eta_hint(n_shards, 4.5)
    eta_v100 = _eta_hint(n_shards, 6.0)
    eta_t4   = _eta_hint(n_shards, 14.0)
    eta_cpu  = _eta_hint(n_shards, 90.0)

    # Walltime suggestion: 1.5× the slowest plausible GPU estimate, rounded up
    walltime_h = max(2, int(n_shards * 14.0 * 1.5 / 60) + 1)
    suggested_walltime = f"{walltime_h:02d}:00:00"

    subs = {
        "__NUM_SHARDS__":         str(n_shards),
        "__BUNDLE_SIZE_MB__":     f"{bundle_size_mb}",
        "__SHARD_SIZE__":         str(max_seq_count),
        "__PHOLD_VERSION__":      args.phold_version,
        "__PIPELINE_COMMIT__":    args.pipeline_commit,
        "__STRUCTURES_SCOPE__":   args.structures_scope,
        "__SUGGESTED_WALLTIME__": suggested_walltime,
        "__ETA_H100__":           eta_h100,
        "__ETA_A100__":           eta_a100,
        "__ETA_V100__":           eta_v100,
        "__ETA_T4__":             eta_t4,
        "__ETA_CPU__":            eta_cpu,
        "__CREATED_AT__":         datetime.now(timezone.utc).isoformat(timespec="seconds"),
    }

    # Render templates
    (out_root / "README.md").write_text(_fill(args.readme_template.read_text(), subs))
    runner = out_root / "run_phold_predict.sh"
    runner.write_text(_fill(args.script_template.read_text(), subs))
    runner.chmod(0o755)

    # Manifest — used by the source machine on resume to validate
    # completeness before resuming the workflow
    manifest = {
        "metagear_pipeline_commit": args.pipeline_commit,
        "phold_version":            args.phold_version,
        "structures_scope":         args.structures_scope,
        "created_at":               subs["__CREATED_AT__"],
        "n_shards":                 n_shards,
        "total_bytes":              total_bytes,
        "shards":                   manifest_shards,
        "expected_outputs_per_shard": [
            "phold_aa.fasta",
            "phold_3di.fasta",
            "phold_prostT5_3di_mean_probabilities.csv",
            "phold_prostT5_3di_all_probabilities.json",
        ],
    }
    (out_root / "manifest.json").write_text(json.dumps(manifest, indent=2))

    print(f"[package_phold_offsite] wrote {n_shards}-shard bundle to {out_root}")
    print(f"[package_phold_offsite]   total shard bytes: {bundle_size_mb} MB")
    print(f"[package_phold_offsite]   max sequences per shard: {max_seq_count}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
