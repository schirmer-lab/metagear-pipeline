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
    --sbatch-template     Path to submit_phold_predict.sbatch.in (templated).
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
    __SUGGESTED_WALLTIME__ e.g. "24:00:00" (capped at 48h — LRZ-HGX allocation max)
    __ETA_H100__, __ETA_H100_B16__, __ETA_A100__, __ETA_V100__, __ETA_T4__, __ETA_CPU__
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
    p.add_argument("--sbatch-template", required=True, type=Path)
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

    # ETA estimates per platform. Per-shard minutes scale roughly with
    # max_seq_count once the model load (~30-60 s per `phold proteins-predict`
    # invocation) is amortised. The figures below are calibrated for the
    # default shard size of 5000 proteins/shard and rescaled here so smaller
    # or larger shards get realistic projections without a model rewrite.
    scale = max(0.2, max_seq_count / 5000.0)
    eta_h100      = _eta_hint(n_shards, 5.0  * scale)   # H100 batch=32 (~4-6 min/shard at 5k)
    eta_h100_b16  = _eta_hint(n_shards, 6.5  * scale)   # H100 batch=16 (~5-8 min/shard at 5k)
    eta_a100      = _eta_hint(n_shards, 10.0 * scale)   # A100 batch=16 (~8-12 min/shard at 5k)
    eta_v100      = _eta_hint(n_shards, 16.0 * scale)   # V100 batch=8  (~12-20 min/shard at 5k)
    eta_t4        = _eta_hint(n_shards, 35.0 * scale)   # T4 / 3070 batch=2 (~25-45 min/shard at 5k)
    eta_cpu       = _eta_hint(n_shards, 300.0 * scale)  # CPU batch=1 (~4-6 h/shard at 5k)

    # Walltime suggestion: 2× the H100-batch=16 estimate (safety margin for
    # cold caches, contention, occasional retries), rounded up to the next
    # whole hour. Capped at 48h — the LRZ-HGX per-allocation maximum and a
    # reasonable ceiling on most GPU schedulers. Users can edit the sbatch
    # header if their site allows longer (or shorter) windows.
    LRZ_MAX_HOURS = 48
    raw_hours = max(2.0, n_shards * 6.5 * scale * 2.0 / 60.0)
    walltime_h = min(LRZ_MAX_HOURS, int(raw_hours + 0.999))   # ceil
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
        "__ETA_H100_B16__":       eta_h100_b16,
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
    sbatch = out_root / "submit_phold_predict.sbatch"
    sbatch.write_text(_fill(args.sbatch_template.read_text(), subs))
    sbatch.chmod(0o755)

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
    print(f"[package_phold_offsite]   suggested SLURM walltime: {suggested_walltime}")
    print(f"[package_phold_offsite]   estimated H100 batch=32 runtime: {eta_h100}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
