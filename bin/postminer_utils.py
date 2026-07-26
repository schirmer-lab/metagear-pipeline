#!/usr/bin/env python
import logging
import os
import sys
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

import click
import numpy as np
import pandas as pd
from pyfaidx import Fasta, FetchError

LOG_FORMAT = "%(asctime)s | %(levelname)s | %(message)s"
LOG_DATE_FORMAT = "%Y-%m-%d %H:%M:%S"


def _setup_logger(
    logger_name: str, log_file: Path | None = None, log_level: str = "INFO"
) -> logging.Logger:
    """Configure and return a logger that writes to stdout and optionally a file."""

    level = getattr(logging, log_level.upper(), logging.INFO)
    logger = logging.getLogger(logger_name)
    if logger.handlers:
        logger.setLevel(level)
        return logger

    logger.setLevel(level)
    logger.propagate = False

    stream_handler = logging.StreamHandler(stream=sys.stdout)
    stream_handler.setFormatter(logging.Formatter(LOG_FORMAT, LOG_DATE_FORMAT))
    logger.addHandler(stream_handler)

    if log_file:
        log_file = Path(log_file)
        log_file.parent.mkdir(parents=True, exist_ok=True)
        file_handler = logging.FileHandler(log_file, mode="w", encoding="utf-8")
        file_handler.setFormatter(logging.Formatter(LOG_FORMAT, LOG_DATE_FORMAT))
        logger.addHandler(file_handler)

    return logger


def _load_gene_map(all_msps_fp: str, sel_category: set[str]) -> pd.DataFrame:
    """Return a mapping dataframe containing gene_name -> msp_name for selected categories."""

    df = pd.read_csv(
        all_msps_fp,
        sep="\t",
        usecols=["msp_name", "gene_category", "gene_name"],
        dtype={"msp_name": "string", "gene_category": "string", "gene_name": "string"},
    )

    df = df[df["gene_category"].isin(sel_category)].dropna(
        subset=["gene_name", "msp_name"]
    )
    df = df.drop_duplicates(subset=["gene_name", "msp_name"])
    df = df.rename(columns={"gene_name": "gene_name", "msp_name": "msp_name"})
    return df[["gene_name", "msp_name"]]


def _read_csv_chunks(csv_path: str, chunk_size: int, **kwargs):
    """Yield chunks from a tabular file, preferring the pyarrow engine when available."""

    for engine in ("pyarrow", "c"):
        try:
            return pd.read_csv(csv_path, engine=engine, chunksize=chunk_size, **kwargs)
        except (
            ValueError
        ) as exc:  # pragma: no cover - engine availability is environment specific
            if "engine" in str(exc) and engine == "pyarrow":
                continue
            raise

    raise RuntimeError(f"Unable to create a CSV chunk iterator for {csv_path}")


@click.group("application")
def main():
    pass


@main.group("helper")
def helper():
    """
    Helper scripts for assembly-based metagenomic pipelines.
    """
    print("psot mspminer helper functions")
    pass


def _load_msp_gc_id(
    all_msps_fp, sel_category={"core", "accessory", "shared_core", "shared_accessory"}
):
    """
    Load gene_id lists for each msp_name, filtered by gene_category == 'core'.

    Parameters:
    - all_msps_fp (str): Path to the input TSV file "all_msps.tsv".

    Returns:
    - dict: Mapping of msp_name -> list of core gene_ids.
    """
    mapping_df = _load_gene_map(all_msps_fp, set(sel_category))
    return mapping_df.groupby("msp_name")["gene_name"].apply(list).to_dict()


def _calculate_sample_means(
    data_file, id_list, output_file=None, sep="\t", chunk_size=500000
):
    """
    Calculate column-wise means for selected rows from a large file.

    Parameters:
    - data_file (str): Path to the large data file.
    - id_list (list or set): List or set of row IDs to extract.
    - output_file (str): Optional path to write the result.
    - sep (str): Field separator (default: tab-delimited).
    - chunk_size (int): Use chunks to handle large files, adjust based on memory

    Returns:
    - pandas.Series: Mean values for each column (sample).
    """
    # Convert ID list to set for fast lookup
    target_ids = set(id_list)

    # Use chunks to handle large files

    chunks = pd.read_csv(data_file, sep=sep, index_col=0, chunksize=chunk_size)

    mean_sum = None
    count = 0

    for chunk in chunks:
        # Filter rows by target IDs
        filtered = chunk[chunk.index.isin(target_ids)]
        if not filtered.empty:
            if mean_sum is None:
                mean_sum = filtered.sum()
            else:
                mean_sum += filtered.sum()
            count += len(filtered)

    if count == 0:
        raise ValueError("No matching IDs found in the file.")

    mean_values = mean_sum / count

    if output_file:
        mean_values.to_csv(output_file, sep=sep)

    return mean_values


def _calculate_sample_medians(
    data_file, id_list, output_file=None, sep="\t", chunk_size=5000000
):
    """
    Calculate column-wise medians for selected rows from a large file.

    Parameters:
    - data_file (str): Path to the large data file.
    - id_list (list or set): List or set of row IDs to extract.
    - output_file (str): Optional path to write the result.
    - sep (str): Field separator (default: tab-delimited).
    - chunk_size (int): Chunk size for reading large files.

    Returns:
    - pandas.Series: Median values for each column (sample).
    """
    target_ids = set(id_list)
    collected_rows = []

    chunks = pd.read_csv(data_file, sep=sep, index_col=0, chunksize=chunk_size)

    for chunk in chunks:
        filtered = chunk[chunk.index.isin(target_ids)]
        if not filtered.empty:
            collected_rows.append(filtered)

    if not collected_rows:
        raise ValueError("No matching IDs found in the file.")

    # Concatenate all collected rows
    all_data = pd.concat(collected_rows)

    # Compute column-wise median
    median_values = all_data.median(axis=0)

    if output_file:
        median_values.to_csv(output_file, sep=sep)

    return median_values


def _extract_fasta_by_ids(
    fasta_source,
    header_ids,
    output_path,
    threads: int = 20,
    logger: logging.Logger | None = None,
):
    """Extract sequences from *fasta_source* for the provided IDs and write them to *output_path*."""

    owns_fasta = not isinstance(fasta_source, Fasta)
    fasta = (
        fasta_source
        if not owns_fasta
        else Fasta(fasta_source, rebuild=False, as_raw=True)
    )
    deduped_ids = list(dict.fromkeys(header_ids))
    max_workers = max(1, threads)
    missing_ids: list[str] = []

    def fetch_record(seq_id: str) -> str | None:
        try:
            seq = fasta[seq_id]
            return f">{seq.name}\n{str(seq)}\n"
        except FetchError:
            missing_ids.append(seq_id)
            return None

    with open(output_path, "w", encoding="utf-8") as out_f:
        if max_workers == 1:
            for seq_id in deduped_ids:
                record = fetch_record(seq_id)
                if record:
                    out_f.write(record)
        else:
            with ThreadPoolExecutor(max_workers=max_workers) as executor:
                futures = {
                    executor.submit(fetch_record, seq_id): seq_id
                    for seq_id in deduped_ids
                }
                for future in as_completed(futures):
                    record = future.result()
                    if record:
                        out_f.write(record)

    if missing_ids:
        message = (
            f"Missing {len(missing_ids)} sequence ids"
            if len(missing_ids) > 5
            else f"Missing sequence ids: {', '.join(missing_ids)}"
        )
        if logger:
            logger.warning(message)
        else:
            print(message)

    if owns_fasta and hasattr(fasta, "close"):
        fasta.close()


# use the median value of core genes to estimate the abundance of each MSPminer
@helper.command(name="get-msp-abd")
@click.option(
    "--rpkm-fp",
    default=False,
    type=str,
    help="file path to the gene abundance file [rpkm]",
)
@click.option(
    "--all-msps-fp",
    default=False,
    type=str,
    help="file path to the major mspminer output file, default name: [all_msps.tsv]",
)
@click.option("--save-fp", default=False, type=str, help="path to save the result")
@click.option(
    "--method",
    default="median",
    type=str,
    help="method for MSP abundance cacluation: [median or mean] or core genes",
)
@click.option(
    "--chunk-size",
    default=250000,
    show_default=True,
    type=int,
    help="Number of genes to load per chunk from the RPKM matrix.",
)
@click.option(
    "--dtype",
    default="float32",
    show_default=True,
    type=click.Choice(["float32", "float64"], case_sensitive=False),
    help="Floating-point precision used for intermediate calculations.",
)
@click.option(
    "--log-file",
    default=None,
    type=click.Path(dir_okay=False, resolve_path=True),
    help="Optional path to a log file for detailed progress output.",
)
@click.option(
    "--log-level",
    default="INFO",
    show_default=True,
    type=click.Choice(
        ["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"], case_sensitive=False
    ),
    help="Verbosity of log output.",
)
@click.option(
    "--progress-every",
    default=10,
    show_default=True,
    type=int,
    help="Emit a progress message every N processed chunks.",
)
def get_msp_abd(
    rpkm_fp,
    all_msps_fp,
    save_fp,
    method,
    chunk_size,
    dtype,
    log_file,
    log_level,
    progress_every,
):
    method = method.lower()
    if method not in {"median", "mean"}:
        raise click.ClickException("Invalid method. Choose either 'median' or 'mean'.")

    if not rpkm_fp or not all_msps_fp or not save_fp:
        raise click.ClickException(
            "'--rpkm-fp', '--all-msps-fp', and '--save-fp' are required inputs."
        )

    logger = _setup_logger(
        "postminer.get_msp_abd", Path(log_file) if log_file else None, log_level
    )
    logger.info("Starting MSP abundance calculation using the %s method", method)

    chunk_size = max(1, int(chunk_size))
    progress_every = max(1, int(progress_every))
    dtype = np.float32 if dtype.lower() == "float32" else np.float64

    logger.info("Loading MSP definitions from %s", all_msps_fp)
    gene_map_df = _load_gene_map(all_msps_fp, {"core"})
    gene_map_df["gene_name"] = gene_map_df["gene_name"].astype(str)
    gene_map_df["msp_name"] = gene_map_df["msp_name"].astype(str)

    if gene_map_df.empty:
        raise click.ClickException("No core genes found in the MSP definitions file.")

    core_gene_count = gene_map_df["gene_name"].nunique()
    total_msps = gene_map_df["msp_name"].nunique()
    logger.info(
        "Found %d MSPs covering %d unique core genes", total_msps, core_gene_count
    )

    gene_name_set = set(gene_map_df["gene_name"].tolist())
    msp_sums: dict[str, np.ndarray] = {}
    msp_counts: dict[str, np.ndarray] = {}
    msp_median_buffers: dict[str, list[np.ndarray]] = defaultdict(list)

    chunk_iter = _read_csv_chunks(rpkm_fp, chunk_size, sep="\t", index_col=0)
    sample_columns: list[str] | None = None
    total_matched_rows = 0

    for chunk_index, chunk in enumerate(chunk_iter, start=1):
        if sample_columns is None:
            sample_columns = chunk.columns.tolist()

        filtered_chunk = chunk.loc[chunk.index.intersection(gene_name_set)]
        del chunk

        if filtered_chunk.empty:
            if chunk_index == 1 or chunk_index % progress_every == 0:
                logger.info("Chunk %d: no overlapping genes", chunk_index)
            continue

        filtered_chunk = filtered_chunk.astype(dtype, copy=False)
        filtered_chunk = filtered_chunk.reset_index()
        index_col_name = filtered_chunk.columns[0]
        if index_col_name != "gene_name":
            filtered_chunk = filtered_chunk.rename(
                columns={index_col_name: "gene_name"}
            )
        merged = filtered_chunk.merge(
            gene_map_df, on="gene_name", how="inner", copy=False
        )
        matched_rows = merged.shape[0]
        total_matched_rows += matched_rows

        if chunk_index == 1 or chunk_index % progress_every == 0:
            logger.info(
                "Chunk %d: matched %d gene rows (running total %d)",
                chunk_index,
                matched_rows,
                total_matched_rows,
            )

        grouped = merged.groupby("msp_name", sort=False)
        for msp_id, frame in grouped:
            values = frame[sample_columns].to_numpy(dtype=dtype, copy=True)

            if method == "mean":
                totals = np.nansum(values, axis=0, dtype=np.float64)
                counts = np.sum(~np.isnan(values), axis=0, dtype=np.int64)

                if msp_id not in msp_sums:
                    msp_sums[msp_id] = np.zeros(len(sample_columns), dtype=np.float64)
                    msp_counts[msp_id] = np.zeros(len(sample_columns), dtype=np.int64)

                msp_sums[msp_id] += totals
                msp_counts[msp_id] += counts
            else:
                msp_median_buffers[msp_id].append(values)

        del filtered_chunk
        del merged

    if sample_columns is None:
        raise click.ClickException("The RPKM file appears to be empty or malformed.")

    if method == "mean":
        if not msp_sums:
            raise click.ClickException(
                "No overlapping core genes were found between the inputs."
            )

        result_matrix = {}
        for msp_id, totals in msp_sums.items():
            counts = msp_counts[msp_id]
            with np.errstate(divide="ignore", invalid="ignore"):
                mean_values = np.divide(
                    totals, counts, out=np.zeros_like(totals), where=counts != 0
                )
            result_matrix[msp_id] = mean_values
    else:
        if not msp_median_buffers:
            raise click.ClickException(
                "No overlapping core genes were found between the inputs."
            )

        result_matrix = {}
        for msp_id, matrices in msp_median_buffers.items():
            stacked = matrices[0] if len(matrices) == 1 else np.vstack(matrices)
            result_matrix[msp_id] = np.nanmedian(stacked, axis=0)

    result_df = pd.DataFrame.from_dict(
        result_matrix, orient="index", columns=sample_columns
    )
    result_df.index.name = "msp_name"
    result_df.sort_index(inplace=True)

    output_path = Path(save_fp)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    result_df.to_csv(output_path, sep="\t", index=True, header=True)

    logger.info(
        "Finished MSP abundance calculation for %d MSPs across %d samples. Output written to %s",
        result_df.shape[0],
        result_df.shape[1],
        output_path,
    )


# get pangenome sequences of each msp
@helper.command(name="get-msp-pangenome")
@click.option(
    "--gene-catalog-fp",
    default=False,
    type=str,
    help="file path to the gene catalog sequences [fasta]",
)
@click.option(
    "--all-msps-fp",
    default=False,
    type=str,
    help="file path to the major mspminer output file, default name: [all_msps.tsv]",
)
@click.option(
    "--msp-pangenome-dir",
    default=False,
    type=str,
    help="folder to save the output pangenome files [msp_id+.pangenome.fasta]",
)
@click.option(
    "--threads",
    default=20,
    show_default=True,
    type=int,
    help="Number of threads used to fetch sequences from the FASTA catalog.",
)
@click.option(
    "--log-file",
    default=None,
    type=click.Path(dir_okay=False, resolve_path=True),
    help="Optional path to a log file for detailed progress output.",
)
@click.option(
    "--log-level",
    default="INFO",
    show_default=True,
    type=click.Choice(
        ["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"], case_sensitive=False
    ),
    help="Verbosity of log output.",
)
@click.option(
    "--progress-every",
    default=100,
    show_default=True,
    type=int,
    help="Emit a progress message every N processed MSPs.",
)
def get_msp_pangenome(
    gene_catalog_fp,
    all_msps_fp,
    msp_pangenome_dir,
    threads,
    log_file,
    log_level,
    progress_every,
):
    if not gene_catalog_fp or not all_msps_fp or not msp_pangenome_dir:
        raise click.ClickException(
            "'--gene-catalog-fp', '--all-msps-fp', and '--msp-pangenome-dir' are required inputs."
        )

    logger = _setup_logger(
        "postminer.get_msp_pangenome", Path(log_file) if log_file else None, log_level
    )

    threads = max(1, int(threads))
    progress_every = max(1, int(progress_every))

    output_dir = Path(msp_pangenome_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    logger.info("Loading MSP gene map from %s", all_msps_fp)
    gene_map_df = _load_gene_map(
        all_msps_fp,
        {"core", "accessory", "shared_core", "shared_accessory"},
    )

    if gene_map_df.empty:
        raise click.ClickException(
            "No genes were found for the requested MSP categories."
        )

    gene_map_df["gene_name"] = gene_map_df["gene_name"].astype(str)
    gene_map_df["msp_name"] = gene_map_df["msp_name"].astype(str)
    gene_map_df = gene_map_df.sort_values("msp_name")

    total_msps = gene_map_df["msp_name"].nunique()
    total_genes = len(gene_map_df)
    logger.info(
        "Preparing pangenome FASTA files for %d MSPs covering %d genes using %d threads",
        total_msps,
        total_genes,
        threads,
    )

    fasta = Fasta(gene_catalog_fp, rebuild=False, as_raw=True)
    try:
        for idx, (msp_id, group_df) in enumerate(
            gene_map_df.groupby("msp_name", sort=False), start=1
        ):
            output_path = output_dir / f"{msp_id}.pangenome.fasta"
            _extract_fasta_by_ids(
                fasta,
                group_df["gene_name"].tolist(),
                output_path,
                threads=threads,
                logger=logger,
            )

            if idx == 1 or idx % progress_every == 0 or idx == total_msps:
                percent_complete = 100.0 * idx / total_msps
                logger.info(
                    "Processed %d/%d MSPs (%.1f%%)",
                    idx,
                    total_msps,
                    percent_complete,
                )
    finally:
        if hasattr(fasta, "close"):
            fasta.close()

    logger.info(
        "Completed MSP pangenome sequence extraction. Output folder: %s", output_dir
    )


if __name__ == "__main__":
    main()
