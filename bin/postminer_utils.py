#!/usr/bin/env python
import click
import pandas as pd
from collections import defaultdict
import os
from pyfaidx import Fasta, FetchError
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

from sklearn import datasets, linear_model
from sklearn.metrics import mean_squared_error, r2_score


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
    # Read only necessary columns
    df = pd.read_csv(
        all_msps_fp,
        sep="\t",
        usecols=["msp_name", "gene_category", "gene_id", "gene_name"],
    )

    # Filter the selected category
    sel_df = df[df["gene_category"].isin(sel_category)]

    # Group by msp_name and aggregate gene_id into lists
    gene_dict = sel_df.groupby("msp_name")["gene_name"].apply(list).to_dict()

    return gene_dict


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


def _extract_fasta_by_ids(fasta_path, header_ids, output_path, threads=20):
    """
    Efficiently extract sequences from a FASTA file by header ID using multithreading.

    Parameters:
    - fasta_path (str or Path): Path to the input FASTA file.
    - header_ids (list or set): Sequence IDs to extract.
    - output_path (str or Path): File to write extracted sequences.
    - threads (int): Number of threads to use.
    """
    fasta = Fasta(fasta_path, rebuild=False, as_raw=True)
    header_ids = list(set(header_ids))  # Remove duplicates and allow indexing

    def fetch_record(seq_id):
        try:
            seq = fasta[seq_id]
            return f">{seq.name}\n{str(seq)}\n"
        except FetchError:
            return f"[Warning] ID not found: {seq_id}\n"

    with ThreadPoolExecutor(max_workers=threads) as executor:
        futures = {
            executor.submit(fetch_record, seq_id): seq_id for seq_id in header_ids
        }

        with open(output_path, "w") as out_f:
            for future in as_completed(futures):
                record = future.result()
                if not record.startswith("[Warning]"):
                    out_f.write(record)
                else:
                    print(record.strip())


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
def get_msp_abd(rpkm_fp, all_msps_fp, save_fp, method):
    # check if input method is valid
    if method not in {"median", "mean"}:
        print("invalid method: {0}, please select from [median or mean]".format(method))
        return

    core_gene_dic = _load_msp_gc_id(all_msps_fp, sel_category={"core"})
    msp_id_lst = list(core_gene_dic.keys())
    msp_id_lst.sort()

    merged_abd = pd.DataFrame()
    for ii, cur_msp_id in enumerate(msp_id_lst):
        if ii % 100 == 0:
            print(ii, round(100 * ii / len(msp_id_lst)))
        # cur_abd: a Series or single-column DataFrame with sample IDs as index
        if method == "median":
            cur_abd = _calculate_sample_medians(rpkm_fp, core_gene_dic[cur_msp_id])
        elif method == "mean":
            cur_abd = _calculate_sample_means(rpkm_fp, core_gene_dic[cur_msp_id])
        else:
            print(
                "invalid method: {0}, please select from [median or mean]".format(
                    method
                )
            )
            return

        # Ensure cur_abd is a Series and name it with msp_id
        if isinstance(cur_abd, pd.DataFrame):
            cur_abd = cur_abd.iloc[:, 0]
        cur_abd.name = cur_msp_id

        # Concatenate along columns (axis=1)
        merged_abd = pd.concat([merged_abd, cur_abd], axis=1)

    rotated_merged_abd = merged_abd.T

    # Save to file, including row and column names
    rotated_merged_abd.to_csv(save_fp, sep="\t", index=True, header=True)
    return


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
def get_msp_pangenome(gene_catalog_fp, all_msps_fp, msp_pangenome_dir):
    all_gene_dic = _load_msp_gc_id(
        all_msps_fp,
        sel_category={"core", "accessory", "shared_core", "shared_accessory"},
    )
    msp_id_lst = list(all_gene_dic.keys())
    msp_id_lst.sort()
    for ii, cur_msp in enumerate(msp_id_lst):
        if ii % 100 == 0:
            print("{}%% done".format(round(100.0 * ii / len(msp_id_lst))))
        cur_sfp = os.path.join(msp_pangenome_dir, cur_msp + ".pangenome.fasta")
        cur_gc_lst = all_gene_dic[cur_msp]
        _extract_fasta_by_ids(gene_catalog_fp, cur_gc_lst, cur_sfp)
    return


if __name__ == "__main__":
    main()
    # all_msps_fp = "/nfs/arxiv/shen/CLD_KCH_2025/analysis/gene_profile/results/mspminer/raw/all_msps.tsv"
    # rpkm_fp = "/nfs/arxiv/shen/CLD_KCH_2025/analysis/gene_profile/results/merge/gene_profile_rpkm_merged.tsv"
    # msp_abd_sfp = "/nfs/arxiv/shen/CLD_KCH_2025/analysis/gene_profile/results/mspminer/msp_abundance.median.RPKM.txt"
    # gc_seq_fp = "/nfs/arxiv/shen/CLD_KCH_2025/analysis/gene_call/results/cdhit/merged_genes.nr_95_90.fa"
    # msp_pangenome_dir = "/nfs/arxiv/shen/CLD_KCH_2025/analysis/gene_profile/results/mspminer/pangenome"
