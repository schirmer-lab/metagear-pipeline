#!/usr/bin/env python
"""
Standalone script to merge Genomad and CheckV outputs into a single summary table.
"""
import os
import click
import pandas as pd


def find_closest_match(taxonomy, ictv_df):
    taxonomy_levels = taxonomy.split(";")
    matching_rows = []

    for i in range(len(taxonomy_levels), 0, -1):
        partial_taxonomy = ";".join(taxonomy_levels[:i])
        match_rows = ictv_df[ictv_df["Family"] == partial_taxonomy]
        if not match_rows.empty:
            matching_rows.extend(match_rows.iterrows())

    if matching_rows:
        return matching_rows[0][1][["Genome type", "Host type"]]
    else:
        return pd.Series({"Genome type": "Unknown", "Host type": "Unknown"})


def taxon_to_genome_type(taxonomy, ictv_taxonomy_df):
    genome_host_info = find_closest_match(taxonomy, ictv_taxonomy_df)
    return genome_host_info["Genome type"]


def taxon_to_host_type(taxonomy, ictv_taxonomy_df):
    genome_host_info = find_closest_match(taxonomy, ictv_taxonomy_df)
    return genome_host_info["Host type"]


def rename_sequence(sequence_name):
    # sequence_name = sequence_name.replace("provirus_", "")
    # sequence_name = sequence_name.rsplit("/", 1)[0]
    # sequence_name = (
    #     sequence_name.replace("|", "_")
    #     .replace("/", "_")
    #     .replace(":", "_")
    #     # .replace("-", "_")
    # )
    return sequence_name


@click.command()
@click.option(
    "--sample-name", "-s", required=True, help="Name of the sample (used in output)."
)
@click.option(
    "--viral-checkv",
    required=True,
    type=click.Path(exists=True),
    help="Path to the virus CheckV quality_summary.tsv file.",
)
@click.option(
    "--provirus-checkv",
    required=False,
    type=click.Path(exists=False),
    default=None,
    help="Path to the provirus CheckV quality_summary.tsv file.",
)
@click.option(
    "--viral-genomad",
    required=True,
    type=click.Path(exists=True),
    help="Path to the virus Genomad summary TSV file.",
)
@click.option(
    "--provirus-genomad",
    required=False,
    type=click.Path(exists=False),
    default=None,
    help="Path to the provirus Genomad summary TSV file.",
)
@click.option(
    "--ictv-taxonomy",
    required=True,
    type=click.Path(exists=True),
    help="Full path to the ICTV_Taxonomy_List.tsv file.",
)
@click.option(
    "--viral-min-genes",
    type=int,
    default=1,
    show_default=True,
    help="Minimum number of viral genes for filtering.",
)
@click.option(
    "--host-viral-genes-ratio",
    type=int,
    default=1,
    show_default=True,
    help="Maximum ratio of host to viral genes for filtering.",
)
@click.option(
    "--output-file",
    "-o",
    required=True,
    type=click.Path(),
    help="Path to write the merged output TSV.",
)
def main(
    sample_name,
    viral_checkv,
    provirus_checkv,
    viral_genomad,
    provirus_genomad,
    ictv_taxonomy,
    viral_min_genes,
    host_viral_genes_ratio,
    output_file,
):
    """
    Merge CheckV and Genomad outputs for viruses and proviruses, annotate genome and host types, and write a merged summary TSV.
    Also apply filtering based on viral_genes and host_genes/viral_genes ratio.
    """
    # Load ICTV taxonomy
    ictv_df = pd.read_csv(ictv_taxonomy, sep="\t")

    # Load CheckV tables
    checkv_cols = [
        "contig_id",
        "contig_length",
        "provirus",
        "proviral_length",
        "gene_count",
        "viral_genes",
        "host_genes",
        "checkv_quality",
        "miuvig_quality",
        "completeness",
        "completeness_method",
        "kmer_freq",
    ]

    df_virus_CHECKV = pd.read_csv(
        viral_checkv,
        sep="\t",
        usecols=checkv_cols,
    )

    if provirus_checkv and os.path.exists(provirus_checkv) and os.path.getsize(provirus_checkv) > 0:
        df_provirus_CHECKV = pd.read_csv(
            provirus_checkv,
            sep="\t",
            usecols=checkv_cols,
        )
    else:
        click.echo("Provirus CheckV file missing or empty. Continuing without it.")
        df_provirus_CHECKV = pd.DataFrame(columns=checkv_cols)

    # Filter out provirus lines from virus CHECKV
    df_virus_CHECKV = df_virus_CHECKV[df_virus_CHECKV["provirus"] != "Yes"]

    # Merge virus and provirus CheckV DataFrames
    checkv_merged_df = pd.concat(
        [df_virus_CHECKV, df_provirus_CHECKV], ignore_index=True
    )

    # Load Genomad tables
    genomad_cols = [
        "seq_name",
        "topology",
        "n_genes",
        "genetic_code",
        "virus_score",
        "n_hallmarks",
        "marker_enrichment",
        "taxonomy",
        "coordinates",
    ]

    df_virus_GENOMAD = pd.read_csv(
        viral_genomad,
        sep="\t",
        usecols=genomad_cols,
    )

    if provirus_genomad and os.path.exists(provirus_genomad) and os.path.getsize(provirus_genomad) > 0:
        df_provirus_GENOMAD = pd.read_csv(
            provirus_genomad,
            sep="\t",
            usecols=genomad_cols,
        )
    else:
        click.echo("Provirus Genomad file missing or empty. Continuing without it.")
        df_provirus_GENOMAD = pd.DataFrame(columns=genomad_cols)

    # Merge virus and provirus Genomad DataFrames
    genomad_merged_df = pd.concat(
        [df_virus_GENOMAD, df_provirus_GENOMAD], ignore_index=True
    )

    # Merge CheckV and Genomad on contig_id = seq_name
    merged_df = pd.merge(
        checkv_merged_df,
        genomad_merged_df,
        left_on="contig_id",
        right_on="seq_name",
        how="left",
    )
    # Rename and drop redundant columns
    merged_df = merged_df.rename(
        columns={"contig_id": "virus_id", "contig_length": "virus_length"}
    )
    if "proviral_length" in merged_df.columns:
        merged_df = merged_df.drop(columns=["proviral_length"])

    merged_df["topology"] = merged_df["topology"].astype(str)
    merged_df["virus_id"] = merged_df["virus_id"].astype(str)

    # Copy n_genes to gene_count if present
    if "n_genes" in merged_df.columns:
        merged_df["gene_count"] = merged_df["n_genes"]
        merged_df = merged_df.drop(columns=["n_genes"])

    # Update provirus flag based on topology
    mask_provirus = merged_df["topology"].str.contains("Provirus", na=False)
    merged_df.loc[mask_provirus, "provirus"] = "Yes"
    merged_df.loc[merged_df["provirus"] == "Yes", "topology"] = merged_df.loc[
        merged_df["provirus"] == "Yes", "topology"
    ].replace("No terminal repeats", "Provirus")

    merged_df = merged_df.fillna("NA")

    # Annotate genome and host types
    merged_df["Genome type"] = merged_df["taxonomy"].apply(
        lambda x: taxon_to_genome_type(x, ictv_df)
    )
    merged_df["Host type"] = merged_df["taxonomy"].apply(
        lambda x: taxon_to_host_type(x, ictv_df)
    )

    # Exclude unwanted columns
    columns_to_exclude = ["seq_name", "topology"]
    merged_df = merged_df.drop(
        columns=[col for col in columns_to_exclude if col in merged_df.columns]
    )

    # Add Sample column
    merged_df.insert(0, "Sample", sample_name)

    # Rename virus_id values
    merged_df["virus_id"] = merged_df["virus_id"].apply(rename_sequence)

    # Ensure output directory exists
    output_dir = os.path.dirname(output_file)
    if output_dir and not os.path.exists(output_dir):
        os.makedirs(output_dir)

    # Write unfiltered TSV
    merged_df.to_csv(output_file, sep="\t", index=False)
    click.echo(f"Merged summary written to {output_file}")

    # Apply filtering based on viral_genes and host/viral ratio
    # Ensure numeric types for filtering
    merged_df["viral_genes"] = (
        pd.to_numeric(merged_df["viral_genes"], errors="coerce").fillna(0).astype(int)
    )
    merged_df["host_genes"] = (
        pd.to_numeric(merged_df["host_genes"], errors="coerce").fillna(0).astype(int)
    )

    filtered_df = merged_df.loc[merged_df["viral_genes"] >= viral_min_genes]
    filtered_df = filtered_df.loc[
        (filtered_df["host_genes"] / filtered_df["viral_genes"])
        <= host_viral_genes_ratio
    ]

    # Construct filtered output file name
    root, ext = os.path.splitext(output_file)
    filtered_output_file = root + ".filtered" + ext

    # Write filtered TSV
    filtered_df.to_csv(filtered_output_file, sep="\t", index=False)
    click.echo(f"Filtered summary written to {filtered_output_file}")


if __name__ == "__main__":
    main()
