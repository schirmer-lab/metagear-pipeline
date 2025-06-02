#!/usr/bin/env python

import sys
import logging
from pathlib import Path

import click
import pandas as pd

# Configure logging
title = Path(__file__).name
logger = logging.getLogger(title)
logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")


@click.command(context_settings={"help_option_names": ["-h", "--help"]})
@click.argument("files", type=click.Path(exists=True), nargs=-1)
@click.option(
    "-o",
    "--output",
    type=click.Path(),
    default="merged_abundance.tsv",
    help="Output filename (TSV).",
)
def merge_abundance(files, output):
    """
    Merge CoverM contig abundance tables (count or rpkm) by concatenating sample columns from each batch file.

    FILES: list of partial abundance TSV files (all count or all rpkm).
    """
    if not files:
        logger.error("No input files provided. Use -h for help.")
        sys.exit(1)

    merged = None
    for f in files:
        path = Path(f)
        logger.info(f"Reading file: {f}")
        try:
            df = pd.read_csv(path, sep="\t", index_col=0)
        except Exception as e:
            logger.error(f"Failed to read {f}: {e}")
            sys.exit(1)

        # Merge by index, preserving all sample columns
        merged = df if merged is None else merged.join(df, how="outer")

    out_path = Path(output)
    logger.info(f"Writing merged output to {out_path}")
    try:
        merged.to_csv(out_path, sep="\t")
    except Exception as e:
        logger.error(f"Failed to write output: {e}")
        sys.exit(1)

    logger.info(f"Successfully merged {len(files)} files into '{output}'")


if __name__ == "__main__":
    merge_abundance()
