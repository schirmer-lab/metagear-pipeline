#!/usr/bin/env python

import logging
import sys
from pathlib import Path

import click
import pandas as pd
import plotly.express as px

logger = logging.getLogger()


def qc_summary_plot(file_in, file_out):
    stats = pd.read_csv(
        file_in,
        usecols=["Sample", "Trimming", "Contamination", "Final"],
        names=["Sample", "Raw", "Trimming", "Contamination", "Final"],
        header=0,
    ).sort_values(by=["Final"], ascending=False)
    fig = px.bar(
        stats,
        x="Sample",
        y=["Final", "Contamination", "Trimming"],
        labels={
            "variable": "Stage",
            "value": "Read count",
        },
    )
    fig.write_html(file_out)


@click.command()
@click.argument("file_in", type=click.Path(exists=True))
@click.argument("file_out", type=click.Path())
@click.option(
    "-l",
    "--log-level",
    default="WARNING",
    help="The desired log level (default WARNING).",
    type=click.Choice(
        ["CRITICAL", "ERROR", "WARNING", "INFO", "DEBUG"], case_sensitive=False
    ),
)
def main(file_in, file_out, log_level):
    """Coordinate argument parsing and program execution."""
    logging.basicConfig(level=log_level, format="[%(levelname)s] %(message)s")

    if not Path(file_in).is_file():
        logger.error(f"The given input file {file_in} was not found!")
        sys.exit(2)

    Path(file_out).parent.mkdir(parents=True, exist_ok=True)
    qc_summary_plot(file_in, file_out)


if __name__ == "__main__":
    main()
