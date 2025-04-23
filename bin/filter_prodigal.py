import gzip
from Bio import SeqIO

import logging
import sys
from pathlib import Path

import click

logger = logging.getLogger()


def process_prodigal_fasta(
    fasta_file: str,
    output_file: str,
    skip_partial: bool = True,
    remove_termination_marker: bool = True,
    convert_header: bool = True,
):
    with open(output_file, "w") as output:
        with gzip.open(fasta_file, "rt") as handle:
            for fasta_record in SeqIO.parse(handle, "fasta"):
                seq_id, start, end, strand, note = fasta_record.description.split(" # ")

                if convert_header:
                    tmp_lst = seq_id.split("_")

                    rel_pos = tmp_lst[-1]
                    contig = "_".join(tmp_lst[:-1])
                    strand_symbol = "+" if strand == "1" else "-"
                    new_id = f"{contig}::{rel_pos}::{start}::{end}::{strand_symbol}"
                    fasta_record.id = new_id
                    fasta_record.description = ""

                if skip_partial:
                    if "partial=00" not in note:
                        continue

                if remove_termination_marker:
                    if fasta_record.seq[-1] == "*":
                        fasta_record.seq = fasta_record.seq[:-1]

                SeqIO.write(fasta_record, output, "fasta")


@click.command()
@click.argument("fasta", type=click.Path(exists=True))
@click.argument("out", type=click.Path())
@click.option(
    "-l",
    "--log-level",
    default="WARNING",
    help="The desired log level (default WARNING).",
    type=click.Choice(
        ["CRITICAL", "ERROR", "WARNING", "INFO", "DEBUG"], case_sensitive=False
    ),
)
def main(fasta, out, log_level):
    """Coordinate argument parsing and program execution."""
    logging.basicConfig(level=log_level, format="[%(levelname)s] %(message)s")

    if not Path(fasta).is_file():
        logger.error(f"The given input file {fasta} was not found!")
        sys.exit(2)

    Path(out).parent.mkdir(parents=True, exist_ok=True)

    process_prodigal_fasta(
        fasta_file=fasta,
        output_file=out,
        skip_partial=True,
        remove_termination_marker=True,
        convert_header=True,
    )


if __name__ == "__main__":
    main()
