#!/usr/bin/env python


"""Provide a command line tool to validate and transform tabular samplesheets."""

import click
import csv
import logging
import sys
import shutil
from pathlib import Path
from urllib.parse import urlparse

from abc import ABC, abstractmethod
import types

logger = logging.getLogger()

# Input types constants
INPUT_TYPES = types.SimpleNamespace()
INPUT_TYPES.CONTIG_READS = "contig_reads"
INPUT_TYPES.CONTIG = "contig"
INPUT_TYPES.BLAST_SEQUENCES = "blast_seqs"
INPUT_TYPES.READS = "reads"
INPUT_TYPES.GROUPED_READS = "grouped_reads"

# Conditions constants
CONDITIONS = types.SimpleNamespace()
CONDITIONS.Required = "Required"
CONDITIONS.SimpleString = "SimpleString"
CONDITIONS.FileExists = "FileExists"
CONDITIONS.FolderExists = "FolderExists"
CONDITIONS.Optional = "Optional"
CONDITIONS.ValidSequenceFormats = "ValidSequenceFormats"

# Optional Status
OPTIONAL = types.SimpleNamespace()
OPTIONAL.Present = "Present"
OPTIONAL.Absent = "Absent"


def create_condition(condition):
    match condition:
        case CONDITIONS.Required:
            return Required()
        case CONDITIONS.FileExists:
            return FileExists()
        case CONDITIONS.SimpleString:
            return Required()
        case CONDITIONS.Optional:
            return Optional()
        case CONDITIONS.FolderExists:
            return FolderExists()
        case _:
            return None


class Condition(ABC):
    def __init__(self):
        super().__init__()
        self._condition_type = ""
        self._result = False
        self._message = ""

    @abstractmethod
    def eval(self, value):
        return self._result, self._message


class Required(Condition):
    def __init__(self):
        super().__init__()
        self._condition_type = CONDITIONS.Required

    def eval(self, value):
        if len(value) > 0:
            self._result = True
            self._message = ""
        else:
            self._result = False
            self._message = "Field is required."

        return self._result, self._message


class Optional(Condition):
    def __init__(self):
        super().__init__()
        self._condition_type = CONDITIONS.Optional

    def eval(self, value):
        if len(value) > 0:
            self._result = True
            self._message = OPTIONAL.Present
        else:
            self._result = True
            self._message = OPTIONAL.Absent

        return self._result, self._message


class FileExists(Condition):
    def __init__(self):
        super().__init__()
        self._condition_type = CONDITIONS.FileExists

    def eval(self, value):
        # Skip if URL (Nextflow handles it)
        result = urlparse(value)
        if all([result.scheme, result.netloc]) or Path(value).is_file():
            self._result = True
            self._message = ""
        else:
            self._result = False
            self._message = f"File {value} doesn't exist."

        return self._result, self._message


class FolderExists(Condition):
    def __init__(self):
        super().__init__()
        self._condition_type = CONDITIONS.FolderExists

    def eval(self, value):
        if Path(value).is_dir():
            self._result = True
            self._message = ""
        else:
            self._result = False
            self._message = f"Directory {value} doesn't exist."

        return self._result, self._message


class FieldValidator:
    def __init__(self, field_name, src_conditions):
        self._field_name = field_name
        self._src_conditions = src_conditions
        self._conditions = []

        for src_condition in self._src_conditions:
            self._conditions.append(create_condition(src_condition))

    def validate(self, field_value):
        # Evaluate optional condition and only apply other conditions if value is present.
        conditions = self._conditions.copy()
        optional_condition = next(
            iter(
                filter(
                    lambda condition: condition._condition_type == CONDITIONS.Optional,
                    conditions,
                )
            ),
            None,
        )

        if optional_condition is not None:
            result, message = optional_condition.eval(field_value)
            # Don't apply conditions if value is absent (optional)
            if message == OPTIONAL.Absent:
                conditions = []

        # Evaluate conditions on value
        for condition in conditions:
            result, message = condition.eval(field_value)
            if not result:
                raise AssertionError(f"{self._field_name}: {message}")


class RowValidator:
    def __init__(self, field_validators):
        self._field_validators = field_validators

    def validate(self, row):
        for validator in self._field_validators:
            validator.validate(row[validator._field_name])


def read_head(handle, num_lines=10):
    """Read the specified number of lines from the current position in the file."""
    lines = []
    for idx, line in enumerate(handle):
        if idx == num_lines:
            break
        lines.append(line)
    return "".join(lines)


def sniff_format(handle):
    """
    Detect the tabular format.

    Args:
        handle (text file): A handle to a `text file`_ object. The read position is
        expected to be at the beginning (index 0).

    Returns:
        csv.Dialect: The detected tabular format.

    .. _text file:
        https://docs.python.org/3/glossary.html#term-text-file

    """
    peek = read_head(handle)
    handle.seek(0)
    sniffer = csv.Sniffer()
    dialect = sniffer.sniff(peek)
    return dialect


def check_input(file_in, validation_type, file_out):
    # Common validators
    sample_validator = FieldValidator(
        "sample", [CONDITIONS.Required, CONDITIONS.SimpleString]
    )
    fasta1_validator = FieldValidator(
        "fastq_1", [CONDITIONS.Required, CONDITIONS.FileExists]
    )
    fasta2_validator = FieldValidator(
        "fastq_2", [CONDITIONS.Optional, CONDITIONS.FileExists]
    )

    contig_validator = FieldValidator(
        "contig", [CONDITIONS.Required, CONDITIONS.FileExists]
    )
    group_validator = FieldValidator(
        "group", [CONDITIONS.Required, CONDITIONS.SimpleString]
    )
    tag_validator = FieldValidator(
        "tag", [CONDITIONS.Required, CONDITIONS.SimpleString]
    )

    analysis_validator = FieldValidator(
        "analysis", [CONDITIONS.Required, CONDITIONS.SimpleString]
    )
    search_sequence_validator = FieldValidator(
        "query_sequence", [CONDITIONS.Required, CONDITIONS.FileExists]
    )
    search_database_validator = FieldValidator(
        "search_database", [CONDITIONS.Required, CONDITIONS.FileExists]
    )

    row_validator = None

    match validation_type:
        case INPUT_TYPES.READS:
            row_validator = RowValidator(
                [sample_validator, fasta1_validator, fasta2_validator]
            )
        case INPUT_TYPES.CONTIG_READS:
            row_validator = RowValidator(
                [sample_validator, contig_validator, fasta1_validator, fasta2_validator]
            )
        case INPUT_TYPES.CONTIG:
            row_validator = RowValidator([sample_validator, contig_validator])
        case INPUT_TYPES.BLAST_SEQUENCES:
            row_validator = RowValidator(
                [
                    analysis_validator,
                    search_sequence_validator,
                    search_database_validator,
                ]
            )
        case INPUT_TYPES.GROUPED_READS:
            row_validator = RowValidator(
                [
                    sample_validator,
                    group_validator,
                    tag_validator,
                    fasta1_validator,
                    fasta2_validator,
                ]
            )
        case _:
            row_validator = None
            raise Exception("Validation type not supported...")

    required_columns = set()
    for field_validator in row_validator._field_validators:
        if CONDITIONS.Required in field_validator._src_conditions:
            required_columns.add(field_validator._field_name)

    # See https://docs.python.org/3.9/library/csv.html#id3 to read up on `newline=""`.
    with file_in.open(newline="") as in_handle:
        reader = csv.DictReader(in_handle, dialect=sniff_format(in_handle))
        # Validate the existence of the expected header columns.
        if not required_columns.issubset(reader.fieldnames):
            req_cols = ", ".join(required_columns)
            logger.critical(
                f"The sample sheet **must** contain these column headers: {req_cols}."
            )
            sys.exit(1)
        # Validate each row.
        for i, row in enumerate(reader):
            try:
                row_validator.validate(row)
            except AssertionError as error:
                logger.critical(f"{str(error)} On line {i + 2}.")
                sys.exit(1)

    # Copy file if there are not error in validation
    shutil.copy(file_in, file_out)


@click.command()
@click.option("--input", "-i", type=click.Path(exists=True), help="Input file.")
@click.option(
    "--validation_type",
    "-t",
    type=click.Choice(
        ["reads", "grouped_reads", "contig_reads", "contig", "blast_seqs"]
    ),
    default="reads",
    help="Validation type.",
)
@click.option("--output", "-o", type=click.Path(), help="Output file.")
def main(input, validation_type, output):
    """Coordinate argument parsing and program execution."""
    if not Path(input).is_file():
        logger.error(f"The given input file {input} was not found!")
        sys.exit(2)
    Path(output).parent.mkdir(parents=True, exist_ok=True)
    check_input(Path(input), validation_type, output)


if __name__ == "__main__":
    main()
