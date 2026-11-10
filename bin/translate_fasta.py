#! /usr/bin/env python

import logging
from typing import List
import sys
import tqdm
from Bio import SeqIO
from Bio.Data.CodonTable import TranslationError
from Bio.SeqRecord import SeqRecord

# NCBI translation tables, tried in order: 11 is the Bacterial, Archaeal and Plant Plastid
# code, 4 the Mold/Protozoan/Coelenterate Mitochondrial and Mycoplasma/Spiroplasma code.
DEFAULT_TRANSLATION_TABLES = [11, 4]


def translate(fna: str, out: str, translation_tables: List[int]):
    """
    Translate a FASTA file with nucleotide sequences to amino acid sequences.
    """
    with open(out, "w") as faa:
        record: SeqRecord
        for record in tqdm.tqdm(SeqIO.parse(fna, "fasta")):
            failed = True
            for translation_table in translation_tables:
                try:
                    translated: SeqRecord = record.translate(
                        table=translation_table, stop_symbol="", cds=True
                    )
                    translated.id = record.id
                    translated.description = ""
                    SeqIO.write(translated, faa, "fasta")
                    failed = False
                    break
                except TranslationError:
                    pass

            if failed:
                logging.warning(f"Translation failed for {record.id}")


# usage:
# python /data/translate_fasta.py ${input_fp} ${output_fp}

if __name__ == "__main__":
    input_fp = sys.argv[1]  # fasta format DNA sequence
    sfp = sys.argv[2]  # fasta format DNA sequence
    translate(input_fp, sfp, DEFAULT_TRANSLATION_TABLES)
