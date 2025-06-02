#! /usr/bin/env python

import logging
from typing import List
import sys
import tqdm
from Bio import SeqIO
from Bio.Data.CodonTable import TranslationError
from Bio.SeqRecord import SeqRecord

DEFAULT_TRANSLATION_TABLES = [11, 4]  # default translation table for microbiome data

"""
1. The Standard Code
2. The Vertebrate Mitochondrial Code
3. The Yeast Mitochondrial Code
4. The Mold, Protozoan, and Coelenterate Mitochondrial Code and the Mycoplasma/Spiroplasma Code
5. The Invertebrate Mitochondrial Code
6. The Ciliate, Dasycladacean and Hexamita Nuclear Code
9. The Echinoderm and Flatworm Mitochondrial Code
10. The Euplotid Nuclear Code
11. The Bacterial, Archaeal and Plant Plastid Code
12. The Alternative Yeast Nuclear Code
13. The Ascidian Mitochondrial Code
14. The Alternative Flatworm Mitochondrial Code
15. Blepharisma Nuclear Code
16. Chlorophycean Mitochondrial Code
21. Trematode Mitochondrial Code
22. Scenedesmus obliquus Mitochondrial Code
23. Thraustochytrium Mitochondrial Code
24. Rhabdopleuridae Mitochondrial Code
25. Candidate Division SR1 and Gracilibacteria Code
26. Pachysolen tannophilus Nuclear Code
27. Karyorelict Nuclear Code
28. Condylostoma Nuclear Code
29. Mesodinium Nuclear Code
30. Peritrich Nuclear Code
31. Blastocrithidia Nuclear Code
32. Balanophoraceae Plastid Code
33. Cephalodiscidae Mitochondrial UAA-Tyr Code
"""


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
