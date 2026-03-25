"""
Filter hifiasm contigs by minimum length.

Used to reduce the number of contigs passed to downstream steps
(e.g., BLAST, MITOS2, MitoFinder).
"""

import logging
from Bio import SeqIO


def filter_contigs_by_length(in_fasta, out_fasta, min_len):
    """
    Filters contigs from a FASTA file based on minimum length.

    Args:
        in_fasta (str): input FASTA file
        out_fasta (str): output FASTA file
        min_len (int): minimum contig length

    Returns:
        str: output FASTA filename
    """

    kept = 0
    total = 0

    with open(out_fasta, "w") as out_handle:
        for record in SeqIO.parse(in_fasta, "fasta"):
            total += 1
            if len(record.seq) >= min_len:
                SeqIO.write(record, out_handle, "fasta")
                kept += 1

    logging.info(f"[filter_hifiasm] Kept {kept}/{total} contigs >= {min_len} bp")

    if kept == 0:
        logging.warning("[filter_hifiasm] No contigs passed the filter!")

    return out_fasta
