"""This script iterates over a multifasta sequence file and returns all proteins that 
contains stop codons in the middle of their sequences.

"""

import re
import sys
from Bio import SeqIO

def find_frameshifts_mitos(in_proteins_fasta):
    """Takes a multifasta protein file and returns information on proteins that contain frameshifts. """
    
    frameshift_genes = []
    for record in SeqIO.parse(in_proteins_fasta, "fasta"):
        seq = str(record.seq)
        description = record.description
        gene_name = description.split(";")[-1].strip()
        print(f"{gene_name}")
        if re.search("[A-Z]\*[A-Z]", seq):
            print("Gene {} contains frameshift".format(gene_name))
            frameshift_genes.append(gene_name)
    return frameshift_genes 

if __name__ == "__main__":
    in_fasta = sys.argv[1]
    print(find_frameshifts_mitos(in_fasta))
