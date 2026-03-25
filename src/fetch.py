"""This script allows the calculation of statistics over MitoHiFi execution.

The get_num_seqs() function is used to calculate the number of sequences in a FASTA file. 
The get_ref_tRNA() function allows the definition of the reference tRNA to be used to rotate 
the potential mito contigs.

"""

from Bio import SeqIO
import os 

#def get_num_seqs(in_fasta):
#    """Gets the number of sequences in a FASTA file.
#     
#    Args:
#        in_fasta (str): input FASTA file 
#    
#    Returns: 
#        int: number of sequences in FASTA file
#    """
#
#    c = 0
#    for rec in SeqIO.parse(in_fasta, "fasta"):
#        c += 1
#    return c

def get_num_seqs(file_path):
    """
    Count number of sequences in FASTA or FASTQ file.
    """

    count = 0
    with open(file_path) as f:
        first_char = f.read(1)
        f.seek(0)

        # FASTA
        if first_char == ">":
            for line in f:
                if line.startswith(">"):
                    count += 1

        # FASTQ
        elif first_char == "@":
            for i, line in enumerate(f):
                pass
            count = (i + 1) // 4

        else:
            raise ValueError("Unknown sequence format")

    return count

def get_ref_tRNA():
    """Defines the reference tRNA to be used for rotating contigs.   
    
    Returns: 
        str: reference tRNA
    """

    tRNAs = {}
    for curr_file in os.listdir('.'):
        if curr_file.endswith('.trnas'):
            with open(curr_file, "r") as infile:
                for line in infile:
                    tRNA = line.split("\t")[0]
                    if tRNA not in tRNAs:
                        tRNAs[tRNA] = 1
                    else:
                        tRNAs[tRNA] += 1
    # if any contig has a tRNA-Phe, use it as the reference gene for rotation
    if 'tRNA-Phe' in tRNAs:
        reference_tRNA = 'tRNA-Phe'
    else:
        reference_tRNA = max(tRNAs, key=tRNAs.get)
    return reference_tRNA

