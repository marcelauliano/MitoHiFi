from Bio import SeqIO
import os 
import logging
import sys

def get_num_seqs(in_fasta):
    """Gets the number of sequences in a FASTA file.
     
    Args:
        in_fasta (str): input FASTA file 
    
    Returns: 
        int: number of sequences in FASTA file
    """

    c = 0
    for rec in SeqIO.parse(in_fasta, "fasta"):
        c += 1
    return c

def get_ref_tRNA():
    """Defines the reference tRNA to be used for rotating contigs.   
    
    Returns: 
        str: reference tRNA
    """
    
    FORMAT='%(asctime)s [%(levelname)s] %(message)s'
    logging.basicConfig(level=logging.DEBUG, stream=sys.stdout,
                        format=FORMAT, datefmt='%Y-%m-%d %H:%M:%S')
    
    tRNAs = {}
    for curr_file in os.listdir('.'):
        if curr_file.endswith('.trnas'):
            curr_seq_tRNAs = set()
            with open(curr_file, "r") as infile:
                for line in infile:
                    tRNA = line.split("\t")[0][:4] #skip tRNA number if it exists
                    # adds tRNA to shared tRNA dict if it does not exist yet
                    if tRNA not in tRNAs:
                        tRNAs[tRNA] = 1
                        curr_seq_tRNAs.add(tRNA)
                    # it only adds one to tRNA count if it has not been counted
                    # in the local .trna file yet
                    elif tRNA not in curr_seq_tRNAs:
                        tRNAs[tRNA] += 1
                        curr_seq_tRNAs.add(tRNA)
    
    #logging.info(f"tRNAs: {tRNAs}") # debug
    #return tRNAs

    # if any contig has a tRNA-Phe, use it as the reference gene for rotation
    if 'trnF' in tRNAs:
        reference_tRNA = 'trnF'
    else:
        reference_tRNA = max(tRNAs, key=tRNAs.get)
    return reference_tRNA
