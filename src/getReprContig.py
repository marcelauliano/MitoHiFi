"""This script chooses between all potential mito contigs the representative one that is 
going to be considered the final mitogenome. 

The get_repr_contig() is the function originally referred to by the main mitohifi.py script.  
Over the process of choosing the final mitogenome, the get_circularization_info() function is 
called to check if the contig was circularized and the get_repr_contig_info() function calculates
some other metrics for all potential contigs and does the actual choice of the final mitogenome,
retriving its ID and its related stats.

The 
"""

import logging
import subprocess
import numpy
import warnings
from Bio import SeqIO
import sys #remove after debugging

def get_circularization_info(seq_id):
    """Retrieves information if contig was circularized
    
    Args:
        seq_id (str): identifier of the target sequence (contig)

    Returns: 
        bool: returns True if contig was circularized and False otherwise
    """

    with open(f"{seq_id}.circularisationCheck.txt", "r") as f:
        for line in f:
            if line == "(False, -1, -1)":
                return False
            else:
                return True

def get_repr_contig_info(cdhit_clstr_file, rel_mito_len, rel_mito_num_genes, rel_mito_perc=0.35, debug=False):
    
    def get_frameshift_info(seq_id):
        """Retrieves information of frameshifts from *.individual.stats file

        Args:
            seq_id (str): identifier of the target sequence (contig)

        Returns:
            str: information about frameshifts present on target sequence genes
        """

        with open(f"{seq_id}.individual.stats", "r") as f: 
            for line in f:
                frameshift_info = line.split("\t")[1]
        return frameshift_info
    
    def get_number_of_genes_diff(seq_id, rel_mito_num_genes):
        """Retrieves the number of genes from *.individual.stats file

        Args:
            seq_id (str): identifier of the target sequence (contig)
            rel_mito_num_genes (int): number of genes in related mito
        
        Returns:
            int: number of genes present on target sequence
        """
        with open(f"{seq_id}.individual.stats", "r") as f:
            for line in f:
                num_genes = int(line.split("\t")[4])
                num_genes_diff = abs(rel_mito_num_genes - num_genes)
            return num_genes_diff

    FORMAT='[%(asctime)s %(levelname)s] %(message)s'
    if debug: # If in debug mode
        logging.basicConfig(level=logging.DEBUG, stream=sys.stdout,
                            format=FORMAT, datefmt='%Y-%m-%d %H:%M:%S')

    # Set maximum size targeted for the representative contig 
    # based on the size of the related mitogenome
    rel_mito_len = int(rel_mito_len)
    rel_mito_upper_lim = (1+rel_mito_perc)*rel_mito_len
    logging.debug(f"Maximum size target for repr contig: {rel_mito_upper_lim}")

    # create `seqs` list of lists to hold information (ID, length, frameshifts, cluster)
    # for each sequence (contig)
    seqs = []
    with open(cdhit_clstr_file, "r") as f:
        for line in f:
            if line[0] == ">":
                curr_cluster = line.strip().replace(">", "") # save current cluster ID
            else:
                # retrieves the ID from sequences that were reverse complemented (rc)
                if "_rc_" in line: 
                    seq_id = line.split()[2].replace(">","").split("_rc_rotated")[0]
                # retrieves the ID from sequences that were **not** reverse complemented
                else:
                    seq_id = line.split()[2].replace(">","").split("_rotated")[0]
                # retrieves sequence lenght and frameshifts information
                seq_len = int(line.split()[1].replace("nt,", ""))
                seq_frameshifts = get_frameshift_info(seq_id)
                seq_circ = get_circularization_info(seq_id)
                seq_num_genes_diff = get_number_of_genes_diff(seq_id, rel_mito_num_genes)
                # appends sequence information (as a list) to the `seqs` list
                seqs.append([seq_id, seq_len, seq_frameshifts, curr_cluster, seq_circ, seq_num_genes_diff])

    # Sorts sequences based on the difference between
    # the number of genes of contigs and the related mito
    # (with contigs with lower differences placed first)
    seqs.sort(key=lambda x: x[5])
    logging.debug(f"Sorted seqs list: {seqs}")

    # Now we'll define the representative contig
    repr_contig = ""
    # First we search if there's any contig for which we have both i) a size less than or equal 
    # to an upper size limit defined based on the related mito length; ii) no frameshifts
    # found; and iii) that was detected as circularized by `circularizationCheck.py` and then
    # had the duplicated end removed
    for seq in seqs:
        if seq[1] <= rel_mito_upper_lim and seq[2] == "No frameshift found" and seq[4] == True:
            repr_contig, repr_contig_cluster = seq[0], seq[3]
            break
    
    # If the first condition is not met, then we search for a contig whose size is less than or
    # equal to the upper size limit and it was circularized        
    if not repr_contig:
        for seq in seqs:
            if seq[1] <= rel_mito_upper_lim and seq[4] == True:
                repr_contig, repr_contig_cluster = seq[0], seq[3]
                warnings.warn("Warning: representative contig contains frameshifts", stacklevel=2)
                break

    # If the second condition is not met, then we search for a contig whose size is less than or
    # equal to the upper size limit and it has no frameshifts        
    if not repr_contig:
        for seq in seqs:
            if seq[1] <= rel_mito_upper_lim and seq[2] == "No frameshift found":
                repr_contig, repr_contig_cluster = seq[0], seq[3]
                warnings.warn("Warning: representative contig wasn't circularized", stacklevel=2)
                break

    # If the third condition is not met, then we search for a contig which has no frameshifts and
    # it was circularized        
    if not repr_contig:
        for seq in seqs:
            if seq[2] == "No frameshift found" and seq[4] == True:
                repr_contig, repr_contig_cluster = seq[0], seq[3]
                warnings.warn("Warning: representative contig may be too large", stacklevel=2)
                break

    # If the fourth condition is not met, then we search for a contig whose size is less than or
    # equal to the upper size            
    if not repr_contig:
        for seq in seqs:
            if seq[1] <= rel_mito_upper_lim:
                repr_contig, repr_contig_cluster = seq[0], seq[3]
                warnings.warn(f"Warning: representative wasn't circularized and it has frameshifts", stacklevel=2)
                break

    # If the fifth condition is not met, then we search for a contig which was circularized            
    if not repr_contig:
        for seq in seqs:
            if seq[4] == True:
                repr_contig, repr_contig_cluster = seq[0], seq[3]
                warnings.warn(f"Warning: representative has frameshifts and it may be too large", stacklevel=2)
                break

    # If the sixtieth condition is not met, then we search for a contig which has no frameshifts            
    if not repr_contig:
        for seq in seqs:
            if seq[2] == "No frameshift found":
                repr_contig, repr_contig_cluster = seq[0], seq[3]
                warnings.warn(f"Warning: representative wasn't circularized and it may be too large", stacklevel=2)
                break

    # If none condition is met, we return the contig with lowest difference in the number of genes regarding the related mito            
    if not repr_contig:
        repr_contig, repr_contig_cluster = seqs[0][0], seqs[0][3]
        warnings.warn("Warning: representative contig contains frameshifts, wasn't circularized and it may be too large", stacklevel=2)

    return (repr_contig, repr_contig_cluster)   

def get_repr_contig(contigs_fasta, rel_mito_len, rel_mito_num_genes, threads="1", debug=False):
    """Gets representative contig from a multifasta file

    Args:
        contigs_fasta (str): file containing all sequences 
        threads (str): number of threads to be used when running CDHIT

    Returns:
        tuple: Representative contig ID, CDHIT cluster where the representative contig came from
    """

    c_threshold = "0.8"
    wordsize = "4"
    cdhit_out = "cdhit.out"
    cdhit_out_clstr = cdhit_out + ".clstr"
    cdhit_cmd = ["cd-hit-est", "-i", contigs_fasta, "-d", "0", "-c", c_threshold, "-n", wordsize, "-o", cdhit_out, "-T", str(threads), "-M", "0"]
    subprocess.run(cdhit_cmd, shell=False, stderr=subprocess.DEVNULL, stdout=subprocess.DEVNULL)

    repr_contig_id, repr_contig_cluster = get_repr_contig_info(cdhit_out_clstr, rel_mito_len, rel_mito_num_genes, debug=debug)
    return (repr_contig_id, repr_contig_cluster)

if __name__ == "__main__":
    repr_contig_id, repr_contig_cluster = get_repr_contig(sys.argv[1], sys.argv[2])
    print("Representative contig is {} that belongs to {}".format(repr_contig_id, repr_contig_cluster))
