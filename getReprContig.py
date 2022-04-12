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

def get_repr_contig_info(cdhit_clstr_file, rel_mito_len, rel_mito_perc=0.10, debug=False):
    
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
                # appends sequence information (as a list) to the `seqs` list
                seqs.append([seq_id, seq_len, seq_frameshifts, curr_cluster, seq_circ])

    # Sorts (descending order) sequences based on second item of list,
    # which represents the sequence lengths
    seqs.sort(key=lambda x: x[1], reverse=True)
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
                warnings.warn("Warning: representative contig contains frameshifts")
                break

    # If the second condition is not met, then we search for a contig whose size is less than or
    # equal to the upper size limit and it has no frameshifts        
    if not repr_contig:
        for seq in seqs:
            if seq[1] <= rel_mito_upper_lim and seq[2] == "No frameshift found":
                repr_contig, repr_contig_cluster = seq[0], seq[3]
                warnings.warn("Warning: representative contig wasn't circularized")
                break

    # If the third condition is not met, then we search for a contig which has no frameshifts and
    # it was circularized        
    if not repr_contig:
        for seq in seqs:
            if seq[2] == "No frameshift found" and seq[4] == True:
                repr_contig, repr_contig_cluster = seq[0], seq[3]
                warnings.warn("Warning: representative contig may be too large")
                break

    # If the fourth condition is not met, then we search for a contig whose size is less than or
    # equal to the upper size            
    if not repr_contig:
        for seq in seqs:
            if seq[1] <= rel_mito_upper_lim:
                repr_contig, repr_contig_cluster = seq[0], seq[3]
                warnings.warn(f"Warning: representative wasn't circularized and it has frameshifts")
                break

    # If the fifth condition is not met, then we search for a contig which was circularized            
    if not repr_contig:
        for seq in seqs:
            if seq[4] == True:
                repr_contig, repr_contig_cluster = seq[0], seq[3]
                warnings.warn(f"Warning: representative has frameshifts and it may be too large")
                break

    # If the sixtieth condition is not met, then we search for a contig which has no frameshifts            
    if not repr_contig:
        for seq in seqs:
            if seq[2] == "No frameshift found":
                repr_contig, repr_contig_cluster = seq[0], seq[3]
                warnings.warn(f"Warning: representative wasn't circularized and it may be too large")
                break

    # If none condition is met, we return the smallest contig available            
    if not repr_contig:
        repr_contig, repr_contig_cluster = seqs[-1][0], seqs[-1][3]
        warnings.warn("Warning: representative contig contains frameshifts, wasn't circularized and it may be too large")

    return (repr_contig, repr_contig_cluster)   

def get_repr_contig(contigs_fasta, rel_mito_len, threads="1", debug=False):
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

    repr_contig_id, repr_contig_cluster = get_repr_contig_info(cdhit_out_clstr, rel_mito_len, debug=debug)
    return (repr_contig_id, repr_contig_cluster)

if __name__ == "__main__":
    repr_contig_id, repr_contig_cluster = get_repr_contig(sys.argv[1], sys.argv[2])
    print("Representative contig is {} that belongs to {}".format(repr_contig_id, repr_contig_cluster))
