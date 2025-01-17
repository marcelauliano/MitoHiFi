"""This script iterates over a multifasta sequence file and returns all proteins that 
contains stop codons in the middle of their sequences. 

The find_frameshifts() function is used for the (default) annotation using `MitoFinder`,
while find_frameshifts_mitos() is used when `MITOS` was selected as the annotation tool.

"""

from Bio import SeqIO
import re
import sys
import logging 

def find_frameshifts(in_gb):
    frameshift_genes = set()
    for record in SeqIO.parse(in_gb, "gb"):
        for feature in record.features:
            if feature.type == 'CDS':
                gene_name = feature.qualifiers['gene'][0]
                gene_amino = feature.qualifiers['translation'][0]
                if re.search(r"[A-Z]\*[A-Z]", gene_amino):
                    print("Gene {} contains frameshift".format(gene_name))
                    frameshift_genes.add(gene_name)
    return frameshift_genes 

def find_frameshifts_mitos(in_proteins_fasta):
    """Takes a multifasta protein file and returns information on proteins that contain frameshifts. """

    frameshift_genes = []
    for record in SeqIO.parse(in_proteins_fasta, "fasta"):
        seq = str(record.seq)
        description = record.description
        gene_name = description.split(";")[-1].strip()
        if re.search(r"[A-Z]\*[A-Z]", seq):
            print("Gene {} contains frameshift".format(gene_name))
            frameshift_genes.append(gene_name)
    return frameshift_genes 

def get_gb_stats(in_gb):
    for record in SeqIO.parse(in_gb, "gb"):
        gb_len = len(record.seq)
        num_genes = 0
        for feature in record.features:
            if feature.type == 'gene':
                num_genes += 1
    return (str(gb_len), str(num_genes))


def get_mitos_stats(in_gff, in_fasta):
    """Takes in a GFF file produced by MITOS and returns annotation statistics."""
    
    logging.info(f"started get_mitos_stats()") # debug
    logging.info(f"in_gff: {in_gff} | in_fasta: {in_fasta}") # debug


    for record in SeqIO.parse(in_fasta, "fasta"):
        seq_len = len(record.seq)
    
    logging.info(f"seq_len: {seq_len}") # debug

    protein_coding_genes = []
    ncRNA_genes = []

    with open(in_gff, "r") as f:
        for line in f.readlines():
            feat_type = line.split()[2]
            if feat_type == "gene":
                feat_name = line.split()[-1].split(";")[-1].replace("gene_id=", "")
                protein_coding_genes.append(feat_name) # debug
            elif feat_type == "ncRNA_gene":
                feat_name = line.split()[-1].split(";")[-1].replace("gene_id=", "")
                ncRNA_genes.append(feat_name)
    
    logging.info(f"protein_coding_genes: {protein_coding_genes}") # debug
    logging.info(f"ncRNA_genes: {ncRNA_genes}") # debug

    num_genes = len(protein_coding_genes) + len(ncRNA_genes)
    return (seq_len, num_genes)

def main():
    frameshifts_found = find_frameshifts(sys.argv[1])
    if frameshifts_found:
        print("Frameshifts were found in genes:", frameshifts_found)
    else:
        print("No frameshift found")
    print("retrieving stats...")
    gb_len, num_genes = get_gb_stats(sys.argv[1])
    print(f"gb_len: {gb_len} bp; num_genes: {num_genes}")

if __name__ == "__main__":
    main()
