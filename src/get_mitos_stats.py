from Bio import SeqIO
import sys

def get_mitos_stats(in_gff, in_fasta):
    """Takes in a GFF file produced by MITOS and returns annotation statistics."""
    
    for record in SeqIO.parse(in_fasta, "fasta"):
        seq_len = len(record.seq)

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
    
    num_genes = len(protein_coding_genes) + len(ncRNA_genes)
    return (seq_len, num_genes)
    return (protein_coding_genes, ncRNA_genes)

if __name__ == "__main__":
    in_gff = sys.argv[1]
    in_fasta = sys.argv[2]
    get_mitos_stats(in_gff, in_fasta)
