from Bio.SeqUtils import seq1
import sys
from getGenesList import get_genes_list

def get_genes_counts(genes_list):
    """
    Takes a list of genes and returns a dictionary containing the name of the 
    gene and the number of occurrences.
    """
    
    genes_count = {}
    for gene in genes_list:
        if gene not in genes_count:
            genes_count[gene] = 1
        else:
            genes_count[gene] += 1

    return genes_count

def get_clean_gene(in_gene):
    """
    Takes a raw gene name and returns its clean format.
    """
    if in_gene.lower().startswith("trn"):
        gene_clean = in_gene.replace("tRNA-", "trn")
        amino_code = gene_clean[3:]
        if len(amino_code) == 3:
            gene_clean = "trn" + seq1(gene_clean[3:])
        else:
            gene_clean = "trn" + amino_code
    else:
            gene_clean = in_gene.replace("-", "").lower()

def compare_genes_dicts(genes1, genes2):
    """
    Takes two lists of genes and return genes shared and specific.
    """

    genes1_counts = get_genes_counts(genes1)
    genes2_counts = get_genes_counts(genes2)
    
    genes1_counts_clean = {}
    for gene1, n in genes1_counts.items():
        if gene1.lower().startswith("trn"):
            gene1_clean = gene1.replace("tRNA-", "trn")
            print(f"gene1_clean: {gene1_clean}")
            amino_code = gene1_clean[3:]
            if len(amino_code) == 3:
                gene1_clean = "trn" + seq1(gene1_clean[3:])
            else:
                gene1_clean = "trn" + amino_code
        else:
            gene1_clean = gene1.replace("-", "").lower()
        genes1_counts_clean[gene1_clean] = n
    
    genes2_counts_clean = {}
    for gene2, n in genes2_counts.items():
        if gene2.lower().startswith("trn"):
            gene2_clean = gene2.replace("tRNA-", "trn")
            amino_code = gene2_clean[3:]
            if len(amino_code) == 3:
                gene2_clean ="trn" + seq1(gene2_clean[3:])
            else:
                gene2_clean = "trn" + amino_code
        else:
            gene2_clean = gene2.replace("-", "").lower()
        genes2_counts_clean[gene2_clean] = n
    
    #genes1_counts_clean = {k.replace("-", "").replace("tRNA", "trn").lower(): v for k, v in genes1_counts.items()}
    #genes2_counts_clean = {k.replace("-", "").replace("tRNA", "trn").lower(): v for k, v in genes2_counts.items()}
    print("### genes1_counts:")
    print(genes1_counts)
    print("### genes1_counts_clean:")
    print(genes1_counts_clean)
    print("### genes2_counts")
    print(genes2_counts)
    print("### genes2_counts_clean:")
    print(genes2_counts_clean)
    
    shared_genes = {}
    genes1_specific_genes = {}
    genes2_specific_genes = {}

    for gene1 in genes1_counts:
        gene1_clean = gene1.replace("-", "").lower()
        if gene1_clean in genes2_counts_clean:
            shared_genes[gene1] = [genes1_counts[gene1], genes2_counts_clean[gene1_clean]]
        else:
            genes1_specific_genes[gene1] = [genes1_counts[gene1], 0]
    
    for gene2 in genes2_counts:
        gene2_clean = gene2.replace("-", "").lower()
        if gene2_clean not in genes1_counts_clean:
            genes2_specific_genes[gene2] = [0, genes2_counts[gene2]]

    print(f"shared genes: {shared_genes}")
    print(f"specific to file1: {genes1_specific_genes}")
    print(f"specific to file2: {genes2_specific_genes}")

if __name__ == "__main__":
    genes_list1 = get_genes_list(sys.argv[1], sys.argv[2])
    genes_list2 = get_genes_list(sys.argv[3], sys.argv[4])
    #print(genes_list1)
    #print(genes_list2)
    compare_genes_dicts(genes_list1, genes_list2)
    #print(get_genes_counts(genes_list))
