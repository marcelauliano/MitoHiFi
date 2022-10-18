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
        if gene.lower().startswith("trn") and gene[-1].isdigit():
            gene = gene[:-1]
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
        # if last character is digit, than it represents the copy number and we want to get rid of it
        if in_gene[-1].isdigit(): 
            in_gene = in_gene[:-1]
        
        gene_clean = in_gene.replace("tRNA-", "trn") # set trn as standard for tRNA genes
        # turn three letter coded aminoacids into one letter code
        amino_code = gene_clean[3:]
        if len(amino_code) == 3:
            gene_clean = "trn" + seq1(gene_clean[3:])
        else:
            gene_clean = "trn" + amino_code
    # NAD and ND are synonyms
    elif in_gene.lower().startswith("nad"):
        gene_clean = in_gene.lower().replace("nad", "nd")
    # CYTB and COB are synonyms
    elif in_gene.lower() == "cytb":
        gene_clean = "cob"
    # 16S and rrnL are synonyms
    elif in_gene.lower().startswith("16s"):
        gene_clean = "rrnl"
    # 12S and rrnS are synonyms
    elif in_gene.lower().startswith("12s"):
        gene_clean = "rrns"
    # for all other genes, remove hiphens and turn to lower case
    else:
        gene_clean = in_gene.replace("-", "").lower()

    return gene_clean

def compare_genes_dicts(genes1, genes2, alphabetically_sorted=False):
    """
    Takes two lists of genes and return genes shared and specific.
    """

    genes1_counts = get_genes_counts(genes1)
    genes2_counts = get_genes_counts(genes2)
    
    genes1_counts_clean = {}
    for gene1, n in genes1_counts.items():
        gene1_clean = get_clean_gene(gene1)
        genes1_counts_clean[gene1_clean] = n
    
    genes2_counts_clean = {}
    for gene2, n in genes2_counts.items():
        gene2_clean = get_clean_gene(gene2)
        genes2_counts_clean[gene2_clean] = n
    
    #print("### genes1_counts:")
    #print(genes1_counts)
    #print("### genes1_counts_clean:")
    #print(genes1_counts_clean)
    #print("### genes2_counts")
    #print(genes2_counts)
    #print("### genes2_counts_clean:")
    #print(genes2_counts_clean)
    
    shared_genes = {}
    genes1_specific_genes = {}
    genes2_specific_genes = {}

    for gene1 in genes1_counts:
        gene1_clean = get_clean_gene(gene1)
        if gene1_clean in genes2_counts_clean:
            shared_genes[gene1] = [genes1_counts[gene1], genes2_counts_clean[gene1_clean]]
        else:
            genes1_specific_genes[gene1] = [genes1_counts[gene1], 0]
    
    for gene2 in genes2_counts:
        gene2_clean = get_clean_gene(gene2)
        if gene2_clean not in genes1_counts_clean:
            genes2_specific_genes[gene2] = [0, genes2_counts[gene2]]
    
    if alphabetically_sorted:
        shared_genes_sorted = dict(sorted(shared_genes.items()))
        genes1_specific_genes_sorted = dict(sorted(genes1_specific_genes.items()))
        genes2_specific_genes_sorted = dict(sorted(genes2_specific_genes.items()))
        return (shared_genes_sorted, genes1_specific_genes_sorted, genes2_specific_genes_sorted)
    else:
        return (shared_genes, genes1_specific_genes, genes2_specific_genes)

if __name__ == "__main__":
    if sys.argv[1] in ["-h", "-help", "--help"]:
        print("""Run usage:
        python compare <annotationFile1> <file1Format>[genbank|gff] <annotationFile2> <file2Format>[genbank|gff] <returnGenesInAlphabeticalOrder>[True|False]\n""")
        sys.exit(0)

    annotation1 = sys.argv[1]
    annotation1_format = sys.argv[2]
    annotation2 = sys.argv[3]
    annotation2_format = sys.argv[4]
    alphabetically_sorted_string = sys.argv[5]
    if alphabetically_sorted_string == "True":
        alphabetically_sorted = True
    elif alphabetically_sorted_string == "False":
        alphabetically_sorted = False
    else:
        print("alphabetically_sorted argument must be either True or False")
        sys.exit(1)
    genes_list1 = get_genes_list(annotation1, annotation1_format)
    genes_list2 = get_genes_list(annotation2, annotation2_format)
    shared, annotation1_specific, annotation2_specific = compare_genes_dicts(genes_list1, genes_list2, alphabetically_sorted)
    print(f"Shared genes:\n{shared}\n\n")
    print(f"Specific to {annotation1}:\n{annotation1_specific}\n\n")
    print(f"Specific to {annotation2}:\n{annotation2_specific}\n\n")
