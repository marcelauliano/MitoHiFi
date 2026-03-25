"""This script returns a list of annotated genes. 

It recognizes annotations in either genbank or GFF format.

"""

import sys
from   Bio import SeqIO


def get_genes_list(in_annotation, format="genbank"):
    """
    Takes a GFF/Genbank file and returns a list of annotated genes.
    """
    
    genes = []
    
    if format == "genbank":
        name_replacement = {
            # long names
            "NADH dehydrogenase subunit 6": "ND6", "NADH dehydrogenase subunit 4L": "ND4L", "12S ribosomal RNA": "s-rRNA", "NADH dehydrogenase subunit 1": "ND1", "ATP synthase F0 subunit 6": "ATP6", "NADH dehydrogenase subunit 2": "ND2", "cytochrome b": "CYTB", "cytochrome c oxidase subunit III": "COIII", "NADH dehydrogenase subunit 4": "ND4", "cytochrome c oxidase subunit I": "COI", "cytochrome c oxidase subunit II": "COII", "16S ribosomal RNA": "l-rRNA", "NADH dehydrogenase subunit 3": "ND3", "NADH dehydrogenase subunit 5": "ND5",

           # gene symbols (common variants)
            "cox1": "COI", "COX1": "COI", "cox2": "COII", "cox3": "COIII", "cob": "CYTB", "cytb": "CYTB", "nad1": "ND1", "nad2": "ND2", "nad3": "ND3", "nad4": "ND4", "nad4l": "ND4L", "nad5": "ND5", "nad6": "ND6", "rrnS": "s-rRNA", "rrnL": "l-rRNA"
        }

        for record in SeqIO.parse(in_annotation, format):
            for feat in record.features:
                if feat.type in ["tRNA", "rRNA", "CDS"]:
                
                    # --- SAFE extraction of gene name ---
                    gene_name = None

                    if 'product' in feat.qualifiers:
                        gene_name = feat.qualifiers['product'][0]
                    elif 'gene' in feat.qualifiers:
                        gene_name = feat.qualifiers['gene'][0]
                    elif 'locus_tag' in feat.qualifiers:
                        gene_name = feat.qualifiers['locus_tag'][0]
                    else:
                        continue  # skip if nothing usable

                    # --- apply replacement safely ---
                    gene_name = name_replacement.get(gene_name, gene_name)

                    genes.append(gene_name)
    
    elif format == "gff":
        with open(in_annotation, "r") as f:
            for line in f:
                feat_type = line.split("\t")[2]
                if feat_type in ["ncRNA_gene", "gene"]:
                    gene_name = line.split("\t")[-1].split(";")[-1].replace("gene_id=", "").strip()
                    genes.append(gene_name)
    else:
        sys.exit("Incompatible format, please specify either genbank or gff")
        
    return genes 

if __name__ == "__main__":
    in_annotation = sys.argv[1]
    format = sys.argv[2]
    genes = get_genes_list(in_annotation, format)
    print(genes)
