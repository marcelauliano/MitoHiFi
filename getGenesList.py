import sys
from   Bio import SeqIO


def get_genes_list(in_annotation, format="genbank"):
    """
    Takes a GFF/Genbank file and returns a list of annotated genes.
    """
    
    genes = []
    
    if format == "genbank":
        for record in SeqIO.parse(in_annotation, format):
            for feat in record.features:
                if feat.type in ["tRNA", "rRNA"]:
                    genes.append(feat.qualifiers['product'][0])
                    #genes.append([feat.qualifiers['product'], feat.location])
                elif feat.type == "CDS":
                    genes.append(feat.qualifiers['gene'][0])
                    #genes.append([feat.qualifiers['gene'], feat.location])
    
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
