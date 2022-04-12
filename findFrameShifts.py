from Bio import SeqIO
import re
import sys

def find_frameshifts(in_gb):
    frameshift_genes = set()
    for record in SeqIO.parse(in_gb, "gb"):
        for feature in record.features:
            if feature.type == 'CDS':
                gene_name = feature.qualifiers['gene'][0]
                gene_amino = feature.qualifiers['translation'][0]
                if re.search("[A-Z]\*[A-Z]", gene_amino):
                    print("Gene {} contains frameshift".format(gene_name))
                    frameshift_genes.add(gene_name)
    return frameshift_genes 

def get_gb_stats(in_gb):
    for record in SeqIO.parse(in_gb, "gb"):
        gb_len = len(record.seq)
        num_genes = 0
        for feature in record.features:
            if feature.type == 'gene':
                num_genes += 1
    return (str(gb_len), str(num_genes))

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
