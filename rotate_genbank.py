from Bio import SeqIO
from Bio.SeqFeature import SeqFeature, FeatureLocation, ExactPosition
import sys

def get_feat_info(feat):

    feat_start = feat.location.start
     
    if "product" in feat.qualifiers:
        feat_name = feat.qualifiers["product"][0]
    elif "gene" in feat.qualifiers:
        feat_name = feat.qualifiers["gene"][0]
    else:
        print(f"Error. Couldn't find feature name for feature:\n{feat}.")
        exit(1)
    
    return (feat_start, feat_name)

def rotate_genbank(in_gbk, ref_gene, out_gbk):
    
    ref_start = None
    for seq_record in SeqIO.parse(in_gbk, "genbank"):
        complete_seq_len = len(seq_record.seq)
        for feat in seq_record.features:
            if feat.type == "gene":
                feat_start, feat_name = get_feat_info(feat)
                if feat_name == ref_gene:
                    ref_start = feat.location.start
    
    if ref_start is None:
        print(f"Error. Reference gene {ref_gene} not found")
        sys.exit(1)
    
    # read original record and shift (rotate) it using the 
    # start of the reference gene as the rotation point 
    # (i.e. sequence will now start at ref_start)
    record = SeqIO.read(in_gbk, "genbank")
    shifted_record = record[ref_start:] + record[:ref_start]
    
    #shifted_record.id = f"{record.id}_rotated".replace("_draft", "")
    #shifted_record.name = f"{record.name}_rotated".replace("_draft", "")
    
    shifted_record.id = f"{record.id}_rotated"
    shifted_record.name = f"{record.name}_rotated"

    with open(out_gbk, "w") as f:
        SeqIO.write(shifted_record, f, "genbank")

    SeqIO.convert(out_gbk, "genbank", out_gbk.replace(".gb", ".fa"), "fasta")


if __name__ == "__main__":
        
    if sys.argv[1] == "-h":
        print("""python rotate_genbank.py <in_gbk> <ref_gene> <out_gbk>            

    in_gbk: input genbank file to be modified
    species_name: name of gene to be used as reference for rotation
    out_gbk: name of new, modified genbank file
                """)
        sys.exit(0)
    
    in_gbk = sys.argv[1]
    ref_gene = sys.argv[2]
    out_gbk = sys.argv[3]

    rotate_genbank(in_gbk, ref_gene, out_gbk)
