import sys
from Bio import SeqIO

def fix_header(in_gb, out_gb):
    
    record = SeqIO.read(in_gb, "genbank")
    record.id = record.description
    record.name = record.description
    with open(out_gb, "w") as f:
        SeqIO.write(record, f, "genbank")

if __name__ == "__main__":
    in_gb = sys.argv[1]
    out_gb = sys.argv[2]
    fix_header(in_gb, out_gb)
