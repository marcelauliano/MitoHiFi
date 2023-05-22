"""This script creates the reverse complement of a sequence.

The reverse_complement() function is the default function to create the reverse complement.
The reverse_complement_mitos() function is the one that should be used when dealing with MITOS annotations.
The reverse_complement_annotation() function adjusts the annotation coordinates of features to match the ones from the reverse complemented sequence.

"""

import sys
from Bio import SeqIO

def reverse_complement(in_gb, out_gb):
    
    # open original record and convert it to 
    # its reverse complement (rc_record)
    record = SeqIO.read(in_gb, "genbank")
    rc_record = record.reverse_complement(id=record.id + "_rc")
    # needs to set molecule type as DNA because Biopython deletes
    # this info when parsing
    rc_record.annotations["molecule_type"] = "DNA"
    
    with open(out_gb, "w") as f:
        SeqIO.write(rc_record, f, "genbank")
    
    out_fasta = out_gb.replace(".gb", ".fasta")
    with open(out_fasta, "w") as f:
        SeqIO.write(rc_record, f, "fasta")

def reverse_complement_mitos(in_fasta, out_fasta):
    
    # open original record and convert it to 
    # its reverse complement (rc_record)
    record = SeqIO.read(in_fasta, "fasta")
    rc_record = record.reverse_complement(id=record.id + "_rc")
    
    with open(out_fasta, "w") as f:
        SeqIO.write(rc_record, f, "fasta")

def reverse_complement_annotation(mitogenome_annotation, mitogenome_fasta, mitogenome_annotation_rc):
    
    for record in SeqIO.parse(mitogenome_fasta, "fasta"):
        mitogenome_len = len(record.seq)

    with open(mitogenome_annotation, "r") as f:
        for line in f:
            fields = line.strip().split()
            new_start = mitogenome_len - int(fields[3]) + 1
            new_end = mitogenome_len - int(fields[4]) + 1
            old_strand = fields[6]
            
            # switches the signal
            if old_strand == "+":
                new_strand = "-"
            else:
                new_strand = "+"
            
            # sets first coord as start
            # and second coord as end
            if new_start > new_end:
                b = new_start
                new_start = new_end
                new_end = b

            fields[3] = str(new_start)
            fields[4] = str(new_end)
            fields[6] = new_strand
            with open(mitogenome_annotation_rc, "a") as f2:
                f2.write("\t".join(fields) + "\n")

if __name__ == "__main__":
    in_gb = sys.argv[1]
    out_gb = sys.argv[2]
    reverse_complement(in_gb, out_gb)
