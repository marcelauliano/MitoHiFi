#Credit: Marcela Uliano-Silva
from Bio import SeqIO
import sys
out1 = {}
for sequence in SeqIO.parse(sys.argv[1], "fasta"):
    id = sequence.id
    length = len(sequence)
    out1[id]=length
sort_out1 = sorted(out1.items(), key=lambda x: x[1], reverse=True)
largest=sort_out1[0]
print(largest[0])
