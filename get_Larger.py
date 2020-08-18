from Bio import SeqIO
import sys
#output=open("/lustre/scratch116/vr/projects/vgp/user/mu2/mito-automation/test-pipe/fasta_to_circu", "w")
out1 = {}
for sequence in SeqIO.parse(sys.argv[1], "fasta"):
    id = sequence.id
    length = len(sequence)
    out1[id]=length
sort_out1 = sorted(out1.items(), key=lambda x: x[1], reverse=True)
largest=sort_out1[0]
print(largest[0])
