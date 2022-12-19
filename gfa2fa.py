import os
import sys

def gfa2fa(gfa_input, fasta_output=""):
    
    if not fasta_output:
        fasta_output = gfa_input.replace(".gfa", ".fa")
    
    # delete out fasta file if it already exists
    # to override its content as opposed to appending
    if os.path.exists(fasta_output):
        os.remove(fasta_output)

    with open(gfa_input) as f1:
        for line in f1:
            if line.startswith("S"):
                seqID, sequence = line.split("\t")[1:3]
                with open(fasta_output, "a") as f2:
                    f2.write(f">{seqID}\n{sequence}\n")

if __name__ == "__main__":
    gfa_input = sys.argv[1]
    gfa2fa(gfa_input)
