"""This script fixes headers from contigs given as input to MitoHiFi in contigs mode (`-c`).

The fix removes/changes suffixes from the contigs IDs in case they conflict with suffixes that
are given by MitoHiFi while the pipeline runs.

"""

import sys 


def fix_headers(fasta_in, fasta_out):
    with open (fasta_in, "r") as f:
        with open(fasta_out, "a") as f2:
            for line in f:
                if ">" in line: 
                    original_id = line[1:].split()[0]
                    if any(x in original_id for x in ["_rc", "_rotated"]):
                        fixed_id = original_id.replace("_rc", ".rc").replace("_rotated", ".rotated")
                        f2.write(f">{fixed_id}\n")
                        print(f'Header substitution: {original_id} replaced by {fixed_id}')
                    else:
                        f2.write(line.strip() + "\n")
                else:
                    f2.write(line.strip() + "\n")

if __name__ == "__main__":
    fix_headers(sys.argv[1], sys.argv[2])
