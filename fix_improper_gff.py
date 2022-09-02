import logging
import os
import sys


def fix_improper_gff(in_gff):
    """Takes a GFF file and fixes all features that have start position greater than end position."""
    
    was_fixed = False
    fixed_filename = in_gff.replace(".gff", ".fixed.gff")
    
    if os.path.exists(fixed_filename):
        logging.info(f"Deleting existing fixed file: {fixed_filename}")
        os.remove(fixed_filename)

    
    with open(in_gff, "r") as f:
        for line in f:
            start = int(line.split("\t")[3])
            end = int(line.split("\t")[4])
            if start > end:
                feature_description = line.split("\t")[-1]
                logging.warning(f"Feature with start position > end position is improper GFF notation. Ignoring it: {feature_description}")
                was_fixed = True
            else:
                with open(fixed_filename, "a") as f2:
                    f2.write(line)

    # if `in_gff` was fixed, rename it by adding a `.old` suffix 
    # and give its original name to its fixed version 
    if was_fixed:
        os.rename(in_gff, in_gff.replace(".gff", ".old.gff"))
        os.rename(fixed_filename, in_gff)
    else:
        os.remove(fixed_filename)

if __name__ == "__main__":
    in_gff = sys.argv[1]
    fix_improper_gff(in_gff)
