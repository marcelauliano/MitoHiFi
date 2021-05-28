#!/usr/bin/env python


'''
    Copyright 2021 Joāo Ferreira Nunes
    This script is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.
    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.
    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
    '''




import subprocess
import numpy
from Bio import SeqIO
import sys #to remove after debugging


def get_largest_cluster(cdhit_clstr_file):
    #clusters = {}
    largest_cluster_len = 0
    largest_cluster = ""
    curr_sequences = []
    curr_cluster = []
    with open(cdhit_clstr_file, "r") as f:
        for line in f:
            if line[0] == ">":
                if len(curr_sequences) > largest_cluster_len:
                    largest_cluster = curr_cluster
                    largest_cluster_len = len(curr_sequences)
                    largest_cluster_seqs = curr_sequences
                curr_cluster = line.strip().replace(">","")
                curr_sequences = []
                #clusters[cluster_id] = 0
            else:
                curr_sequences.append(line.strip())
        # catch the last cluster
        if len(curr_sequences) > largest_cluster_len:
            largest_cluster = curr_cluster
            largest_cluster_len = len(curr_sequences)
            largest_cluster_seqs = curr_sequences
    print("[x] Current sequences: {}".format(curr_sequences))
    for sequence in largest_cluster_seqs:
        if sequence[-1] == "*":
            representative_seq = sequence

    return (largest_cluster, representative_seq)

def get_repr_contig(contigs_fasta, threads="1"):
    c_threshold = "0.8"
    wordsize = "4"
    cdhit_out = "cdhit.out"
    cdhit_out_clstr = cdhit_out + ".clstr"
    cdhit_cmd = ["cd-hit-est", "-i", contigs_fasta, "-d", "0", "-c", c_threshold, "-n", wordsize, "-o", cdhit_out, "-T", str(threads)]
    subprocess.run(cdhit_cmd, shell=False)

    repr_contig_cluster, repr_contig_info = get_largest_cluster(cdhit_out_clstr)
    if "_rc_" in repr_contig_info:
        repr_contig_id = repr_contig_info.split()[2].replace(">","").split("_rc_rotated")[0]
    else:
        repr_contig_id = repr_contig_info.split()[2].replace(">","").split("_rotated")[0]
    return (repr_contig_id, repr_contig_cluster)

if __name__ == "__main__":
    repr_contig_id, repr_contig_cluster = get_repr_contig(sys.argv[1], sys.argv[2])
    print("Representative contig is {} that belongs to {}".format(repr_contig_id, repr_contig_cluster))
