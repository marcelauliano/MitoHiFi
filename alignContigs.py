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


import os
import subprocess
from Bio import SeqIO
from Bio.Seq import Seq

def concatenate_contigs(contigs_list, out_file="all_mitogenomes.rotated.fa"):
    # create a list of BioSeqs that contains all contigs
    contigs = []
    for contig in contigs_list:
        with open(contig) as handle:
            for record in SeqIO.parse(handle, "fasta"):
                contigs.append(record)
    # write concatenated FASTA file
    with open(out_file, "w") as output_handle:
        SeqIO.write(contigs, output_handle, "fasta")

    return out_file

def mafft_align(multifasta_file, threads='1', out_file="all_mitogenomes.rotated.aligned.fa", clustal_format=False):
    if clustal_format:
        print("MAFFT output will be in clustal format")
        mafft = "mafft --clustalout --thread {} {} > {}".format(threads, multifasta_file, out_file)
    else:
        mafft = "mafft --thread {} {} > {}".format(threads, multifasta_file, out_file)
    print("MAFFT alignment will be called with:\n" + mafft + "\n")
    subprocess.run(mafft, shell=True)

if __name__ == '__main__':
    contigs_files = []
    for curr_file in os.listdir('.'):
        if curr_file.endswith('mitogenome.rotated.fa'):
            contigs_files.append(curr_file)
    concatenate_contigs(contigs_files)
