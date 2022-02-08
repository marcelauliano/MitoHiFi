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


def get_mito_length(mito_file):
    """
    Returns the length of a mitogenome

    Keyword arguments:
    mito_file -- the input mitogenome FASTA file
    """

    from Bio import SeqIO
    num_sequences = len(list(SeqIO.parse(mito_file, "fasta")))
    if num_sequences > 1:
        raise Exception("failed because mitogenome file has more than one sequence. {} sequences were found".format(str(num_sequences)))
    for index, record in enumerate(SeqIO.parse(mito_file, "fasta")):
        mito_len = len(record.seq)
    if not mito_len:
        raise Exception("Could not find the length of the mitogenome. Check if it is in proper FASTA format")
    return mito_len

def main():
    import argparse 

    parser = argparse.ArgumentParser(description="Get length of mitogenome")
    parser.add_argument('mito_file', type=str, help="Mitogenome FASTA file")

    args = parser.parse_args()
    
    mito_len = get_mito_length(args.mito_file)

    print("Length of mitogenome {} is {} bp".format(args.mito_file, mito_len))

if __name__ == "__main__":
    main()
