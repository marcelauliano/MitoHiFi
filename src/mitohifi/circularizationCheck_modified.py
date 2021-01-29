#!/usr/bin/env python
# Version: 1.0
# Author: Alex Schomaker - alexschomaker@ufrj.br
# LAMPADA - IBQM - UFRJ

"""
Copyright (c) 2014 Alex Schomaker Bastos - LAMPADA/UFRJ

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
"""
import os
import shlex
import sys
from subprocess import Popen

from Bio import SearchIO, SeqIO


def circularizationCheck(resultFile, circularSize, circularOffSet):
    """
    Check, with blast, if there is a match between the start and the end of a sequence.
    Returns a tuple with (True, start, end) or False, accordingly.
    """
    refSeq = SeqIO.read(resultFile, "fasta")
    sizeOfSeq = len(refSeq)
    id = refSeq.id
    print("Formatting database for blast (circularization check)...")
    try:
        # if blastFolder == 'installed':
        #     command = "formatdb -in " + resultFile + " -p F" # need to formatdb refseq first
        # else:
        command = f"makeblastdb -in {resultFile} -dbtype nucl"  # need to formatdb refseq first
        args = shlex.split(command)
        formatDB = Popen(args, stdout=open(os.devnull, "wb"))
        formatDB.wait()
    except Exception:
        print("\nformatDB during circularization check failed...\n")
        return (False, -1, -1)

    print("Running blast of input against itself to check for circularization...")

    with open("circularization_check.blast.xml", "w") as blastResultFile:
        # if blastFolder == 'installed':
        # command = "blastall -p blastn -d " + resultFile + " -i " + resultFile + " -m 7" #call BLAST with XML output
        # else:
        command = f"blastn -task blastn -db {resultFile} -query {resultFile} -outfmt 5"  # call BLAST with XML output
        args = shlex.split(command)
        blastAll = Popen(args, stdout=blastResultFile)
        blastAll.wait()

    blastparse = SearchIO.parse(
        "circularization_check.blast.xml", "blast-xml"
    )  # get all queries

    """
    Let's loop through all blast results and see if there is a circularization.
    Do it by looking at all HSPs in the parse and see if there is an alignment of the ending of the sequence
    with the start of that same sequence. It should have a considerable size, you don't want to say it circularized
    if only a couple of bases matched.
    Returns True or False, x_coordinate, y_coordinate
    x coordinate = starting point of circularization match
    y coordinate = ending point of circularization match
    """
    for qresult in blastparse:  # in each query...
        # loop through all HSPs looking for a circularization
        # (perceived as a hsp with start somewhat close to the query finish)
        for hsp in qresult.hsps:
            if (
                hsp.query_range[0] >= 0
                and hsp.query_range[0] <= circularOffSet
                and hsp.hit_range[0] >= sizeOfSeq - hsp.aln_span - circularOffSet
                and hsp.hit_range[0] <= sizeOfSeq + circularOffSet
                and hsp.aln_span >= circularSize
                and hsp.aln_span < sizeOfSeq * 0.90
            ):
                if hsp.hit_range[0] < hsp.query_range[0]:
                    return (
                        str(id),
                        True,
                        hsp.hit_range[0],
                        hsp.hit_range[1],
                    )  # it seems to have circularized, return True
                else:
                    return (str(id), True, hsp.query_range[0], hsp.query_range[1])

    # no circularization was observed in the for loop, so we exited it, just return false
    return (str(id), False, -1, -1)


if __name__ == "__main__":
    # dirpath = os.getcwd()
    # outfile = os.path.join(dirpath, 'circularizationCheck.out')
    output = open("circularisationCheck.txt", "w")

    if sys.argv[1] == "-h" or sys.argv[1] == "--help":
        print("Usage: fasta_file")
    else:
        run = circularizationCheck(
            sys.argv[1],
            220,
            40,
        )
        output.write("%s, %s, %i, %i" % run)
    output.close()
