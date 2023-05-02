#!/usr/bin/python
#Version: 1.0
#Author: Alex Schomaker - alexschomaker@ufrj.br
#LAMPADA - IBQM - UFRJ

'''
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
'''

from Bio import SeqIO, SearchIO
import argparse
from subprocess import Popen
import shlex, sys, os

def circularizationCheck(resultFile, contig_id, circularSize=220, circularOffSet=40):
	'''
	Check, with blast, if there is a match between the start and the end of a sequence.
	Returns a tuple with (True, start, end) or False, accordingly.
	'''
	refSeq = SeqIO.read(resultFile, "fasta")
	sizeOfSeq = len(refSeq)

	try:
		command = "makeblastdb -in " + resultFile + " -dbtype nucl" #need to formatdb refseq first
		args = shlex.split(command)
		formatDB = Popen(args, stdout=open(os.devnull, 'wb'))
		formatDB.wait()
	except:
		print ('')
		print ("formatDB during circularization check failed...")
		print ('')
		return (False,-1,-1)
		
	with open(f"{contig_id}.circularization_check.blast.tsv",'w') as blastResultFile:
		command = "blastn -task blastn -db " + resultFile + " -query " + resultFile + " -outfmt 6" #call BLAST with TSV output
		args = shlex.split(command)
		blastAll = Popen(args, stdout=blastResultFile)
		blastAll.wait()

	blastparse = SearchIO.parse(f'{contig_id}.circularization_check.blast.tsv', 'blast-tab') #get all queries

	'''
	Let's loop through all blast results and see if there is a circularization.
	Do it by looking at all HSPs in the parse and see if there is an alignment of the ending of the sequence 
    with the start of that same sequence. It should have a considerable size, you don't want to say it circularized
	if only a couple of bases matched.
	Returns True or False, x_coordinate, y_coordinate
	x coordinate = starting point of circularization match
	y coordinate = ending point of circularization match
	'''
	for qresult in blastparse: #in each query...
		for hsp in qresult.hsps: #loop through all HSPs looking for a circularization (perceived as a hsp with start somewhat close to the query finish)
			if (hsp.query_range[0] >= 0 and hsp.query_range[0] <= circularOffSet) and (hsp.hit_range[0] >= sizeOfSeq - hsp.aln_span - circularOffSet and hsp.hit_range[0] <= sizeOfSeq + circularOffSet) and hsp.aln_span >= circularSize and hsp.aln_span < sizeOfSeq * 0.90:
				if hsp.hit_range[0] < hsp.query_range[0]:
					return (True,hsp.hit_range[0],hsp.hit_range[1]) #it seems to have circularized, return True
				else:
					return (True,hsp.query_range[0],hsp.query_range[1])

	#no circularization was observed in the for loop, so we exited it, just return false
	return (False,-1,-1)

def get_circo_mito(contig_id, circular_size, circular_offset):
    """It takes a contig ID and circularizes it, returning the circularization points.

    Args:
    	contig_id (str): ID from contig to be circularized
    	circular_size (int): size to consider when checking for circularization
    	circular_offset (int): offset from start and finish to consider when looking for circularization

    Returns:
    	list: circularization history for the input contig
    """
    def cut_coords(contig_fasta, circularization_position, fasta_out):
        #print(f"cut_coords function called!\ncontig_fasta: {contig_fasta}; circularization_position: {circularization_position}; fasta_out: {fasta_out}") # for debugging
        record=SeqIO.read(contig_fasta, "fasta")
        id = record.id
        get= record[circularization_position:]
        with open(fasta_out, "w") as f:
            f.write(get.format('fasta'))

    # find circularization point
    circularization_history = []
    contig_fasta = "".join([contig_id, ".mito.fa"])
    mitogenome_fasta_out = "".join([contig_id, ".mitogenome.fa"])
    circularization_checks = "".join([contig_id, ".circularisationCheck.txt"])

    circularization_info = circularizationCheck(resultFile=contig_fasta, circularSize=circular_size, circularOffSet=circular_offset, contig_id=contig_id)

    # writes circularization information to '[contig_id].circularisationCheck.txt' file
    with open(circularization_checks, "w") as f:
        f.write(str(circularization_info))
        circularization_history.append(str(circularization_info))
    
    is_circularizable = circularization_info[0]
    circularization_position = int(circularization_info[2])

    # if the contig is not circularizable, then create the "mitogenome.fasta" directly from it
    if not is_circularizable:
        record=SeqIO.read(contig_fasta, "fasta")
        id = record.id
        with open(mitogenome_fasta_out, "w") as f:  
            f.write(record.format('fasta'))        
    else:
        # if the contig is circularizable, then run circularization iteratively until the "mitogenome.fasta"
        # is no longer circularizable    
        while is_circularizable:
            cut_coords(contig_fasta, circularization_position, mitogenome_fasta_out)
            contig_fasta = mitogenome_fasta_out
            circularization_info = circularizationCheck(resultFile=contig_fasta, circularSize=circular_size, circularOffSet=circular_offset, contig_id=contig_id)
            is_circularizable = circularization_info[0]
            circularization_position = int(circularization_info[2])
            with open(circularization_checks, "a") as f:
                f.write("\n" + str(circularization_info))
                circularization_history.append(str(circularization_info))
    
    return circularization_history

if __name__ == '__main__':
	parser = argparse.ArgumentParser(description='Find circularization points')
	parser.add_argument('--result-file', type=str, help='Input FASTA file to be processed')
	parser.add_argument('--circular-size', type=int, default=220,
						help='Size to consider when checking for circularization')
	parser.add_argument('--circular-offset', type=int, default=40,
						help='Offset from start and finish to consider when looking for circularization')

	args = parser.parse_args()
	
	if sys.argv[1] == '-h' or sys.argv[1] == '--help':
		print ('Usage: fasta_file')
	else:
		print(circularizationCheck(args.result_file, args.circular_size, args.circular_offset))
