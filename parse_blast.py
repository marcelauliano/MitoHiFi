#! /usr/bin/env python

'''
    Copyright 2020 Marcela Uliano-Silva
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

import logging            
import pandas as pd
import sys

def parse_blast(query_perc=50, min_query_perc=80, max_query_len=5):
    logging.info("Filtering thresholds applied:")
    logging.info(f"Minimum query percentage = {query_perc}")
    logging.info(f"Minimum query length = {min_query_perc}% subject length")
    logging.info(f"Maximum query length = {max_query_len} times subject length")
    my_names = ["qseqid", "sseqid", "pident", "alilength" , "mismatch", "gapopen", "qstart", "qend", 
                "sstart", "send", "evalue" , "bitscore", "leng_query", "s_length",]

    blast_cov = pd.read_csv("contigs.blastn", 
                            sep="\t", names = my_names)

    #Get the percentage of the query in the blast aligment
    blast_cov['alilength']*100 / (blast_cov['leng_query'])
    blast_cov['%q_in_match'] = blast_cov['alilength']*100 / (blast_cov['leng_query'])

    #sum percentages of query sequence in blast match based on column id
    a= blast_cov.groupby('qseqid')['%q_in_match'].sum().to_frame().rename(columns={'qseqid':'%q_in_match'}).reset_index()

    #get size of query and subject and drop duplicates
    seqsizes = blast_cov[['qseqid', 'leng_query', 's_length']].drop_duplicates(subset='qseqid')

    #merge 'a' and 'seqsizes' dataframes by 'qseqid'
    result = pd.merge(a, seqsizes, on='qseqid')
    result
    # Now let's filter the blast matches
    # if the lenght of the query is `max_query_len` times (or greater than)
    # the size of the subject (close-related mitogenome), drop it. 
    # (As its likely the match belongs to a NUMT)
    n_times = (result['s_length'] * max_query_len)
    result1 = result[(result['leng_query'] < n_times)].sort_values(by='%q_in_match', ascending=False)

    # if the lenght of the query is 80% smaller than the length of the subject, drop it. It's unlikely you will have a complete mitogenome.
    slen=result1['s_length']
    result1['perc'] = result1["leng_query"]*100/(result1["s_length"])
    ac=result1[result1['perc'] > min_query_perc].sort_values(by='%q_in_match')

    # if the % of the query in the blast match is smaller than 70%, drop it
    #ac[(ac['%q_in_match'] > 60)].sort_values(by='%q_in_match', ascending=False)
    ac[(ac['%q_in_match'] > 0)].sort_values(by='%q_in_match', ascending=False).to_csv("parsed_blast_all.txt", index=False, sep="\t")
    ac[(ac['%q_in_match'] > query_perc)].sort_values(by='%q_in_match', ascending=False).to_csv("parsed_blast.txt", index=False, sep="\t")

def get_contigs_ids(blast_output):
    """
    args:
    blast_line is a blast output, that can be either parsed_blast.txt or parsed_blast_all.txt
    
    returns:
    The ID from each contig from blast_output, i.e., the BLAST queries
    """
    contigs_ids = set()

    num_lines = sum(1 for line in open(blast_output, "r") if line.strip()) # counts the number of lines in the input file
    if num_lines >= 2:
        with open(blast_output, "r") as f:
            next(f) # skips blast header
            for line in f:
                contigs_ids.add(line.split()[0])    

    return contigs_ids

def main():
    parse_blast()

if __name__ == '__main__':
    main()

