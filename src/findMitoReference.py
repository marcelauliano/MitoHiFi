#!/usr/bin/env python

"""
This script finds mitogenomes from related species in public databases. 

License:
    Copyright 2022 Ksenia Krasheninnikova
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

"""

from argparse import ArgumentParser
import os
from Bio import Entrez
from Bio import SeqIO
from io import StringIO

def get_lineage(species):
    handle = Entrez.esearch(db="taxonomy", term=species, idtype="acc")
    record = Entrez.read(handle)
    if not len(record['IdList']) :
        raise Exception('No such species in NCBI!') 
    if len(record['IdList']) > 1 :
        raise Exception('More than one appropriate entires in NCBI!') 
    handle = Entrez.efetch(db="Taxonomy", id=record['IdList'][0], retmode="xml")
    records = Entrez.read(handle)
    handle.close()
    for e in records[0]["Lineage"].split(' ')[::-1]:
        yield e

def find_full_mito(group, outfolder, length_threshold, considered, org_type='mitochondrion', n=1):
    term = "(\""+ group +"\"[Organism] AND complete " +\
                "genome[All Fields]) AND "+ org_type +"[filter]]"
    handle = Entrez.esearch(db="nucleotide",term=term, idtype="acc")
    record = Entrez.read(handle)
    if record['IdList'] :
        for ncbi_code in record['IdList']:
            handle = Entrez.efetch(db="nucleotide", id=ncbi_code, rettype="gb", retmode="text")
            record_ = handle.read()
            handle.close()            
            for seqrecord in  SeqIO.parse(StringIO(record_), "gb") :
                if seqrecord.id in considered:
                    continue
                considered.add(seqrecord.id)
                if len(seqrecord) > length_threshold:
                    if len(seqrecord.features) < 10:
                        print('Not enough features in gb file! skipping ' + seqrecord.id)
                        continue
                    with open(os.path.join(outfolder, ncbi_code+'.gb') , "w") as out:
                        out.write(record_)
                    handle = Entrez.efetch(db="nucleotide", id=ncbi_code, rettype="fasta", \
                                                                        retmode="text")
                    record_ = handle.read()
                    handle.close()
                    with open(os.path.join(outfolder, ncbi_code+'.fasta') , "w") as out:
                        out.write(record_)
                    print("output is written to " + os.path.join(outfolder,ncbi_code) + ".[gb,fasta]")
                    n -= 1
                    if n == 0:
                        return considered, n
    return considered, n

if __name__ == '__main__':
    parser = ArgumentParser()
    parser.add_argument('--species', required=True, help='latin name')
    parser.add_argument('--email', nargs='?', default="")
    parser.add_argument('--outfolder', nargs='?', default="")
    parser.add_argument('--type', choices=['mitochondrion','chloroplast'], default='mitochondrion', \
                                  help='Specify the type of organelle')
    parser.add_argument('--min_length', nargs='?', type=int, default=0, \
                                  help='Minimal appropriate length')
    parser.add_argument('-n', nargs='?',  type=int, default=1, \
                                  help='Number of genomes to report. Reported in order of identification')
    parser.add_argument('--ncbi-api-key', nargs='?', help='Set NCBI API key')
    args = parser.parse_args()    
    n = args.n
    if n < 1:
        print('Number of genomes to report must be at least 1 (default)')
        exit(1)
    Entrez.email = args.email
    if args.ncbi_api_key:
        Entrez.api_key = args.ncbi_api_key
    if not os.path.isdir(args.outfolder):
        os.mkdir(args.outfolder)
    print('Looking for ' + args.type +' for ' + args.species)
    considered = set()
    for g in [args.species] + list(get_lineage(args.species)):
        if n > 0:
            print('Looking for an appropriate organelle among ' + g)
            considered, n = find_full_mito(g, args.outfolder, args.min_length, considered, args.type, n)
    if n == args.n:            
        print("No appropriate mitogenome found")

