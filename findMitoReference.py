#!/usr/bin/env python

'''
    Copyright 2021 Ksenia Krasheninnikova
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

def query_same_species(species, outfolder, min_length, org_type='mitochondrion'):
    term = "\""+ species +"\" [Organism] AND "  +\
            org_type+"[All Fields]"
    handle = Entrez.esearch(db="nucleotide",term=term, idtype="acc")
    record = Entrez.read(handle)
    rs = {}
    max_length = min_length
    max_id = "" 
    if record['IdList'] :
        for ncbi_code in record['IdList']:
            handle = Entrez.efetch(db="nucleotide", id=ncbi_code, rettype="gb", retmode="text")
            record_ = handle.read()
            handle.close()
            for seqrecord in  SeqIO.parse(StringIO(record_), "gb") :
                if len(record_) > max_length:
                    if len(seqrecord.features) < 10:
                        print('Not enough features in gb file! skipping..')
                        continue
                    rs[ncbi_code] = record_
                    max_length = len(seqrecord)
                    max_id = ncbi_code
    if max_id:
        ncbi_code = max_id
        with open(os.path.join(outfolder, ncbi_code+'.gb') , "w") as out:
            out.write(rs[ncbi_code].format("fasta"))
        handle = Entrez.efetch(db="nucleotide", id=ncbi_code, rettype="fasta", \
                                                                        retmode="text")
        record_ = handle.read()
        handle.close()
        with open(os.path.join(outfolder, ncbi_code+'.fasta') , "w") as out:
            out.write(record_.format("genbank"))
        print("output is written to " + os.path.join(outfolder,ncbi_code) + ".[gb,fasta]")
        return 0
    return 1

def find_full_mito(group, outfolder, length_threshold, org_type='mitochondrion'):
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
                if len(seqrecord) > length_threshold:
                    if len(seqrecord.features) < 10:
                        print('Not enough features in gb file! skipping..')
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
                    return 0
    return 1

if __name__ == '__main__':
    parser = ArgumentParser()
    parser.add_argument('--species', required=True, help='latin name')
    parser.add_argument('--email', nargs='?', default="")
    parser.add_argument('--outfolder', nargs='?', default="")
    parser.add_argument('-t', choices=['mitochondrion','chloroplast'], default='mitochondrion', \
                                  help='specify the type of organelle')
    parser.add_argument('--min_length', nargs='?', type=int  , default=0, \
                                  help='minimal appropriate length')
    args = parser.parse_args()    
    Entrez.email = args.email
    if not os.path.isdir(args.outfolder):
        os.mkdir(args.outfolder)
    print('Looking for '+args.t+' for '+args.species)
    ret = query_same_species(args.species, args.outfolder, args.min_length, args.t)
    if ret == 0:
        exit(0)
    print('Mito for the same species is not found')
    print('Looking among close species')
    if ret == 1:
        for g in get_lineage(args.species):
            if find_full_mito(g, args.outfolder, args.min_length, args.t) == 0:
                exit(0) 
    print("No appropriate mitogenome found")
    exit(ret)
