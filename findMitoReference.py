#!/usr/bin/env python

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

def query_same_species(species, outfolder, min_length):
    term = "\""+ species +"\" [Organism] AND "  +\
            "mitochondrion[All Fields]"
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

def find_full_mito(group, outfolder, length_threshold):
    term = "(\""+ group +"\"[Organism] AND complete " +\
                "genome[All Fields]) AND mitochondrion[filter]]"
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
    parser.add_argument('--email', required=True)
    parser.add_argument('--outfolder', nargs='?', default="closestMitoReference", \
                        help='Folder to house the downloaded reference. Default: closestMitoReference')
    parser.add_argument('-s', action='store_true', help='search for an exact species')
    parser.add_argument('--min_length', nargs='?', type=int  , default=0, \
                                  help='minimal appropriate length')
    args = parser.parse_args()    
    Entrez.email = args.email
    if not os.path.isdir(args.outfolder):
        os.mkdir(args.outfolder)
    if args.s:
        ret = query_same_species(args.species, args.outfolder, args.min_length)
        if ret == 1:
            print('No appropriate mitogenome found')
        exit(ret)
    for g in get_lineage(args.species):
        if find_full_mito(g, args.outfolder, args.min_length) == 0:
           exit(0) 
    print("No appropriate mitogenome found")
    exit(1)
