#!/usr/bin/env python

from argparse import ArgumentParser
import os
from Bio import SeqIO
from Bio.Seq import Seq
import subprocess
from collections import deque
import itertools

def get_phe_pos(path):
    with open(path) as f:
        for _, record in enumerate(SeqIO.parse(path, "genbank")):
            for fea in record.features:
                if fea.type == 'gene' and 'product' in fea.qualifiers and \
                    'tRNA-Phe' in fea.qualifiers['product'] :
                    return fea.location.start, fea.location.strand
    return None, None    
                
def make_rc(path, rc_path):
    record = SeqIO.read(path, "fasta")
    rc_record = record.reverse_complement(id=record.id + "_rc")
    with open(rc_path, 'w') as f:
        SeqIO.write(rc_record, f, 'fasta')

def annotate(workdir, path, ref_gb, contig_id, o_code, max_contig_size, threads):
    cur = os.path.abspath(os.getcwd())
    out_id = contig_id + "_RC.annotation"
    if workdir:
        os.chdir(workdir)
    command = "mitofinder --max-contig-size {} -j {} -a {} -r {} -o {} -p {}".format(max_contig_size, out_id, path, ref_gb, o_code, threads) 
    #command = ' '.join(['mitofinder -j', out_id, '-a', path, '-r', ref_gb, '-o', o_code])
    print("Running mitofinder for RC: ", command)
    subprocess.call(command, shell=True)
    os.chdir(cur)
    return os.path.join(workdir, out_id, \
                        out_id + '_MitoFinder_mitfi_Final_Results', \
                        out_id + '_mtDNA_contig.gb')

def rotate(genome, start, contig_id):
    record = SeqIO.read(genome, "fasta")
    d = deque(record.seq)
    d.rotate(-(start+1))
    record.seq = Seq(''.join(d))
    record.id = record.id + '_rotated'
    name = os.path.join(os.path.dirname(genome), contig_id + '.mitogenome.rotated.fa')
    with open(name, 'w') as f:
        SeqIO.write(record, f, 'fasta')

if __name__ == '__main__':
    parser = ArgumentParser()
    parser.add_argument('--gb', help='path to genbank annotation file')
    parser.add_argument('--ref-gb', help='path to genbank of the reference genome')
    parser.add_argument('--mito', help='mitogenome fasta')
    args = parser.parse_args()
    print('Extracting initial position of tRNA-Phe...')
    start, strand = get_phe_pos(args.gb)
    genome = args.mito
    if start == None:
        raise Exception('tRNA-Phe is not present in file '+args.gb)
    new_gb = None
    if strand == -1:
        print('Making reverse complement of mitogenome...')
        rc = os.path.join(os.path.dirname(args.mito), 'RC_mito.fa')
        make_rc(args.mito, rc)
        print('Annotating reverse complement sequence...')
        new_gb = annotate(os.path.dirname(args.mito), os.path.abspath(rc), 
                             os.path.abspath(args.ref_gb))
        print('Extracting initial position of tRNA-Phe...')
        start, strand = get_phe_pos(new_gb)
        genome = rc
    print('Final rotation...')
    rotate(genome, start)
    print(' '.join(['Rotated to tRNA-Phe genome is at ', \
            os.path.join(os.path.dirname(genome), 'mitogenome.rotated.fa')]))
    if new_gb:
        print('Mitogenome annotation is at ', new_gb)
    print('Done')
        

