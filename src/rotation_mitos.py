#!/usr/bin/env python

from argparse import ArgumentParser
import os
from Bio import SeqIO
from Bio.Seq import Seq
import subprocess
from collections import deque
import itertools
import logging
import re

def get_trna_pos(path):
    """Gets the position for each tRNA in an input file

    Args:
        path (str): input file to be processed

    Returns:
        dict: Positions for each tRNA gene found
    """
    
    trnas = {}
    with open(path, "r") as f:
        for line in f:
            feat_type = line.split()[2]
            if feat_type == "tRNA":
                tRNA_start = line.split()[3]
                tRNA_strand = line.split()[6]
                feat_info = line.split()[-1].strip()
                tRNA_name = re.search("gene_id=.*", feat_info).group(0).split("=")[-1]
                trnas[tRNA_name] = [tRNA_start, tRNA_strand]

    if len(trnas) > 0:
        return trnas
    
    return None

def get_phe_pos(path):
    """Gets the position of the tRNA-Phe gene

    Args:
        path (str): path of input genbank file

    Returns:
        tuple: Start position of tRNA-Phe, Strand
    """

    with open(path) as f:
        for _, record in enumerate(SeqIO.parse(path, "genbank")):
            for fea in record.features:
                if fea.type == 'gene' and 'product' in fea.qualifiers and \
                    'tRNA-Phe' in fea.qualifiers['product'] :
                    return fea.location.start, fea.location.strand
    return None, None    
                
def make_rc(path, rc_path):
    """Creates a reverse complement.
     
    Args:
        path (str): input FASTA file 
        rc_path (str): new FASTA file created (reverse complement of input) 
    
    Returns: 
        None
    """
    record = SeqIO.read(path, "fasta")
    rc_record = record.reverse_complement(id=record.id + "_rc")
    with open(rc_path, 'w') as f:
        SeqIO.write(rc_record, f, 'fasta')

def annotate(workdir, path, ref_gb, contig_id, o_code, max_contig_size, threads):
    """Annotate reverse complemented genome.
    """

    cur = os.path.abspath(os.getcwd())
    out_id = contig_id + "_RC.annotation"
    if workdir:
        os.chdir(workdir)
    command = "mitofinder --max-contig-size {} -j {} -a {} -r {} -o {} -p {}".format(max_contig_size, out_id, path, ref_gb, o_code, threads) 
    mitofinder_cmd = ["mitofinder", "--new-genes", "--max-contig-size", str(max_contig_size),
                    "-j", out_id, "-a", path, "-r", ref_gb, "-o", o_code, "-p", str(threads),
                    "--circular-size", "8000"] 
    subprocess.run(mitofinder_cmd, stderr=subprocess.DEVNULL, stdout=subprocess.DEVNULL)
    os.chdir(cur)
    return os.path.join(workdir, out_id, \
                        out_id + '_MitoFinder_mitfi_Final_Results', \
                        out_id + '_mtDNA_contig.gb')

def rotate(genome, start, contig_id):
    """Creates new genome file (suffix .mitogenome.rotated.fa) after rotating it.

    Args:
        genome (str): input genome file
        start (int): position at which rotate the input genome
        contig_id (str): identifier representing the original contig 
    
    Returns:
        None
    """
    logging.info("Started rotation_mitos.rotate()") #debug
    logging.info(f"mitogenome_fasta: {genome} | start: {start} | contig_id: {contig_id}") #debug
    record = SeqIO.read(genome, "fasta")
    logging.info("record = SeqIO.read (...) [OK]")
    shifted_record = record[start-1:] + record[:start-1]
    #shifted_record = record[start:] + record[:start]
    logging.info("shifted_record = (...) [OK]")
    shifted_record.id = f"{record.id}_rotated"
    
    name = os.path.join(os.path.dirname(genome), contig_id + '.mitogenome.rotated.fa')
    logging.info(f"name: {name}") # debug
    with open(name, 'w') as f:
        SeqIO.write(shifted_record, f, 'fasta')

def rotate_annotation(mitogenome_annotation, mitogenome_fasta, start, contig_id):
    """Takes a GFF file and the reference gene position and rotates the annotation to have the reference gene at position 1."""
    
    logging.info(f"Started rotate_annotation function. mitogenome_annotation: {mitogenome_annotation} | mitogenome_fasta: {mitogenome_fasta} | start: {start} | contig_id: {contig_id}") #debug
    for record in SeqIO.parse(mitogenome_fasta, "fasta"):
        mitogenome_len = len(record.seq)
    
    features = [] 
    with open(mitogenome_annotation, "r") as f:
       for line in f:
            feature_fields = line.strip().split("\t")
            feature_start = int(feature_fields[3])
            feature_end = int(feature_fields[4])
            
            
            if feature_start == start:
                feature_fields[3] = 1
                feature_fields[4] = feature_end - start + 1
            elif feature_start > start:
                feature_fields[3] = feature_start - start + 1 
                feature_fields[4] = feature_end - start + 1
            elif feature_start < start and feature_end < start:
                feature_fields[3] = mitogenome_len - start + feature_start + 1
                feature_fields[4] = mitogenome_len - start + feature_end + 1
            elif feature_start < start and feature_end > start: # here we have a feature that overlaps the reference
                print("warning: feature being ignore because it overlaps with reference gene and starts before it:")
                print(feature_fields)
                continue
            
#            # if feature start equal to reference gene start
#            if feature_start == start:
#                feature_fields[3] = 1
#            # if feature starts after reference gene
#            elif feature_start > start:
#                feature_fields[3] = feature_start - start
#            # if feature starts before reference gene
#            elif feature_start < start: 
#                feature_fields[3] = mitogenome_len - start + feature_start
#
#            # if feature start equal to reference gene start
#            if feature_end == start:
#                feature_fields[4] = 1
#            # if feature starts after reference gene
#            elif feature_end > start:
#                feature_fields[4] = feature_end - start
#            # if feature starts before reference gene
#            elif feature_end < start: 
#                feature_fields[4] = mitogenome_len - start + feature_end

            features.append(feature_fields)
    logging.info(f"features: {features}") # debug 
    mitogenome_rotated_annotation = mitogenome_annotation.replace(".gff", ".rotated.gff")
    logging.info(f"mitogenome_rotated_annotation: {mitogenome_rotated_annotation}") # debug
    if os.path.exists(mitogenome_rotated_annotation):
        print("deleting existing output file")
        os.remove(mitogenome_rotated_annotation)
    
    logging.info("sorting features...") # debug
    features.sort(key=lambda s: s[3])
    print(f"mito length: {mitogenome_len}")
    for feat in features:
        with open(mitogenome_rotated_annotation, "a") as f:
            f.write("\t".join(str(x) for x in feat) + "\n")

    logging.info("rotate_annotation function done") # debug

if __name__ == '__main__':
    parser = ArgumentParser()
    parser.add_argument('--gff', help='path to annotation file')
    parser.add_argument('--start', help='start of reference gene', type=int)
    parser.add_argument('--mito', help='mitogenome fasta')
    parser.add_argument('--contig', help='contig ID')
    args = parser.parse_args()
    
    rotate_annotation(args.gff, args.mito, args.start, args.contig) 
