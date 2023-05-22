"""This script circularizes, annotates and rotates mito contigs using MITOS for annotation. 

The process_contig_mitos() function does the circularization, i.e. removal of artifactual repeated 
sequences at both ends of the contig and the annotation, i.e. gene prediction.
The process_contig_02_mitos() function rotates the contig, given a reference gene that will be set as the
beginning of the sequence. This function also calculates statistics for the contig, which are saved
to a file named `{contig_id}.individual.stats`.

In the context of the mitohifi.py script, both functions are usually run in parallel for each
potential mito contig. The process_contig_02_mitos() function will only be called after process_contig_mitos() is run for all potential contigs.

"""

import re
import os
import concurrent.futures
import sys
import functools
import subprocess
import getMitoLength
import logging
import warnings 
import rotation_mitos
import shutil
from Bio import SeqIO 
from circularizationCheck import get_circo_mito
from getReprContig import get_circularization_info
from findFrameShifts import find_frameshifts_mitos, get_mitos_stats
import filterfasta
from fix_improper_gff import fix_improper_gff 
from rotate_genbank import rotate_genbank
from reverse_complement import reverse_complement, reverse_complement_mitos, reverse_complement_annotation

def process_contig_mitos(threads_per_contig, circular_size, circular_offset, contigs, max_contig_size, rel_gbk, gen_code, refseq_db, contig_id): 
    """Circularize and annotate a contig.
     
    Args:
        threads_per_contig (int): number of threads to be used
        circular_size (int): size to consider when checking for circularization
        circular_offset (int): offset from start and finish to consider when looking for circularization
        contigs (str): filename of contigs file (containing all contigs)
        max_contig_size (int): maximum contig size allowed
        rel_gbk (str): filename of related mito genbank file
        gen_code (str): species genetic code
        contig_id (str): target contig ID
        refseq_db (str): refseq database to be used

    Returns: 
        None
    """    
    
    logging.info(f"Working with contig {contig_id}") 
    # retrieves the FASTA files for each contig
    filterfasta.filterFasta(idList=[contig_id], inStream=contigs, outPath="".join([contig_id, ".mito.fa"]), log=False)
    # circularizes each contig and saves circularization history to a file
    logging.info(f"Started {contig_id} circularization")
    circularization_history = get_circo_mito(contig_id, circular_size, circular_offset)
    logging.info(f"{contig_id} circularization done. Circularization info saved on ./potential_contigs/{contig_id}/{contig_id}.circularisationCheck.txt")
    # annotates mitogenome(s) using mitofinder
    logging.info(f"Started {contig_id} (MITOS2) annotation")
    annotation_dir = f"{contig_id}.annotation"
    
    # removes annotation dir if it exists
    if os.path.isdir(annotation_dir):
        shutil.rmtree(annotation_dir, ignore_errors=True)
        logging.info(f"Deleted existing annotation directory {annotation_dir}")
    # creates a new annotation dir
    os.mkdir(annotation_dir)
    mitos2_cmd = ["python2", "runmitos.py", "--noplots", "-i", contig_id+".mitogenome.fa", "-c", gen_code,
                    "--refdir", "MITOS/data", "-r", refseq_db, "-o", annotation_dir]
    subprocess.run(mitos2_cmd) 
    #logging.info(f"{contig_id} annotation done. Annotation log saved on ./potential_contigs/{contig_id}/{contig_id}.annotation_MitoFinder.log")
    logging.info(f"{contig_id} annotation done.")
    
    # check if annotation file has been created
    mitogenome_gff = os.path.join(contig_id + ".annotation", "result.gff") 
    if not os.path.isfile(mitogenome_gff):
        warnings.warn("Contig "+ contig_id + " does not have an annotation file, an error has possibly occurred during MITOS2 execution")
        return    
    
    # check and remove if any feature has a start position greater than its end position
    fix_improper_gff(mitogenome_gff)
    
    #rotate the mitogenome
    trnas = rotation_mitos.get_trna_pos(mitogenome_gff)
    if trnas:
        with open(f"{contig_id}.trnas", "w") as outfile:
            for trna in trnas:
                outfile.write("\t".join([trna, trnas[trna][0], trnas[trna][1], "\n"]))
        return
    else:
        warnings.warn(f"No tRNA gene found in {mitogenome_gff}... Skipping contig {contig_id}")
        return

def process_contig_02_mitos(ref_tRNA, threads_per_contig, circular_size, circular_offset, contigs, max_contig_size, rel_gbk, gen_code, contig_id): 
    """Rotate a contig related to a reference tRNA gene and calculate contig statistics. 
     
    Args:
        ref_tRNA (str): tRNA gene to be used as reference for rotation (contig starts at reference tRNA)
        threads_per_contig (int): number of threads to be used
        circular_size (int): size to consider when checking for circularization
        circular_offset (int): offset from start and finish to consider when looking for circularization
        contigs (str): filename of contigs file (containing all contigs)
        max_contig_size (int): maximum contig size allowed
        rel_gbk (str): filename of related mito genbank file
        gen_code (str): species genetic code
        contig_id (str): target contig ID

    Returns: 
        None
    """    
    logging.info(f"Started process_contig_02_mitos function | tRNA_ref: {ref_tRNA}") # debug
    logging.info(f"Started {contig_id} rotation.")
    if not os.path.isfile(f"{contig_id}.trnas"):
        warnings.warn(f"Contig {contig_id} does not have annotated tRNAs, skipping it...")
        return
    
    logging.info(f"Opening {contig_id}.trnas") #debug 
    start = None
    with open(f"{contig_id}.trnas", "r") as infile:
        for line in infile:
            tRNA = line.strip().split("\t")[0][:4] # [:4] skips tRNA number of it exists 
            if tRNA == ref_tRNA:
                start = int(line.strip().split("\t")[1])
                strand = line.strip().split("\t")[2]
    
    if start == None:
        warnings.warn(f"Reference gene {ref_tRNA} is not present in contig {contig_id}. Skipping contig...")
        return

    logging.info(f"ref_tRNA: {ref_tRNA} | start: {start} | strand: {strand}") 
    mitogenome_annotation = os.path.join(contig_id + ".annotation", "result.gff")
    #mitogenome_gb = contig_id + ".mitogenome.gb"
    #shutil.copy(mitogenome_annotation, mitogenome_gb)
    # create copy of annotation file without the _draft modification
    # implemented by MitoFinder at ID and NAME sections
    
#    # comment out (no GB in MITOS mode)
#    record = SeqIO.read(mitogenome_annotation, "genbank")
#    record.id = record.description
#    record.name = record.description
#    with open(mitogenome_gb, "w") as f:
#        SeqIO.write(record, f, "genbank")

    mitogenome_fasta = contig_id + ".mitogenome.fa"
    if strand == "-":
        logging.info(f"{ref_tRNA} is at reverse complement of {contig_id}.mitogenome.fa")
        logging.info(f"For that reason we'll reverse complement {contig_id}.mitogenome.fa before the rotation")
        #mitogenome_gb_rc = mitogenome_gb.replace("")
        
        #genome_rc = contig_id + "_RC.mitogenome."
        #rc = os.path.join(os.path.dirname(genome), genome_rc)
        #rotation.make_rc(genome, rc)
        mitogenome_fasta_rc = re.sub(r".fa$", "_RC.fasta", mitogenome_fasta)
        reverse_complement_mitos(mitogenome_fasta, mitogenome_fasta_rc)
        mitogenome_annotation_rc = re.sub(r".gff$", "_RC.gff", mitogenome_annotation)
        reverse_complement_annotation(mitogenome_annotation, mitogenome_fasta, mitogenome_annotation_rc)
        #logging.info(f"Reverse complement generated: . Starting reverse complement annotation...")
        #mitogenome_gb = rotation.annotate(os.path.dirname(genome), os.path.abspath(genome_rc), os.path.abspath(rel_gbk), contig_id, gen_code, max_contig_size, str(threads_per_contig))
        logging.info(f"Reverse complementation for contig {contig_id} done")
        mitogenome_fasta = mitogenome_fasta_rc
        mitogenome_annotation = mitogenome_annotation_rc
        if not os.path.isfile(mitogenome_fasta_rc):
            warnings.warn("Contig "+ contig_id + " does not have a reverse complemented file, check log.")
            return
    
        # update reference gene start position after reverse complementation
        with open(mitogenome_annotation, "r") as f:
            for line in f:
                if line.split("\t")[2] == "tRNA":
                    tRNA_info = line.strip().split("\t")[-1] 
                    tRNA_ID = tRNA_info.split(";")[-1].split("=")[-1][:4]
                    if tRNA_ID == ref_tRNA:
                        start = int(line.split("\t")[3])
                        break

        logging.info("Finished reverse complementation.")
    #sys.exit(0)
    
    logging.info("Starting rotation_mitos.rotate function")
    rotation_mitos.rotate(mitogenome_fasta, start, contig_id)
    #mitogenome_rotated_fasta = f"{contig_id}.mitogenome.rotated.fa"
    #rotate_genbank(mitogenome_gb, ref_tRNA, mitogenome_rotated_gb)
    
    rotated_file = os.path.join(os.path.dirname(mitogenome_fasta), contig_id + '.mitogenome.rotated.fa')
    mitogenome_fasta = rotated_file 
    try:
        f = open(rotated_file)
    except FileNotFoundError:
        sys.exit(f"""Rotated file {mitogenome_rotated_gb} not found.
        An error may have occurred during the rotation step.""")
    finally:
        f.close()
    
    rotation_mitos.rotate_annotation(mitogenome_annotation, mitogenome_fasta, start, contig_id)
    mitogenome_rotated_annotation = mitogenome_annotation.replace(".gff", ".rotated.gff")
    logging.info(f"mitogenome_rotated_annotation: {mitogenome_rotated_annotation}") #debug
    logging.info(f"Rotation of {contig_id} done. Rotated is at {rotated_file}") 
    # check frameshifts in genes from contig and save findings to 
    # `{contig_id}.individual.stats` intermediate file
    logging.info(f"Started calculating mitocontig stats... {contig_id}") # debug
    mitogenome_faa = mitogenome_annotation.replace("_RC", "").replace(".gff", ".faa")
    logging.info(f"mitogenome_faa: {mitogenome_faa}") #debug
    frameshifts = find_frameshifts_mitos(mitogenome_faa)
    logging.info(f"frameshifts: {frameshifts}") #debug
    seq_len, num_genes = get_mitos_stats(mitogenome_annotation, mitogenome_fasta)
    logging.info(f"get_mitos_stats output: seq_len={seq_len} | num_genes={num_genes}") #debug
    contig_dir = os.path.join("potential_contigs", contig_id)
    logging.info(f"contig_dir: {contig_dir}") #debug
    #mitogenome_location = os.path.join(contig_dir, mitogenome_rotated_gb)
    mitogenome_annotation_location = os.path.join(contig_dir, mitogenome_rotated_annotation)
    logging.info(f"mitogenome_annotation_location: {mitogenome_annotation_location}") #debug
    is_circ = get_circularization_info(contig_id)
    logging.info(f"is_circ: {is_circ}") #debug
    if not frameshifts:
        all_frameshifts = "No frameshift found"
    elif len(frameshifts)==1:
        all_frameshifts = "".join(frameshifts)
    elif len(frameshifts)>1:
        all_frameshifts = ";".join(frameshifts)
    with open(f"{contig_id}.individual.stats", "w") as outfile:
        outfile.write("\t".join([contig_id, all_frameshifts, mitogenome_annotation_location, str(seq_len), str(num_genes), str(is_circ)+"\n"]))

if __name__ == "__main__":
    
    if sys.argv[1] == "-h":
        print("""Usage:
        arg1 = contigs
        arg2 = related_mito_fasta
        arg3 = related_mito_gbk
        arg4 = genetic code
        """)
        sys.exit()

    contigs = sys.argv[1].split()
    print(f"Contigs: {contigs}")


    rel_mito_len = getMitoLength.get_mito_length(sys.argv[2])
    print("Length of related mitogenome is: {} bp".format(rel_mito_len))

    # calculate maximum contig size accepted by mitofinder when annotating the contigs
    max_contig_size = 5*rel_mito_len

    related_gbk = sys.argv[3]
    gen_code = sys.argv[4]

    partial_annotate_mito = functools.partial(annotate_mito, max_contig_size, related_gbk, gen_code)
    with concurrent.futures.ProcessPoolExecutor() as executor:
        executor.map(partial_annotate_mito, contigs)
    print("Finished annotation!")
