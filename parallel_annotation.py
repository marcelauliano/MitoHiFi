import os
import concurrent.futures
import sys
import functools
import subprocess
import getMitoLength
import logging
import warnings 
import rotation
from circularizationCheck import get_circo_mito
from getReprContig import get_circularization_info
import findFrameShifts
import filterfasta

def process_contig(threads_per_contig, circular_size, circular_offset, contigs, max_contig_size, rel_gbk, gen_code, contig_id): 
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
    logging.info(f"Started {contig_id} (MitoFinder) annotation")
    mitofinder_cmd = ["mitofinder", "--new-genes", "--max-contig-size", str(max_contig_size),
                    "-j", contig_id+".annotation", "-a", contig_id+".mitogenome.fa",
                    "-r", rel_gbk, "-o", gen_code, "-p", str(threads_per_contig),
                    "--circular-size", "8000"] 
    subprocess.run(mitofinder_cmd, stderr=subprocess.DEVNULL, stdout=subprocess.DEVNULL)
    logging.info(f"{contig_id} annotation done. Annotation log saved on ./potential_contigs/{contig_id}/{contig_id}.annotation_MitoFinder.log")
    # rotates the mitogenome
    mitogenome_gb = os.path.join(contig_id + ".annotation", contig_id + ".annotation_MitoFinder_mitfi_Final_Results", contig_id + ".annotation_mtDNA_contig.gb") 
    if not os.path.isfile(mitogenome_gb):
        warnings.warn("Contig "+ contig_id + " does not have an annotation file, check MitoFinder's log")
        return
    
    trnas = rotation.get_trna_pos(mitogenome_gb)
    if trnas:
        with open(f"{contig_id}.trnas", "w") as outfile:
            for trna in trnas:
                outfile.write("\t".join([trna, trnas[trna][0], trnas[trna][1], "\n"]))
        return
    else:
        warnings.warn(f"No tRNA gene found in {mitogenome_gb}... Skipping contig {contig_id}")
        return

def process_contig_02(ref_tRNA, threads_per_contig, circular_size, circular_offset, contigs, max_contig_size, rel_gbk, gen_code, contig_id): 
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
    logging.info(f"Started {contig_id} rotation.")
    if not os.path.isfile(f"{contig_id}.trnas"):
        warnings.warn(f"Contig {contig_id} does not have annotated tRNAs, skipping it...")
        return
    
    with open(f"{contig_id}.trnas", "r") as infile:
        for line in infile:
            if line.strip().split("\t")[0] == ref_tRNA:
                start = int(line.strip().split("\t")[1])
                strand = int(line.strip().split("\t")[2])
    
    if start == None:
        warnings.warn(f"Reference gene {ref_tRNA} is not present in contig {contig_id}. Skipping contig...")
        return

    mitogenome_gb = os.path.join(contig_id + ".annotation", contig_id + ".annotation_MitoFinder_mitfi_Final_Results", contig_id + ".annotation_mtDNA_contig.gb")

    genome = contig_id + ".mitogenome.fa"
    if strand == -1:
        logging.info(f"{ref_tRNA} is at reverse complement of {contig_id}.mitogenome.fa")
        logging.info(f"For that reason we'll reverse complement {contig_id}.mitogenome.fa before the rotation")
        genome_rc = contig_id + "_RC.mitogenome.fa"
        rc = os.path.join(os.path.dirname(genome), genome_rc)
        rotation.make_rc(genome, rc)
        logging.info(f"Reverse complement generated: {contig_id}_RC.mitogenome.fa. Starting reverse complement annotation...")
        mitogenome_gb = rotation.annotate(os.path.dirname(genome), os.path.abspath(genome_rc), os.path.abspath(rel_gbk), contig_id, gen_code, max_contig_size, str(threads_per_contig))
        logging.info(f"Annotation of reverse complement for contig {contig_id} done")
        genome = rc
        if not os.path.isfile(mitogenome_gb):
            warnings.warn("Contig "+ contig_id + " does not have a reverse complemented annotation file, check MitoFinder's log")
            return

    rotation.rotate(genome, start, contig_id)
    
    try:
        f = open(mitogenome_gb)
    except FileNotFoundError:
        sys.exit(f"""Annotation file {mitogenome_gb} not found.
        An error may have occurred when annotating contig {contig_id}. Check MitoFinder""")
    finally:
        f.close()

    rotated_file = os.path.join(os.path.dirname(genome), contig_id + '.mitogenome.rotated.fa')
    logging.info(f"Rotation of {contig_id} done. Rotated is at {rotated_file}") 
    # check frameshifts in genes from contig and save findings to 
    # `{contig_id}.individual.stats` intermediate file
    frameshifts = findFrameShifts.find_frameshifts(mitogenome_gb)
    gb_len, num_genes = findFrameShifts.get_gb_stats(mitogenome_gb)
    contig_dir = os.path.join("potential_contigs", contig_id)
    mitogenome_location = os.path.join(contig_dir, mitogenome_gb)
    is_circ = get_circularization_info(contig_id)
    if not frameshifts:
        all_frameshifts = "No frameshift found"
    elif len(frameshifts)==1:
        all_frameshifts = "".join(frameshifts)
    elif len(frameshifts)>1:
        all_frameshifts = ";".join(frameshifts)
    with open(f"{contig_id}.individual.stats", "w") as outfile:
        outfile.write("\t".join([contig_id, all_frameshifts, mitogenome_location, gb_len, num_genes, str(is_circ)+"\n"]))

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
