import logging
import os
import pandas as pd
import sys
import subprocess
import argparse
import alignContigs
import plot_coverage
import plot_coverage_final_mito
import PIL
import numpy as np

def get_contigs_headers(in_fasta):
    """
    Takes a multifasta file and returns a dictionary with the
    contigs_ids as keys and the contigs headers as values
    """
    
    contigs_info = {}
    with open(in_fasta) as f:
        for line in f:
            if line.startswith(">"):
                contig_header = line.strip().split()[0].replace(">", "")
                contig_id = contig_header.replace("_rc", "").replace("_rotated", "")
                contigs_info[contig_id] = contig_header
    return contigs_info

def get_contigs_to_map():
    """
    Iterates over `contigs_stats.tsv` to list all potential contigs
    and returns a list with the filenames of their rotated FASTA files
    """

    try:
        with open("contigs_stats.tsv") as f:
            pass
    except FileNotFoundError:
        sys.exit("No contigs_stats.tsv. Cannot find list of potential contigs to map against")
    
    contigs_to_map = []
    with open("contigs_stats.tsv", "r") as f:
        next(f)
        next(f)
        for line in f:
            contig_id = line.split()[0]
            if contig_id == "final_mitogenome":
                contig_fasta = "final_mitogenome.fasta"
            else:
                contig_fasta = contig_id + ".mitogenome.rotated.fa"
            try:
                with open(contig_fasta) as f2:
                    pass
            except FileNotFoundError:
                    sys.exit(f"File {contig_fasta} not found")
            contigs_to_map.append(contig_fasta)
    
    return contigs_to_map


def map_potential_contigs(in_reads, contigs, threads=1, covMap=20): 
    """
    Args:
        in_reads (list): reads file to be mapped against contigs
        contigs (list): list of contigs FASTA files to concatenate and map the reads against
        threads (int): number of threads to be used for computation
        covMap (int): minimum mapping quality to filter reads when building final coverage plot
    Return:
        (str) Filename of the sorted mapping file in BAM format
    """

    concatenated_contigs = alignContigs.concatenate_contigs(contigs, out_file="all_potential_contigs.fa")

    try:
        with open("all_potential_contigs.fa") as f:
            pass
    except FileNotFoundError:
        sys.exit("""No all_potential_contigs.fa file.
        An error may have occurred when concatenating potential contigs into a single FASTA file.""")
    #logging.info(f"map_potential_contigs done. Concatenated_contigs: {concatenated_contigs}")
    #sys.exit(0)
     
   # with open('concated-mitos.fasta', 'w') as outfile:
   #     for f in lista:
   #         with open(f) as infile:
   #             outfile.write(infile.read())
    
    # map reads against concatenated contigs file
    minimap_cmd = ["minimap2", "-t", str(threads), "--secondary=no", "-ax", "map-pb", concatenated_contigs] + in_reads
    samtools_cmd = ["samtools", "view", "-@", str(threads), "-b", "-F4", "-F", "0x800", "-q", str(covMap), "-o", "HiFi-vs-potential_contigs.bam"]
    logging.info("Reads mapping:")
    logging.info(" ".join(minimap_cmd) + " | " + " ".join(samtools_cmd))
    minimap = subprocess.Popen(minimap_cmd, stdout=subprocess.PIPE, stderr=subprocess.DEVNULL)
    subprocess.run(samtools_cmd, stderr=subprocess.STDOUT, stdin=minimap.stdout)
    minimap.wait()
    minimap.stdout.close()

    # sorting and creating index for the mapping file
    try:
        with open("HiFi-vs-potential_contigs.bam") as f:
            pass
    except FileNotFoundError:
        sys.exit("""No HiFi-vs-potential_contigs.bam file.
        An error may have occurred when mapping reads to potential contigs""")
    
    sort_cmd = ["samtools", "sort", "-@", str(threads), "HiFi-vs-potential_contigs.bam", "-o", "HiFi-vs-potential_contigs.sorted.bam"]
    logging.info("Sorting mapping file:")
    logging.info(" ".join(sort_cmd))
    subprocess.run(sort_cmd, stderr=subprocess.STDOUT)
    try:
        with open("HiFi-vs-potential_contigs.sorted.bam") as f:
            pass
    except FileNotFoundError:
        sys.exit("""No HiFi-vs-potential_contigs.sorted.bam file.
        An error may have occurred when sorting the HiFi-vs-potential_contigs.bam file""")
    
    index_cmd = ["samtools", "index", "HiFi-vs-potential_contigs.sorted.bam"]
    logging.info("Indexing sorted mapping file:")
    logging.info(" ".join(index_cmd))
    subprocess.run(index_cmd, stderr=subprocess.STDOUT)
    try:
        with open("HiFi-vs-potential_contigs.sorted.bam.bai") as f:
            pass
    except FileNotFoundError:
        sys.exit("""No HiFi-vs-potential_contigs.sorted.bam.bai file.
        An error may have occurred when indexing the HiFi-vs-potential_contigs.sorted.bam file""")
    
    return "HiFi-vs-potential_contigs.sorted.bam"

def map_final_mito(in_reads, threads=1, covMap=20): 
    """
    Args:
        in_reads (list): reads file to be mapped against contigs
        threads (int): number of threads to be used for computation
        covMap (int): minimum mapping quality to filter reads when building final coverage plot
    Return:
        (str) Filename of the sorted mapping file in BAM format
    """

    try:
        with open("final_mitogenome.fasta") as f:
            pass
    except FileNotFoundError:
        sys.exit("""No final_mitogenome.fasta file.
        An error may have occurred when choosing the representative contig.""")
    
    # map reads 
    minimap_cmd = ["minimap2", "-t", str(threads), "--secondary=no", "-ax", "map-pb", "final_mitogenome.fasta"] + in_reads
    samtools_cmd = ["samtools", "view", "-@", str(threads), "-b", "-F4", "-F", "0x800", "-q", str(covMap), "-o", "HiFi-vs-final_mitogenome.bam"]
    logging.info("Reads mapping:")
    logging.info(" ".join(minimap_cmd) + " | " + " ".join(samtools_cmd))
    minimap = subprocess.Popen(minimap_cmd, stdout=subprocess.PIPE, stderr=subprocess.DEVNULL)
    subprocess.run(samtools_cmd, stderr=subprocess.STDOUT, stdin=minimap.stdout)
    minimap.wait()
    minimap.stdout.close()

    # sorting and creating index for the mapping file
    try:
        with open("HiFi-vs-final_mitogenome.bam") as f:
            pass
    except FileNotFoundError:
        sys.exit("""No HiFi-vs-final_mitogenome.bam file.
        An error may have occurred when mapping reads to final mitogenome""")
    
    sort_cmd = ["samtools", "sort", "-@", str(threads), "HiFi-vs-final_mitogenome.bam", "-o", "HiFi-vs-final_mitogenome.sorted.bam"]
    logging.info("Sorting mapping file:")
    logging.info(" ".join(sort_cmd))
    subprocess.run(sort_cmd, stderr=subprocess.STDOUT)
    try:
        with open("HiFi-vs-final_mitogenome.sorted.bam") as f:
            pass
    except FileNotFoundError:
        sys.exit("""No HiFi-vs-final_mitogenome.sorted.bam file.
        An error may have occurred when sorting the HiFi-vs-final_mitogenome.bam file""")
    
    index_cmd = ["samtools", "index", "HiFi-vs-final_mitogenome.sorted.bam"]
    logging.info("Indexing sorted mapping file:")
    logging.info(" ".join(index_cmd))
    subprocess.run(index_cmd, stderr=subprocess.STDOUT)
    try:
        with open("HiFi-vs-final_mitogenome.sorted.bam.bai") as f:
            pass
    except FileNotFoundError:
        sys.exit("""No HiFi-vs-final_mitogenome.sorted.bam.bai file.
        An error may have occurred when indexing the HiFi-vs-final_mitogenome.sorted.bam file""")
    
    return "HiFi-vs-final_mitogenome.sorted.bam"

def split_mapping_by_contig(all_contigs_mapping, contigs_headers, threads=1):
    """
    Takes a mapping file with reads mapped to a multifasta file
    and creates individual mapping files for each contig from 
    the `contigs_headers` dictionary. 
    `contigs_headers` contains contigs_ids as values and contigs 
    headers as values
    """
    
    mapped_contigs = [] 
    for contig_id in contigs_headers:
        contig_mapping_file = f"{contig_id}.bam"
        samtools_cmd = ["samtools", "view", "-@", str(threads), "-b", "-o", contig_mapping_file, all_contigs_mapping, contigs_headers[contig_id]]
        logging.info(f"Retrieve BAM for contig {contig_id}:")
        logging.info(" ".join(samtools_cmd))
        subprocess.run(samtools_cmd, stderr=subprocess.STDOUT)
         
        try:
            with open(contig_mapping_file) as f:
                pass
        except FileNotFoundError:
            sys.exit(f"""No {contig_mapping_file} file.
            An error may have occurred when subsetting {all_contigs_mapping} file for {contig_id}""")
        
        contig_mapping_sorted_file = f"{contig_id}.sorted.bam"
        sort_cmd = ["samtools", "sort", "-@", str(threads), contig_mapping_file, "-o", contig_mapping_sorted_file]
        logging.info("Sorting mapping file:")
        logging.info(" ".join(sort_cmd))
        subprocess.run(sort_cmd, stderr=subprocess.STDOUT)
        
        try:
            with open(contig_mapping_sorted_file) as f:
                pass
        except FileNotFoundError:
            sys.exit(f"""No {contig_mapping_sorted_file} file.
            An error may have occurred when sorting the {contig_mapping_file} file""")
        
        contig_mapping_index_file = contig_mapping_sorted_file + ".bai"
        index_cmd = ["samtools", "index", contig_mapping_sorted_file]
        logging.info("Indexing sorted mapping file:")
        logging.info(" ".join(index_cmd))
        subprocess.run(index_cmd, stderr=subprocess.STDOUT)
        try:
            with open(contig_mapping_index_file) as f:
                pass
        except FileNotFoundError:
            sys.exit(f"""No {contig_mapping_index_file} file.
            An error may have occurred when indexing the {contig_mapping_sorted_file} file""")

        mapped_contigs.append(contig_id)
    return mapped_contigs

def merge_images(img_list, out_file):

    imgs = [ PIL.Image.open(i) for i in img_list ]
    imgs_comb = np.vstack([np.asarray(i) for i in imgs])
    imgs_comb = PIL.Image.fromarray(imgs_comb)
    imgs_comb.save(out_file)

def create_coverage_plot(mapped_contigs, winSize, repr_contig, is_final_mito=False):

    print("create_coverage_plot function started:") #debug
    coverage_plots = []
    for contig_id in mapped_contigs:
        if is_final_mito:
            genome_filename = plot_coverage_final_mito.make_genome_file(contig_id)
        else:
            genome_filename = plot_coverage.make_genome_file(contig_id)
        print(f"genome_filename: {genome_filename}") #debug
        genome_windows_filename = plot_coverage.make_genome_windows(genome_filename, winSize)
        print(f"genome_windows_filename: {genome_windows_filename}") #debug
        if is_final_mito:
            windows_depth_filename = plot_coverage_final_mito.get_windows_depth(genome_windows_filename, "HiFi-vs-final_mitogenome.bam")
        else:
            windows_depth_filename = plot_coverage.get_windows_depth(genome_windows_filename, f"{contig_id}.bam")
        print(f"filename: {windows_depth_filename}") #debug
        if is_final_mito:
            coverage_plot_filename = plot_coverage_final_mito.plot_coverage("final_mitogenome", windows_depth_filename, winSize, isFinalMito=True)
            return "final_mitogenome.coverage.png" 
        else:
            if contig_id == repr_contig:
                coverage_plot_filename = plot_coverage.plot_coverage(contig_id, windows_depth_filename, winSize, isFinalMito=True)
            else:
                coverage_plot_filename = plot_coverage.plot_coverage(contig_id, windows_depth_filename, winSize)
        print(f"coverage_plot_filename: {coverage_plot_filename}") #debug
        coverage_plots.append(coverage_plot_filename)

        print(f"coverage_plots: {coverage_plots}")
        
        ### potential identation problem
        merge_images(coverage_plots, "coverage_plot.png")
        try:
            with open("coverage_plot.png", "r") as f:
                pass
        except FileNotFoundError:
           sys.exit("coverage_plot.png file not found.") 

    return "coverage_plot.png"

if __name__ == "__main__":
    in_reads = sys.argv[1]
    concatenate_map(in_reads)

