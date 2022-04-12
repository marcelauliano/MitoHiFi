import concurrent.futures
from concurrent.futures import wait
import logging
import shutil
import subprocess
import time
import warnings
import pandas as pd
import argparse
import parse_blast
from Bio import SeqIO
import sys
import os
import cleanUpCWD
import fetch
import filterfasta
import findFrameShifts
import fixContigHeaders
import functools
import rotation
import getMitoLength
import getReprContig
from parallel_annotation import process_contig, process_contig_02
import shlex
from circularizationCheck import circularizationCheck, get_circo_mito
import alignContigs


def main():
    
    __version__ = '2.14.2'
    start_time = time.time()

    parser = argparse.ArgumentParser(prog='MitoHiFi')
    parser._action_groups.pop()
    required = parser.add_argument_group('required arguments')
    optional = parser.add_argument_group('optional arguments')
    mutually_exclusive_group = required.add_mutually_exclusive_group(required=True)    
    mutually_exclusive_group.add_argument("-r", help= "-r: Pacbio Hifi Reads from your species", metavar='<reads>.fasta')
    mutually_exclusive_group.add_argument("-c", help= "-c: Assembled fasta contigs/scaffolds to be searched to find mitogenome", metavar='<contigs>.fasta')
    required.add_argument("-f", help= "-f: Close-related Mitogenome is fasta format", required = "True", metavar='<relatedMito>.fasta')
    required.add_argument("-g", help= "-k: Close-related species Mitogenome in genebank format", required = "True", metavar='<relatedMito>.gbk')
    required.add_argument("-t", help= "-t: Number of threads for (i) hifiasm and (ii) the blast search", required = "True", type=int, metavar='<THREADS>')    
    optional.add_argument("-d", help="-d: debug mode to output additional info on log", action="store_true")    
    optional.add_argument("-a", help="-a: Choose between animal (default) or plant", default="animal", choices=["animal", "plant"])
    optional.add_argument("-p", help="-p: Percentage of query in the blast match with close-related mito", type=int, default=50, metavar='<PERC>')
    optional.add_argument("-m", help="-m: Number of bits for HiFiasm bloom filter [it maps to -f in HiFiasm] (default = 0)", type=int, default=0, metavar='<BLOOM FILTER>')
    optional.add_argument('--circular-size', help='Size to consider when checking for circularization', type=int, default=220)
    optional.add_argument('--circular-offset', help='Offset from start and finish to consider when looking for circularization', type=int, default=40)
    optional.add_argument('-v', '--version', action='version', version='%(prog)s ' + __version__)
    optional.add_argument("-o", help="""-o: Organism genetic code following NCBI table (for mitogenome annotation):
    1. The Standard Code 2. The Vertebrate MitochondrialCode 3. The Yeast Mitochondrial Code 
    4. The Mold,Protozoan, and Coelenterate Mitochondrial Code and the Mycoplasma/Spiroplasma Code 5. The Invertebrate Mitochondrial Code 
    6. The Ciliate, Dasycladacean and Hexamita Nuclear Code 9. The Echinoderm and Flatworm Mitochondrial Code 10. The Euplotid Nuclear Code 
    11. The Bacterial, Archaeal and Plant Plastid Code 12. The Alternative Yeast Nuclear Code 13. The Ascidian Mitochondrial Code 
    14. The Alternative Flatworm Mitochondrial Code 16. Chlorophycean Mitochondrial Code 21. Trematode Mitochondrial Code 
    22. Scenedesmus obliquus Mitochondrial Code 23. Thraustochytrium Mitochondrial Code 24. Pterobranchia Mitochondrial Code 
    25. Candidate Division SR1 and Gracilibacteria Code 
        """, type=str, default='1', metavar='<GENETIC CODE>')
    args = parser.parse_args()
    
    # Set log message format
    FORMAT='%(asctime)s [%(levelname)s] %(message)s'

    if args.d: # If in debug mode
        logging.basicConfig(level=logging.DEBUG, stream=sys.stdout,
                            format=FORMAT, datefmt='%Y-%m-%d %H:%M:%S')
    else:
        logging.basicConfig(level=logging.INFO, stream=sys.stdout,
                            format=FORMAT, datefmt='%Y-%m-%d %H:%M:%S')
    
    # Welcome message
    logging.info("Welcome to MitoHifi v2. Starting pipeline...")
    logging.debug("Running MitoHiFi on debug mode.")
    
    # Measure the length of the related mitogenome 
    rel_mito_len = getMitoLength.get_mito_length(args.f)
    rel_mito_num_genes = getMitoLength.get_mito_genes(args.g)
    logging.info("Length of related mitogenome is: {} bp".format(rel_mito_len))
    logging.info("Number of genes on related mitogenome: {}".format(rel_mito_num_genes))
    
    # If input are reads, map them to the related mitogenome and assemble the mapped ones
    if args.r:
        logging.info("Running MitoHifi pipeline in reads mode...")
        logging.info("1. First we map your Pacbio HiFi reads to the close-related mitogenome")
        minimap_cmd = ["minimap2", "-t", str(args.t), "--secondary=no", "-ax", "map-pb", args.f] + shlex.split(args.r) 
        samtools_cmd = ["samtools", "view", "-@", str(args.t), "-S", "-b", "-F4", "-F", "0x800"] 
        logging.info(" ".join(minimap_cmd) + " | " + " ".join(samtools_cmd) + " > reads.HiFiMapped.bam")        
        minimap = subprocess.Popen(minimap_cmd, stdout=subprocess.PIPE, stderr=subprocess.DEVNULL)
        mapped_reads_f = open("reads.HiFiMapped.bam", "w")
        subprocess.run(samtools_cmd, stderr=subprocess.STDOUT, stdin=minimap.stdout, stdout=mapped_reads_f)
        minimap.wait()
        minimap.stdout.close()
        
        try:
            f = open("reads.HiFiMapped.bam")
        except FileNotFoundError:
            sys.exit("""No reads.HiFiMapped.bam file.
            An error may have occurred when mapping reads to the close-related mitogenome""")
        finally:
            f.close()

        logging.info("2. Now we filter out any mapped reads that are larger than the reference mitogenome to avoid NUMTS")
        bam2fasta_cmd = ["samtools", "fasta", "reads.HiFiMapped.bam"]
        logging.info("2.1 First we convert the mapped reads from BAM to FASTA format:")
        logging.info(" ".join(bam2fasta_cmd) + " > gbk.HiFiMapped.bam.fasta")
        mapped_fasta_f = open("gbk.HiFiMapped.bam.fasta", "w")
        subprocess.run(bam2fasta_cmd, stdout=mapped_fasta_f, stderr=subprocess.DEVNULL)
        before_filter = fetch.get_num_seqs("gbk.HiFiMapped.bam.fasta")
        logging.info(f"Total number of mapped reads: {before_filter}")

        logging.info(f"2.2 Then we filter reads that are larger than {rel_mito_len} bp")
        filterfasta.filterFasta(minLength=rel_mito_len, neg=True, inStream="gbk.HiFiMapped.bam.fasta",
                                outPath="gbk.HiFiMapped.bam.filtered.fasta", log=False)
        
        try:
            f = open("gbk.HiFiMapped.bam.filtered.fasta")
        except FileNotFoundError:
            sys.exit("""No gbk.HiFiMapped.bam.filtered.fasta file.
            An error may have occurred when filtering reads larger than the reference mitogenome""")
        finally:
            f.close()

        after_filter = fetch.get_num_seqs("gbk.HiFiMapped.bam.filtered.fasta")
        logging.info(f"Number of filtered reads: {after_filter}")

        logging.info("3. Now let's run hifiasm to assemble the mapped and filtered reads!")
        
        hifiasm_cmd = ["hifiasm", "--primary", "-t", str(args.t), "-f", str(args.m), 
                    "-o", "gbk.HiFiMapped.bam.filtered.assembled",
                    "gbk.HiFiMapped.bam.filtered.fasta"]

        logging.info(" ".join(hifiasm_cmd))
        with open("hifiasm.log", "w") as hifiasm_log_f:
            subprocess.run(hifiasm_cmd, stderr=subprocess.STDOUT, stdout=hifiasm_log_f)       
        
        try:
            f1 = open("gbk.HiFiMapped.bam.filtered.assembled.p_ctg.gfa")
            f2 = open("gbk.HiFiMapped.bam.filtered.assembled.a_ctg.gfa")
        except FileNotFoundError:
            sys.exit("""No gbk.HiFiMapped.bam.filtered.assembled.[a/p]_ctg.gfa file(s).
            An error may have occurred when assembling reads with HiFiasm.""")
        finally:
            f1.close()
            f2.close()

        gfa2fa_script = os.path.join(os.path.dirname(os.path.realpath(__file__)),"gfa2fa") # gets path to gfa2fa script
        
        with open("gbk.HiFiMapped.bam.filtered.assembled.p_ctg.fa", "w") as p_ctg_f:
            subprocess.run([gfa2fa_script, "gbk.HiFiMapped.bam.filtered.assembled.p_ctg.gfa"], stdout=p_ctg_f)
        with open("gbk.HiFiMapped.bam.filtered.assembled.a_ctg.fa", "w") as a_ctg_f:
            subprocess.run([gfa2fa_script, "gbk.HiFiMapped.bam.filtered.assembled.a_ctg.gfa"], stdout=a_ctg_f)
        
        with open("hifiasm.contigs.fasta", "w") as hifiasm_f:
            subprocess.run(["cat", "gbk.HiFiMapped.bam.filtered.assembled.p_ctg.fa", "gbk.HiFiMapped.bam.filtered.assembled.a_ctg.fa"], stdout=hifiasm_f)
        
        contigs = "hifiasm.contigs.fasta"
    
    else:
        logging.info("Running MitoHifi pipeline in contigs mode...")
       
        logging.info("1. Fixing potentially conflicting FASTA headers")
        original_contigs = args.c
        fixContigHeaders.fix_headers(original_contigs, "fixed_header_contigs.fasta")
        
        os.remove(original_contigs) # remove original contig file  
        shutil.move("fixed_header_contigs.fasta", original_contigs) # replace original contigs file by the version that has the headers fixed
        
        contigs = original_contigs

    # Set number for the current step (for improving understanding of the log)
    # On reads mode, it should be 4; on contigs mode, 2    
    if args.r:
        step = 4
    else:
        step = 2 

    logging.info(f"{step}. Let's run the blast of the contigs versus the close-related mitogenome")

    makeblastdb_cmd = ["makeblastdb", "-in", args.f, "-dbtype", "nucl"]
    logging.info(f"{step}.1. Creating BLAST database:")
    logging.info(" ".join(makeblastdb_cmd))
    subprocess.run(makeblastdb_cmd, stderr=subprocess.DEVNULL, stdout=subprocess.DEVNULL)
    logging.info("Makeblastdb done.")
    
    blast_cmd = ["blastn", "-query", contigs, "-db", args.f, "-num_threads", str(args.t),
                "-out", "contigs.blastn", "-outfmt", "6 std qlen slen"]
    logging.info(f"{step}.2. Running blast of contigs against close-related mitogenome:")
    logging.info(" ".join(blast_cmd))
    subprocess.run(blast_cmd, stderr=subprocess.STDOUT)
    logging.info("Blast done.")

    try:
        f = open("contigs.blastn")
    except FileNotFoundError:
        sys.exit("""No contigs.blastn file.
        An error may have occurred when BLASTing contigs against close-related mitogenome""")
    finally:
        f.close()    

    step += 1
    logging.info(f"{step}. Filtering BLAST output to select target sequences")

    #the next script parses a series of conditions to exclude blast with NUMTs. 
    if args.a == "plant":
        ## if species is a plant, set minimum query percentage equal to 0% of related mito 
        ## and maximum query lenght 10 times the lenght of the related
        parse_blast.parse_blast(query_perc=args.p, min_query_perc=0, max_query_len=10)
    else:
        ## if species is an animal, set minimum query percentage equal to 80% of related mito
        ## and maximum query length 5 times the length of the related (default values from 
        ## parse_blast function
        parse_blast.parse_blast(query_perc=args.p, max_query_len=5)

    # select contigs to be circularized
    # first look for contigs in parsed_blast.txt
    contigs_ids = parse_blast.get_contigs_ids("parsed_blast.txt")

    # if we don't find contigs in parse_blast.txt 
    # look for contigs in parsed_blast_all.txt
    if len(contigs_ids) == 0:
        contigs_ids = parse_blast.get_contigs_ids("parsed_blast_all.txt")

    # if we can't find any contigs even in parsed_blast_all.txt, then we exit the pipeline
    if len(contigs_ids) == 0:
        sys.exit("""Attention!
'parsed_blast.txt' and 'parsed_blast_all.txt' files are empty.
The pipeline has stopped !! You need to run further scripts to check if you have mito reads pulled to a large NUMT!""")

    logging.info("Filtering BLAST finished. A list of the filtered contigs was saved on ./contigs_filtering/contigs_ids.txt file")

    # records all contigs kept for the downstream steps in a file called 'contigs_ids.txt'
    with open("contigs_ids.txt", "w") as f:
        for contig_id in contigs_ids:
            f.write(contig_id + "\n")

    # removes file that contains history of circularization of it already exists
    try:
        os.remove('all_contigs.circularisationCheck.txt')
    except OSError:
        pass

    step += 1
    logging.info(f"{step}. Now we are going to circularize, annotate and rotate each filtered contig. Those are potential mitogenome(s).")
    
    # Set maximum contig size accepted by mitofinder when annotating the contigs
    max_contig_size = 5*rel_mito_len

    threads_per_contig = 1
    if args.t // len(contigs_ids) > 1:
        threads_per_contig = args.t // len(contigs_ids)

    logging.debug(f"Threads per contig={threads_per_contig}")
    logging.debug(f"Thresholds for circularization: circular size={args.circular_size} | circular offset={args.circular_offset}")
    logging.debug(f"Thresholds for annotation (MitoFinder): maximum contig size={max_contig_size}")
    partial_process_contig = functools.partial(process_contig, threads_per_contig,
                                               args.circular_size, args.circular_offset,
                                               contigs, max_contig_size, args.g, args.o)
    with concurrent.futures.ProcessPoolExecutor() as executor:
        executor.map(partial_process_contig, contigs_ids)
    
    tRNA_ref = fetch.get_ref_tRNA() 
    logging.debug(f"tRNA to be used as reference for rotation: {tRNA_ref}") 
    
    partial_process_contig_02 = functools.partial(process_contig_02, tRNA_ref,
                                                threads_per_contig, args.circular_size,
                                                args.circular_offset, contigs, max_contig_size,
                                                args.g, args.o)
    with concurrent.futures.ProcessPoolExecutor() as executor:
        executor.map(partial_process_contig_02, contigs_ids)

    # concatenates circularization histories from all contigs into 
    # a single concatenate history file 
    circularization_history_files = []
    for curr_file in os.listdir('.'):
        if curr_file.endswith('.circularisationCheck.txt'):
            circularization_history_files.append(curr_file)
    logging.debug(f"circularization_history_files: {circularization_history_files}") # for debugging
    with open("all_contigs.circularisationCheck.txt", "a") as outfile:
        for circ_hist_file in circularization_history_files:
            contig_id = circ_hist_file.split('.')[0]
            with open(circ_hist_file, "r") as infile:
                for line in infile:
                    outfile.write("\t".join([contig_id, line.strip()+"\n"]))
    
    step += 1
    #align final mitogenome rotated contigs
    logging.info(f"{step}. Now the rotated contigs will be aligned")
    # list all final rotated contigs 
    contigs_files = []
    for curr_file in os.listdir('.'):
        if curr_file.endswith('mitogenome.rotated.fa'):
            contigs_files.append(curr_file)
    if not contigs_files:
        sys.exit("""No rotated contigs found. 
        An error has possibly occurred during annotation and/or rotation of your contigs""")
    # first concatenate all rotated contigs into a single multifasta file
    logging.info("List of contigs that will be aligned: " + str(contigs_files) + "\n")
    concat_fasta = alignContigs.concatenate_contigs(contigs_files)
    # then run MAFFT alignment between the rotated contigs using the multifasta as input and clustal as output format
    alignContigs.mafft_align(multifasta_file=concat_fasta, threads=args.t, clustal_format=True)
    logging.info("Alignment done and saved at ./final_mitogenome_choice/all_mitogenomes.rotated.aligned.fa\n")

    try:
        f = open("all_mitogenomes.rotated.aligned.fa")
    except FileNotFoundError:
        sys.exit("""No all_mitogenomes.rotated.aligned.fa file.
        An error may have occurred when aligning the potential contigs""")
    finally:
        f.close()

    step += 1
    logging.info(f"{step}. Now we will choose the most representative contig" + "\n")
    repr_contig_id, repr_contig_cluster = getReprContig.get_repr_contig("all_mitogenomes.rotated.fa", rel_mito_len, args.t, args.d)
    logging.info(f"Representative contig is {repr_contig_id} that belongs to {repr_contig_cluster}. This contig will be our final mitogenome. See all contigs and clusters in cdhit.out.clstr")
    
    repr_contig_fasta = repr_contig_id + ".mitogenome.rotated.fa"
    repr_contig_get_gb = ["mitofinder", "--new-genes", "--max-contig-size",
                        str(max_contig_size), "-j", "final_mitogenome.annotation",
                        "-a", repr_contig_fasta, "-r", args.g, "-o", args.o, "-p", str(args.p),
                        "--circular-size", "8000"]
    subprocess.run(repr_contig_get_gb, stderr=subprocess.DEVNULL, stdout=subprocess.DEVNULL)

    final_fasta = os.path.join("final_mitogenome.annotation", "final_mitogenome.annotation_MitoFinder_mitfi_Final_Results", "final_mitogenome.annotation_mtDNA_contig.fasta")
    final_gbk = os.path.join("final_mitogenome.annotation", "final_mitogenome.annotation_MitoFinder_mitfi_Final_Results", "final_mitogenome.annotation_mtDNA_contig.gb")
    
    # Generate contigs stats
    step += 1
    logging.info(f"""{step}. Calculating final stats for final mitogenome and other potential contigs.
    Stats will be saved on contigs_stats.tsv file.""")
    ## Print first three lines (comment and header)
    frameshifts = findFrameShifts.find_frameshifts(final_gbk)  
    contig_len, num_genes = findFrameShifts.get_gb_stats(final_gbk)
    is_circ = getReprContig.get_circularization_info(repr_contig_id)
    if not frameshifts:
        all_frameshifts = "No frameshift found"
    elif len(frameshifts)==1:
        all_frameshifts = "".join(frameshifts)
    elif len(frameshifts)>1:
        all_frameshifts = ";".join(frameshifts)
    with open("contigs_stats.tsv", "w") as f: 
        f.write(f"# Related mitogenome is {rel_mito_len} bp long and has {rel_mito_num_genes} genes\n")
        f.write("\t".join(["contig_id", "frameshifts_found", "genbank_file", 
        "length(bp)", "number_of_genes", "was_circular\n"]))
        f.write("\t".join(["final_mitogenome", all_frameshifts, "final_mitogenome.gb",
        contig_len, num_genes, str(is_circ)+"\n"]))
    ## Iterate over each contig and print its info (ID, framshifts and genbank file used 
    ## to search for the frameshifts)
    
    # get list with all individual contig stats files 
    contigs_stats_files = []
    for curr_file in os.listdir('.'):
        if curr_file.endswith('.individual.stats'):
            # skips addition of representative contig, which is the 
            # same as the final_mitogenome
            if curr_file.split('.')[0] != repr_contig_id: 
                contigs_stats_files.append(curr_file)
    
    with open("contigs_stats.tsv", "a") as outfile:
        for contig_stats in contigs_stats_files:
            with open(contig_stats, "r") as infile:
                shutil.copyfileobj(infile, outfile)

    # copying final FASTA and GBK to working directory
    shutil.copy(final_fasta, "final_mitogenome.fasta")
    shutil.copy(final_gbk, "final_mitogenome.gb")

    # cleaning up working directory 
    cleanUpCWD.clean_up_work_dir(contigs_ids)

    logging.info("Pipeline finished!" )
    runtime= time.time() - start_time
    logging.info(f"Run time: {runtime:.2f} seconds")

if __name__ == '__main__':
    main()
