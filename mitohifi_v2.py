import shutil
import subprocess
import warnings
import pandas as pd
import argparse
from datetime import date
import parse_blast
from Bio import SeqIO
import sys
import os
import cleanUpCWD
import filterfasta
import rotation
import getMitoLength
import getReprContig
from circularizationCheck import circularizationCheck
import alignContigs

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

def get_circo_mito(contig_id, circular_size, circular_offset):
    """
    It gets a contig ID and circularizes it, registering the circularization points
    args:
    the ID from the contig that we want to circularize
    """

    def cut_coords(contig_fasta, circularization_position, fasta_out):

        record=SeqIO.read(contig_fasta, "fasta")
        id = record.id
        get= record[circularization_position:]
        with open(mitogenome_fasta_out, "w") as f:
            f.write(get.format('fasta'))

    # find circularization point
    circularization_history = []
    contig_fasta = "".join([contig_id, ".mito.fa"])
    mitogenome_fasta_out = "".join([contig_id, ".mitogenome.fa"])
    circularization_checks = "".join([contig_id, ".circularisationCheck.txt"])

    circularization_info = circularizationCheck(resultFile=contig_fasta, circularSize=circular_size, circularOffSet=circular_offset)

    # writes circularization information to '[contig_id].circularisationCheck.txt' file
    with open(circularization_checks, "w") as f:
        f.write(str(circularization_info))
        circularization_history.append(str(circularization_info))
    
    is_circularizable = circularization_info[0]
    circularization_position = int(circularization_info[2])

    # if the contig is not circularizable, then create the "mitogenome.fasta" directly from it
    if not is_circularizable:
        record=SeqIO.read(contig_fasta, "fasta")
        id = record.id
        with open(mitogenome_fasta_out, "w") as f:  
            f.write(record.format('fasta'))        
    else:
        # if the contig is circularizable, then run circularization iteratively until the "mitogenome.fasta"
        # is no longer circularizable    
        while is_circularizable:
            cut_coords(contig_fasta, circularization_position, mitogenome_fasta_out)
            contig_fasta = mitogenome_fasta_out
            circularization_info = circularizationCheck(resultFile=contig_fasta, circularSize=circular_size, circularOffSet=circular_offset)
            is_circularizable = circularization_info[0]
            circularization_position = int(circularization_info[2])
            with open(circularization_checks, "a") as f:
                f.write("\n" + str(circularization_info))
                circularization_history.append(str(circularization_info))

    return circularization_history

def main():

    today = date.today()

    parser= argparse.ArgumentParser(add_help=False)
    mutually_exclusive_group = parser.add_mutually_exclusive_group(required=True)
    mutually_exclusive_group.add_argument("-r", help= "-r: Pacbio Hifi Reads from your species")
    mutually_exclusive_group.add_argument("-c", help= "-c: Assembled fasta contigs/scaffolds to be searched to find mitogenome")
    parser.add_argument("-h", "--help", action="help", default=argparse.SUPPRESS, help= "Print this help message") 
    parser.add_argument("-f", help= "-f: Close-related Mitogenome is fasta format", required = "True")
    parser.add_argument("-g", help= "-k: Close-related species Mitogenome in genebank format", required = "True")
    parser.add_argument("-t", help= "-t: Number of threads for different steps, such as hifiasm, blast, minimap2...", required = "True", type=int)
    parser.add_argument("-p", help="-p: Percentage of query (contig) in the blast match with close-related mito", type=int, default=50)
    parser.add_argument('--circular-size', help='Size to consider when checking for circularization', type=int, default=220)
    parser.add_argument('--circular-offset', help='Offset from start and finish to consider when looking for circularization', type=int, default=40)
    parser.add_argument("-o", help="""-o: Organism genetic code following NCBI table (for mitogenome annotation):
    1. The Standard Code 2. The Vertebrate MitochondrialCode 3. The Yeast Mitochondrial Code 
    4. The Mold,Protozoan, and Coelenterate Mitochondrial Code and the Mycoplasma/Spiroplasma Code 5. The Invertebrate Mitochondrial Code 
    6. The Ciliate, Dasycladacean and Hexamita Nuclear Code 9. The Echinoderm and Flatworm Mitochondrial Code 10. The Euplotid Nuclear Code 
    11. The Bacterial, Archaeal and Plant Plastid Code 12. The Alternative Yeast Nuclear Code 13. The Ascidian Mitochondrial Code 
    14. The Alternative Flatworm Mitochondrial Code 16. Chlorophycean Mitochondrial Code 21. Trematode Mitochondrial Code 
    22. Scenedesmus obliquus Mitochondrial Code 23. Thraustochytrium Mitochondrial Code 24. Pterobranchia Mitochondrial Code 
    25. Candidate Division SR1 and Gracilibacteria Code 
        """, type=str, default='1')
    args = parser.parse_args()

    print("MitoHifi v2.2" + "/n")
    print("Started at:" , today)
    
    # measure the length of the related mitogenome 
    rel_mito_len = getMitoLength.get_mito_length(args.f)
    print("Length of related mitogenome is: {} bp".format(rel_mito_len))
    # calculate maximum contig size accepted by mitofinder when annotating the contigs
    max_contig_size = 5*rel_mito_len

    # if input are reads, map them to the related mitogenome and assemble the mapped ones
    if args.r:
        print("\nRunning MitoHifi pipeline in reads mode\n")
       
        print("\nFirst we map your PacbioHiFi reads to the close-related mitogenome\n")

        minimap = "minimap2 "+ "-t " + str(args.t) + " " + "--secondary=no -ax map-pb " + args.f + " " + args.r + " | " + "samtools view -@ " + str(args.t) + "-S -b -F4 -F 0x800 > reads.HiFiMapped.bam"
        print(minimap)
        subprocess.run(minimap, shell=True)

        print("\nNow we filter out any mapped reads that are larger than the reference mitogenome to avoid NUMTS\n")
        samtools = "samtools fasta " + "reads.HiFiMapped.bam >  gbk.HiFiMapped.bam.fasta"
        print(samtools)
        subprocess.run(samtools, shell=True)

        filterfasta.filterFasta(minLength=rel_mito_len, neg=True, inStream="gbk.HiFiMapped.bam.fasta", outPath="gbk.HiFiMapped.bam.filtered.fasta")

        print("\nNow let's run hifiasm to assemble the mapped and filtered reads!\n")

        hifiasm = "hifiasm -t " + str(args.t) + " -o gbk.HiFiMapped.bam.filtered.assembled gbk.HiFiMapped.bam.filtered.fasta 2>hifiasm.log"
        print(hifiasm)
        subprocess.run(hifiasm, shell=True)

        get_p = "scripts/gfa2fa gbk.HiFiMapped.bam.filtered.assembled.p_ctg.gfa >  gbk.HiFiMapped.bam.filtered.assembled.p_ctg.fa"
        get_a = "scripts/gfa2fa gbk.HiFiMapped.bam.filtered.assembled.a_ctg.gfa >  gbk.HiFiMapped.bam.filtered.assembled.a_ctg.fa"
        cat = "cat gbk.HiFiMapped.bam.filtered.assembled.p_ctg.fa gbk.HiFiMapped.bam.filtered.assembled.a_ctg.fa > hifiasm.contigs.fasta"

        subprocess.run(get_p, shell=True)
        subprocess.run(get_a, shell=True)
        subprocess.run(cat, shell=True)
        
        contigs = "hifiasm.contigs.fasta"
    
    else:
        print("\nRunning MitoHifi pipeline in contigs mode\n")
        contigs = args.c
    
    print("\nLet's run the blast of the contigs versus the close-related mitogenome\n")

    makeblastdb = "makeblastdb -in " + args.f + " -dbtype nucl"
    print(makeblastdb)
    subprocess.run(makeblastdb, shell=True)
    print("\nmakeblastdb done. Running blast with the contigs\n")
    blast = "blastn -query " + contigs + " -db "+ args.f + " -num_threads " + str(args.t) + " -out contigs.blastn -outfmt '6 std qlen slen'"
    print(blast)
    subprocess.run(blast, shell=True)
    print("Blast done!" + "\n")

    #the next script parses a series of conditions to exclude blast with NUMTs. 
    parse_blast.parse_blast(query_perc=args.p)

    #We check for circularisation

    # select contigs to be circularized
    # first look for contigs in parsed_blast.txt
    contigs_ids = get_contigs_ids("parsed_blast.txt")

    # if we don't find contigs in parse_blast.txt 
    # look for contigs in parsed_blast_all.txt
    if len(contigs_ids) == 0:
        contigs_ids = get_contigs_ids("parsed_blast_all.txt")

    # if we can't find any contigs even in parsed_blast_all.txt, then we exit the pipeline
    if len(contigs_ids) == 0:
        sys.exit("""\n Attention! \n The 'parsed_blast.txt' and 'parsed_blast_all.txt' files are empty. The pipeline has stopped !! \n You need to run further scripts to check if you have mito reads pulled to a large NUMT!!""")

    # records all contigs kept for the downstream steps in a file called 'contigs_ids.txt'
    with open("contigs_ids.txt", "w") as f:
        for contig_id in contigs_ids:
            f.write(contig_id + "\n")

    # removes file that contains history of circularization of it already exists
    try:
        os.remove('all_contigs.circularisationCheck.txt')
    except OSError:
        pass
      
    print("\n" + "6-) Now we are going to circularize, annotate and rotate each contig which is a potential mitogenome" + "\n")
    #retrieves the FASTA files for each contig
    for contig_id in contigs_ids:
        print("\n" + "Working with contig: " + contig_id + "\n")
        # retrieves the FASTA files for each contig
        filterfasta.filterFasta(idList=[contig_id], inStream=contigs, outPath="".join([contig_id, ".mito.fa"]))
        # circularizes each contig and saves circularization history to a file
        circularization_history = get_circo_mito(contig_id, args.circular_size, args.circular_offset)
        for circularization_event in circularization_history: 
            with open('all_contigs.circularisationCheck.txt', 'a') as f:
                f.write("\t".join([contig_id, circularization_event, "\n"]))    
        # annotates mitogenome(s) using mitofinder
        print("Running mitofinder with maximum contig size of {} bp".format(max_contig_size))
        annotate_mito = "mitofinder --max-contig-size {} -j {}.annotation -a {}.mitogenome.fa -r {} -o {} -p {}".format(max_contig_size, contig_id, contig_id, args.g, args.o, args.t)
        print(annotate_mito)
        subprocess.run(annotate_mito, shell=True)
        # rotates the mitogenome
        mitogenome_gb = os.path.join(contig_id + ".annotation", contig_id + ".annotation_MitoFinder_mitfi_Final_Results", contig_id + ".annotation_mtDNA_contig.gb") 
        if not os.path.isfile(mitogenome_gb):
            warnings.warn("Contig "+ contig_id + " does not have an annotation file, check MitoFinder's output")
            continue
        start, strand = rotation.get_phe_pos(mitogenome_gb)
        genome = contig_id + ".mitogenome.fa"
        if start == None:
            warnings.warn('tRNA-Phe is not present in file ' + mitogenome_gb + '... Skipping contig ' + contig_id + '\n')
            continue
        new_gb = None
        if strand == -1:
            genome_rc = contig_id + "_RC.mitogenome.fa"
            rc = os.path.join(os.path.dirname(genome), genome_rc)
            rotation.make_rc(genome, rc)
            new_gb = rotation.annotate(os.path.dirname(genome), os.path.abspath(genome_rc), os.path.abspath(args.g), contig_id, args.o, max_contig_size, args.t)
            start, strand = rotation.get_phe_pos(new_gb)
            genome = rc
        rotation.rotate(genome, start, contig_id)
        print(' '.join(['Rotated to tRNA-Phe genome is at', os.path.join(os.path.dirname(genome), contig_id + '.mitogenome.rotated.fa')]))
        if new_gb:
            print("new_gb: ")
            print(new_gb)
            print('Mitogenome annotation is at ', new_gb)

    #align final mitogenome rotated contigs
    print("\n" + "7-) Now the final rotated contigs will be aligned" + "\n")
    # list all final rotated contigs 
    contigs_files = []
    for curr_file in os.listdir('.'):
        if curr_file.endswith('mitogenome.rotated.fa'):
            contigs_files.append(curr_file)
    # first concatenate all rotated contigs into a single multifasta file
    print("List of contigs that will be aligned: " + str(contigs_files) + "\n")
    concat_fasta = alignContigs.concatenate_contigs(contigs_files)
    # then run MAFFT alignment between the rotated contigs using the multifasta as input and clustal as output format
    alignContigs.mafft_align(multifasta_file=concat_fasta, threads=args.t, clustal_format=True)
  
    print("\n" + "8-) Now we will choose the most representative contig" + "\n")
    repr_contig_id, repr_contig_cluster = getReprContig.get_repr_contig("all_mitogenomes.rotated.fa", args.t)
    print("Representative contig is {} that belongs to {}. See all contigs and clusters in cdhit.out.clstr".format(repr_contig_id, repr_contig_cluster))
    
    repr_contig_fasta = repr_contig_id + ".mitogenome.rotated.fa"
    repr_contig_get_gb = ["mitofinder", "--max-contig-size", str(max_contig_size), "-j", "final_mitogenome.annotation", "-a", repr_contig_fasta, "-r", args.g, "-o", args.o, "-p", str(args.p)]
    subprocess.run(repr_contig_get_gb, shell=False)

    print("Pipeline finished!" )
    final_fasta = os.path.join("final_mitogenome.annotation", "final_mitogenome.annotation_MitoFinder_mitfi_Final_Results", "final_mitogenome.annotation_mtDNA_contig.fasta")
    final_gbk = os.path.join("final_mitogenome.annotation", "final_mitogenome.annotation_MitoFinder_mitfi_Final_Results", "final_mitogenome.annotation_mtDNA_contig.gb")
    
    # copying final FASTA and GBK to working directory
    shutil.copy(final_fasta, "final_mitogenome.fasta")
    shutil.copy(final_gbk, "final_mitogenome.gb")

    # cleaning up working directory 
    cleanUpCWD.clean_up_work_dir(contigs_ids)

if __name__ == '__main__':
    main()
