from Bio import SeqIO
import time
import argparse
import pandas as pd
import matplotlib.pyplot as plt
import subprocess
import os
import shutil
import shlex
import logging

def make_genome_file(in_fasta):
    
    try:
        c = 0 # count number of seqs
        for record in SeqIO.parse(in_fasta, "fasta"):
            mito_len = len(record)
            mito_id = record.id
            c += 1
        
        if c > 1:
            raise ValueError("More than one sequence in final_mitogenome.fasta")
    
    except ValueError:
        sys.exit(1)

    with open("final_mitogenome.genome.txt", "w") as f:
        f.write("\t".join([mito_id, str(mito_len)]))

    return "final_mitogenome.genome.txt" 

def make_genome_windows(genome_file, winSize):

    bedtools_cmd = ['bedtools', 'makewindows', '-g', genome_file, '-w', str(winSize)]
    windows_filename = genome_file.replace(".txt", f".{winSize}.txt")
    windows_file = open(windows_filename, "w")
    subprocess.run(bedtools_cmd, stdout=windows_file)

    return windows_filename

def get_windows_depth(windows_file, bam_file):

    windows_depth_filename = windows_file.replace(".txt", ".mean.depth")
    windows_depth_file = open(windows_depth_filename, "w")
    coverage_cmd = ['bedtools', 'coverage', '-a', windows_file, '-b', bam_file, '-mean']
    subprocess.run(coverage_cmd, stdout=windows_depth_file)

    return windows_depth_filename

def move_intermediate_files(files_list):
    
    if not os.path.isdir('final_mitogenome_coverage'):
        os.mkdir('final_mitogenome_coverage')
    
    for f in files_list:
        shutil.move(f, os.path.join(os.getcwd(), "final_mitogenome_coverage"))
        

def plot_coverage(depth_file, winSize):

  df = pd.read_csv(depth_file, sep="\t", names=['sequence', 'start', 'end', 'depth'])  
  df1 = df.astype({'start': 'int', 'end': 'int', 'depth': 'float'})

  df1['position'] = df1['start'] + (df1['end'] - df1['start'])/2
  df2 = df1.astype({'position': 'object'})
  
  fig, ax = plt.subplots(1,1)
  ax.bar(x=df2['position'], height=df2['depth'], width=winSize)
  ax.set_xlabel('Genome position (bp)')
  ax.set_ylabel('Coverage depth')

  plt.savefig("final_mitogenome.coverage.png")

def main():

    parser = argparse.ArgumentParser()
    parser._action_groups.pop()
    required = parser.add_argument_group('required arguments')
    optional = parser.add_argument_group('optional arguments')
    mutually_exclusive_group = optional.add_mutually_exclusive_group(required=False)
    mutually_exclusive_group.add_argument('--maxNumSeqs', help='Maximum number of sequences')
    mutually_exclusive_group.add_argument('--seqs', help='List of sequences to be used')
    required.add_argument('-ref', help='Reference file in FASTA format')
    required.add_argument('-reads', help='Reads file(s) in FASTQ format (accepts .gz compressed file)')
    optional.add_argument('-winSize', help='Windows size', default=300, type=int)
    optional.add_argument('-t', help='number of threads to use', default=1, type=int)
    
    args = parser.parse_args()

    # mapping reads to reference
    out_map_file = f"{args.ref}.hifiMapped.bam"
    minimap_cmd = ["minimap2", "-t", str(args.t), "--secondary=no", "-ax", "map-pb", args.ref] + shlex.split(args.reads) 
    samtools_cmd = ["samtools", "view", "-@", str(args.t), "-S", "-b", "-F4", "-F", "0x800"] 
    logging.info("1. Mapping HiFi reads against final_mitogenome.fasta:")
    logging.info(" ".join(minimap_cmd) + " | " + " ".join(samtools_cmd) + f" > {out_map_file}")        
    minimap = subprocess.Popen(minimap_cmd, stdout=subprocess.PIPE, stderr=subprocess.DEVNULL)
    mapped_reads_f = open(out_map_file, "w")
    subprocess.run(samtools_cmd, stderr=subprocess.STDOUT, stdin=minimap.stdout, stdout=mapped_reads_f)
    minimap.wait()
    minimap.stdout.close()


    # creating coverage plot
    logging.info("2. Creating windows file...")
    genome_filename = make_genome_file(args.ref) 
    genome_windows_filename = make_genome_windows(genome_filename, args.winSize)
    windows_depth_filename = get_windows_depth(genome_windows_filename, out_map_file) 
    plot_coverage(windows_depth_filename, args.winSize)

if __name__ == "__main__":
    main()    
