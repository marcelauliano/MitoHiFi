import sys
import subprocess

def make_genome_file(in_bam, target_seq):
    samtools_cmd = ['samtools', 'view', '-h', in_bam]
    proc = subprocess.Popen(samtools_cmd, stdout=subprocess.PIPE)
    for bytes_line in proc.stdout:
        line = bytes_line.decode("utf-8")
        if line.split()[0] == "@HD":
            pass
        elif line.split()[0] == "@SQ":
            seq_ID = line.split()[1].replace("SN:","")
            seq_len = line.split()[2].replace("LN:", "")
            if seq_ID == target_seq:
                with open(f"{target_seq}.genome.txt", "w") as f:
                    f.write("\t".join([seq_ID, seq_len]) + "\n")
                    break
    return (f"{target_seq}.genome.txt", seq_len)

def make_genome_windows(genome_file, win_size):
    bedtools_cmd = ['bedtools', 'makewindows', '-g', genome_file, '-w', str(win_size)]
    windows_filename = genome_file.replace(".txt", f"{win_size}.txt")
    print(f"Creating windows genome file and saving it as {windows_filename}")
    windows_file = open(windows_filename, "w")
    subprocess.run(bedtools_cmd, stdout=windows_file)

    return windows_filename
    
if __name__ == "__main__":
    in_file = sys.argv[1]
    target_seq = sys.argv[2]
    make_genome_file(in_file, target_seq)
