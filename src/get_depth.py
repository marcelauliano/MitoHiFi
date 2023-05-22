import subprocess

def get_depth(bam_file, genome_file, target_seq):
    #depth_filename = in_bam.replace('.bam', '.depth')
    depth_filename = f"{target_seq}.perBase.depth"
    depth_file = open(depth_filename, 'w')
    print(f'Calculating depth per base for sequence {target_seq} and saving it to {depth_filename}')
    coverage_cmd = ['bedtools', 'coverage', '-a', genome_file, '-b', bam_file ,'-d']
    subprocess.run(coverage_cmd, stdout=depth_file)

def get_windows_depth(windows_file, bam_file, target_seq):
    windows_depth_filename = windows_file.replace(".txt", ".mean.depth")
    print(f"Calculating depth per base for sequence {target_seq} and saving it to {windows_depth_filename}")
    windows_depth_file = open(windows_depth_filename, "w")
    coverage_cmd = ['bedtools', 'coverage', '-a', windows_file, '-b', bam_file, '-mean']
    subprocess.run(coverage_cmd, stdout=windows_depth_file)

    return windows_depth_filename