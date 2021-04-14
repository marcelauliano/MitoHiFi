import os 
import shutil
import sys

def clean_up_work_dir(contigs_list): 
    # create directory that will keep all contigs selected after alignment with related mito
    if not os.path.isdir("potential_contigs"):
        os.mkdir("potential_contigs")
    # iterate over each contig, searching for files that have it as prefix (that should be moved
    # to the related contig folder (contig_dir_path)
    for contig in contigs_list:
        for curr_file in os.listdir('.'):
            if curr_file.endswith((".nhr", ".nin", "nsq")): # delete intermediate files
                curr_file_path = os.path.join(os.getcwd(), curr_file)
                os.remove(curr_file_path)
                continue
            # move mapping and assembly files (generated in read mode) to reads_mapping_and_assembly folder
            if curr_file.startswith("gbk.HiFiMapped") or curr_file=="reads.HiFiMapped.bam":
                if not os.path.isdir("reads_mapping_and_assembly"):
                    os.mkdir("reads_mapping_and_assembly")
                reads_mapping_path = os.path.join(os.getcwd(), "reads_mapping_and_assembly")
                shutil.move(curr_file, reads_mapping_path)
                continue
            if curr_file.startswith(contig): # move files to their corresponding contigs folders
                contig_dir_path = os.path.join(os.getcwd(), "potential_contigs", contig)
                if not os.path.isdir(contig_dir_path):
                    os.mkdir(contig_dir_path)
                shutil.move(curr_file, contig_dir_path)
    
    # moving all files that contain info related to circularization step to `contigs_circularization` folder
    circ_files = ["all_contigs.circularisationCheck.txt", "circularization_check.blast.tsv", "circularization_check.blast.xml"]
    if not os.path.isdir("contigs_circularization"):
        os.mkdir("contigs_circularization")
    for f in circ_files:
        if os.path.isfile(f):
            shutil.move(f, "contigs_circularization")

    # moving files related to the process of choosing the representative final mitogenome (using CDHIT)
    final_choice_files = ["all_mitogenomes.rotated.aligned.fa", "all_mitogenomes.rotated.fa", "cdhit.out", "cdhit.out.clstr"]
    if not os.path.isdir("final_mitogenome_choice"):
        os.mkdir("final_mitogenome_choice")
    for f in final_choice_files:
        if os.path.isfile(f):
            shutil.move(f, "final_mitogenome_choice")

    # if both hifiasm files exist, move them to the reads_mapping_and_assembly folder
    if all(x in os.listdir('.') for x in ["hifiasm.log", "hifiasm.contigs.fasta"]):
        if not os.path.isdir("reads_mapping_and_assembly"):
            os.mkdir("reads_mapping_and_assembly")
        reads_mapping_path = os.path.join(os.getcwd(), "reads_mapping_and_assembly")
        shutil.move("hifiasm.log", reads_mapping_path)
        shutil.move("hifiasm.contigs.fasta", reads_mapping_path)
    
    # move files from contigs filtering step (after BLASTing against related mito) to contigs_filtering folder
    contigs_selection_files = ["contigs.blastn", "contigs_ids.txt", "parsed_blast_all.txt", "parsed_blast.txt"]
    contigs_selection_folder = "contigs_filtering"
    if not os.path.isdir(contigs_selection_folder):
        os.mkdir(contigs_selection_folder)
    for f in contigs_selection_files:
        if os.path.isfile(f):
            shutil.move(f, os.path.join(os.getcwd(), contigs_selection_folder))

def main():
    contigs_list = []
    for curr_file in os.listdir('.'):
        if curr_file.endswith(".mitogenome.rotated.fa"):
            contig_id = curr_file.replace(".mitogenome.rotated.fa", "")
            contigs_list.append(contig_id)
    clean_up_work_dir(contigs_list)

if __name__ == "__main__":
    main()
