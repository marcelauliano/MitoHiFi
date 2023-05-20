#!/bin/bash

# Set base directory
BASE_DIR="."

#This script is going to run MitoHiFi, creating a folder for each run, for as many species as you have in the input file. The input file is described next.
# This is a file that has the exact structure of AdditonalFile1.tsv. It's a file with 4 columns separated by tab, where the first column is the species name, thrid column is the SRA ID of the PacBio HiFi data to be dowloaded and the forth column is the mitochondrial genetic code for that species.

#Set the input file
INPUT_FILE="AdditonalFile1-funghi.tsv"

# Read input file
while IFS=$'\t' read -r species_name _ sra_id mito_code; do
  # Extract genus and species from the species name
  genus=$(echo "$species_name" | awk '{print $1}')
  species=$(echo "$species_name" | awk '{print $2}')

  # Create directory name
  dir_name="${genus}_${species}"
  dir_name_with_spaces="${species_name//_/ }"

  # Create directory path
  dir_path="$BASE_DIR/$dir_name"

  # Create directory if it doesn't exist
  mkdir -p "$dir_path"

  # Change to the directory
  cd "$dir_path" || exit 1

  # Get the first SRA ID from the comma-separated list
  sra_id=$(echo "$sra_id" | awk -F, '{print $1}')

  # Run fastq-dump command
  fastq-dump --fasta -O . --gzip "$sra_id"

  #Get reads file
  gz_file=$(find . -maxdepth 1 -name "*.gz" -type f -printf '%f\n' -quit)

  #Run MBG
  ./MBG -k 1001 -w 250 -a 5 -u 150 -i "$gz_file" -o "$gz_file"_MBGcontigs.gfa -t 4

  #Get contigs
  awk '/^S/ {print ">"$2"\n"$3}' "$gz_file"_MBGcontigs.gfa > "$gz_file"_MBGcontigs.fa

  # Run mitohifi command
  singularity exec --bind /path/to/container_directory:/path/to/container_directory docker://ghcr.io/marcelauliano/mitohifi:master findMitoReference.py --species "$dir_name_with_spaces" --outfolder . --min_length 12000 -n 1
  
  # Get the fasta file
  fasta_file=$(ls *.fasta)

  # Get the prefix of the fasta file
  prefix="${fasta_file%.*}"

  # Find the gb file with the same prefix
  gb_file=$(find . -maxdepth 1 -name "${prefix}*.gb" -type f -print -quit)

  # Run mitohifi command - mitohifi.py
  mbg_file=$(find . -maxdepth 1 -name "*MBGcontigs.fa" -type f -printf '%f\n' -quit)
  
  # Print selected files
  echo "Selected files:"
  echo "fasta_file: $fasta_file"
  echo "gb_file: $gb_file"
  echo "gz_file: $gz_file"
  echo "mito_code: $mito_code"
  echo "mbg_file: $mbg_file"

  singularity exec --bind /path/to/container_directory:/path/to/container_directory docker://ghcr.io/marcelauliano/mitohifi:master mitohifi.py -c "$mbg_file" -f "$fasta_file" -g "$gb_file" -o "$mito_code" -t 4

  # Change back to the previous directory
  cd ..

done < "$INPUT_FILE"
