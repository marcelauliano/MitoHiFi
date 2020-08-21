#!/bin/bash

set -e -o pipefail

#++++                  This script is part of:                    ++++
#++++                  CCS mitogenome                             ++++
#++++                  Darwin Tree of Life Assembly Pipeline      ++++
#++++                  Credit: M Uliano-Silva                     ++++


if [ -z $1 ]; then

	echo "use $0 -h for help"
	exit 0
elif [ $1 == "-h" ]; then

cat<<EOF
	 Usage: '$0 -c contigs.fasta -f close-related_mitogenome.fasta -g close-related_mitogenome.gb -t threads'
	-c: assemnbled fasta contigs/scaffolds to be searched to find mitogenome
	-f: Close-related species mitogenome in fasta format
	-g: Close-related species mitogenome in genbank format 
	-t: Number of threads for the blast search 
	-o: <integer> MitoFinder parameter: Organism genetic code following NCBI table (integer):
                        1. The Standard Code 2. The Vertebrate Mitochondrial
                        Code 3. The Yeast Mitochondrial Code 4. The Mold,
                        Protozoan, and Coelenterate Mitochondrial Code and the
                        Mycoplasma/Spiroplasma Code 5. The Invertebrate
                        Mitochondrial Code 6. The Ciliate, Dasycladacean and
                        Hexamita Nuclear Code 9. The Echinoderm and Flatworm
                        Mitochondrial Code 10. The Euplotid Nuclear Code 11.
                        The Bacterial, Archaeal and Plant Plastid Code 12. The
                        Alternative Yeast Nuclear Code 13. The Ascidian
                        Mitochondrial Code 14. The Alternative Flatworm
                        Mitochondrial Code 16. Chlorophycean Mitochondrial
                        Code 21. Trematode Mitochondrial Code 22. Scenedesmus
                        obliquus Mitochondrial Code 23. Thraustochytrium
                        Mitochondrial Code 24. Pterobranchia Mitochondrial
                        Code 25. Candidate Division SR1 and Gracilibacteria
EOF
	exit 0
fi

#set options

while getopts ":c:f:g:t:o:" opt; do

	case $opt in
	
	c)
		contigs="$OPTARG"
		;;
	f)
		fasta="$OPTARG"
		;;
	g)
		genbank="$OPTARG"
		;;
	t)
		threads="$OPTARG"
		;;

	o)	mitocode="$OPTARG"
		;;
	\?)
		echo ""Usage: sh runMitoHiFi.sh [-c] [-f] [-g] [-t] [-o]
		exit 0
		;;
	esac
done

printf "\n\n++++                        mitoHiFi beta version                       ++++\n"
printf     "++++ Darwin Tree of Life Hifi mitogenome circularisation and annotation ++++\n"
printf     "++++                      Credit: M Uliano-Silva                        ++++\n\n"

printf "\nStarted at at: $(date "+%Y-%m-%d %H-%M-%S")\n"


echo -e "\nFirst let's run the blast with the close-related mitogenome\n"

makeblastdb -in ${fasta} -dbtype nucl
echo -e "\nmakeblastdb done. Running blast with CCS contigs\n"

blastn -query ${contigs} -db ${fasta} -num_threads ${threads} -out contigs.blastn -outfmt '6 std qlen slen'
echo -e "Blast done!\n"


python scripts/parse_blast.py 

python scripts/filterfasta.py -i contig.id ${contigs} > ${contigs}.mito.fa

python script/circularizationCheck.modified.py ${contigs}.mito.fa

#cut the fasta to get only one copy of the mitogenome
python scripts/cut_coords.py ${contigs}.mito.fa  > mitogenome.fasta

#annotate the mitogenome with mitofinder
mitofinder -j mitogenome.annotation -a mitogenome.fasta -r ${genbank} -o ${mitocode}

echo -e "\nPipeline done!!!\n Your mito genome is the file mitogenome.fasta. \n Annotation: Please look inside 'mitogenome.annotation/mitogenome.annotation_Final_Results/' folder to find your mitogenome annotated in genbank format.\n"
printf "\n\nDone!" &&
printf "\n\nCompleted at: $(date "+%Y-%m-%d %H-%M-%S")\n\n"
