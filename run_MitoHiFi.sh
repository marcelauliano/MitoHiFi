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

# This awk line calculates the percentage of the query and subject in the blast match on each line. Percetages become columns 15 (query) and 16 (subject)
#cat ${contigs}.blastn | awk {print$($4*100)/$14)}' > ${contigs}.blastn.cov

#awk bellow sums up the percentages of each query on a blast hit. Column 2 is the percentage, column 3 the is the query size and column 4 the subject size.
#echo -e \nNow we check which contigs have 97% or more of its length in the blast match with the close-related species mitogenome:\n"
#awk {arr[$1]+=$15} END {for (i in arr) {print i,arr[i]}}' ${contigs}.blastn.cov > ${contigs}.blastn.cov.sum

#awk bellow calculates the percentage of the size of the query in relation to the close-related mitogenome. If your contigs are less than 80% of the size of the close-related mito, its likely you don't have the complete mitogenome
#cat ${contigs}.blastn.cov.sum | awk #{print$1(($3*100)/$4)}' > ${contigs}.blastn.cov.sum.perc
#percentage=$(cat ${contigs}.blastn.cov.sum.perc | awk {print$5}' )

#if [[ $percentage < 80 ]]; then
#	echo -e \nIt seems that all of your contigs are 80% smaller than the close-related mitogenome size. Not sure if you have the complete mitogenome assembled here! Exiting...\n'
#	exit 0
#else
#	echo -e \n...\n'
#fi
#get the IDs for the queries that are at least 88% on the blast match. If more than 97% of the query is on the blast hit, its likely its the mitogenome.
#cat ${contigs}.blastn.cov.sum | awk '{if > ${contigs}.blastn.cov.sum.perc.88
#cat ${contigs}.blastn.cov.sum.perc.88 | awk  ${contigs}.blastn.cov.sum.perc.88.id

#echo -e se ones are contig that have 88% or more of their total length in a blast match with the close related mitogenome. Their are likely to  represent the mitogenome for your species!\n"
#cat ${contigs}.blastn.cov.sum.perc.88.id


#echo -e Let's circularize the largest so we can be sure our mitogenone is complete\n"


#get the fasta sequences
python vamos2.py contigs.blastn
python scripts/filterfasta.py -i potentialmito.id ${contigs} > ${contigs}.mito.fa

#sort the fasta by the larger length and save it's id to a file
#python scripts/get_Larger.py ${contigs}.blastn.cov.sum.perc.88.fasta > ${contigs}.blastn.cov.sum.perc.88.Largest.id

#get the fasta for the largest contig. Rembering that all contigs worked so far have 97% or more of its length in a blast match with the close-related mitogenome. So it's unlikely they would represent NUMTs. 
#python scripts/filterfasta.py -i ${contigs}.blastn.cov.sum.perc.88.Largest.id ${contigs}.blastn.cov.sum.perc.88.fasta > ${contigs}.LargerContig.fasta

#find the coordinates where mitogenome cirularises 
python cur.py ${contigs}.mito.fa

#cut the fasta to get only one copy of the mitogenome
python scripts/cut_coords.py ${contigs}.mito.fa  > mitogenomeehehehehe.fasta

#annotate the mitogenome with mitofinder
mitofinder -j mitogenome.annotation -a mitogenomeehehehehe.fasta -r ${genbank} -o ${mitocode}

echo -e "\nPipeline done!!!\n Your mito genome is the file mitogenome.fasta. \n Annotation: Please look inside the mitofinder Final_Result folder to find your mitogenome annotated in genbank format.\n"
printf "\n\nDone!" &&
printf "\n\nCompleted at: $(date "+%Y-%m-%d %H-%M-%S")\n\n"
