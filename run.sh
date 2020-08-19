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
EOF
	exit 0
fi

#set options

while getopts ":c:f:g:t:" opt; do

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
	\?)
		echo ""Usage: cmd [-c] [-f] [-g] [-t]
		;;
	esac
done

printf "\n\n++++                        mitoHiFi beta version            ++++\n"
printf "++++ Darwin Tree of Life Hifi mitogenome circularisation and annotation ++++\n"
printf "++++     Credit: M Uliano-Silva       ++++\n\n"

printf "\nStarted at at: $(date "+%Y-%m-%d %H-%M-%S")\n"
printf "\nWith command:\n"
printf "\n${0}\n"
echo -e "\nFirst let's run the blast with the close-related mitogenome\n"

makeblastdb -in ${fasta} -dbtype nucl
echo -e "\nmakeblastdb done. Running blast with CCS contigs\n"

blastn -query ${contigs} -db ${fasta} -num_threads ${threads} -out ${contigs}.blastn -outfmt '6 std qlen slen'
echo -e "blast done!\n"

# This awk line bellow calculates how much percentage of the query and subject are in the blast match on each line
cat ${contigs}.blastn | awk '{print$1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$11"\t"$12"\t"$13"\t"$14"\t"(($4*100)/$13)"\t"(($4*100)/$14)}' > ${contigs}.blastn.cov

#line bellow excludes contigs smaller than close-related species mitogenome size
cat ${contigs}.blastn.cov | awk '{ if ($13 >= $14) print $0}' > ${contigs}.blastn.cov.NoPartials

#awk bellow sums up the percentages of each query on a blast hit. 
echo -e "\nNow we check which contigs have 97% or more of its length in the blast match and print ids to a file\n"
awk '{arr[$1]+=$15} END {for (i in arr) {print i,arr[i]}}' ${contigs}.blastn.cov.NoPartials > ${contigs}.blastn.cov.NoPartials.sum

#If more than 97% of the query is on the blast hit, its likely its the mitogenome.
cat ${contigs}.blastn.cov.NoPartials.sum | awk '{if ($2 >= 97) print $0}' > ${contigs}.blastn.cov.NoPartials.sum.97

#get the IDs for the queries that are more than 97% on the blast match
cat ${contigs}.blastn.cov.NoPartials.sum.97 | awk '{print $1}' > ${contigs}.blastn.cov.NoPartials.sum.97.id
echo -e "Those ones\n"
cat ${contigs}.blastn.cov.NoPartials.sum.97.id


echo -e "\nAs we have more than one, let's go with the largest contig and let's circularize it so we can be sure our mitogenone is complete\n"

#get the fasta sequences
python scripts/filterfasta.py -i ${contigs}.blastn.cov.NoPartials.sum.97.id ${contigs} > ${contigs}.blastn.cov.NoPartials.sum.97.fasta

#sort the fasta by the larger length and save it's id to a file
python scripts/get_Larger.py ${contigs}.blastn.cov.NoPartials.sum.97.fasta > ${contigs}.blastn.cov.NoPartials.sum.97.LargerContig.fasta.id

#get the ultimate fasta
python scripts/filterfasta.py -i ${contigs}.blastn.cov.NoPartials.sum.97.LargerContig.fasta.id ${contigs}.blastn.cov.NoPartials.sum.97.fasta > ${contigs}.LargerContig.fasta

#find the coordinates where mitogenome cirularises 
python scripts/circularizationCheck.original.py ${contigs}.LargerContig.fasta

#cut the fasta to get only one copy of the mitogenome
python scripts/cut_coords.py ${contigs}.LargerContig.fasta > mitogenome.fasta

#annotate the mitogenome with mitofinder
mitofinder -j mitogenome.annotation -a mitogenome.fasta -r ${genbank} -o 5

echo -e "\nPipeline done!!!\n Your mito genome is the file mitogenome.fasta. \n Annotation: Please look inside the mitofinder Final_Result folder to find your mitogenome annotated in genbank format.\n"
printf "\n\nDone!" &&
printf "\n\nCompleted at: $(date "+%Y-%m-%d %H-%M-%S")\n\n"
