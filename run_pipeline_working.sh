#!/bin/bash

set -e -o pipefail

fasta=$1
mito_fa=$2
mito_gb=$3
threads=$4

#++++                  This script is part of:                    ++++
#++++                       CCS mitogenome                        ++++
#++++                  Darwin Tree of Life Assembly Pipeline      ++++
#++++                   Credit: M Uliano-Silva                    ++++


if [ -z $1 ]; then

	echo "use $0 -h for help"
	exit 0
elif [ $1 == "-h" ]; then

        cat << EOF
        
        Usage: <contigs.fasta>  <close-related_mitogenome.fasta> <close-related mitogenome.gb> <num_threads>"
        "<contigs.fasta>  fasta contigs to search for mitogenome."     
        "<close-related_mitogenome.fasta> Close-related species mitogenome in fasta format"
         "<close-related mitogenome.gb> Close-related species mitogenome in genbank format" 
        "<num_threads> Number of threads for the blast search" 
EOF 

        exit 0
fi

if [ -z $2 ]; then
        echo "No close-related species fasta provided"
        exit 1
else
        echo "Close-related mitogenome fasta is $mito_fa"
fi

if [ -z $3 ]; then
        echo "No close-related species genbank file provided"
        exit 1
else
        echo "Close-related mitogenome genbank is $mito_gb"
fi

echo -e "\nFirst let's run the blast with the close-related mitogenome\n"

blast/bin/makeblastdb -in $mito_fa -dbtype nucl
echo -e "\nmakeblastdb done. Running blast with CCS contigs\n"

blast/bin/blastn -query $fasta -db $mito_fa -num_threads $threads -out $fasta.blastn -outfmt '6 std qlen slen'
echo -e "blast done!\n"

# This awk line bellow calculates how much percentage of the query and subject are in the blast match on each line
cat $fasta.blastn | awk '{print$1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$11"\t"$12"\t"$13"\t"$14"\t"(($4*100)/$13)"\t"(($4*100)/$14)}' > $fasta.blastn.cov

#line bellow excludes contigs smaller than close-related species mitogenome size
cat $fasta.blastn.cov | awk '{ if ($13 >= $14) print $0}' > $fasta.blastn.cov.NoPartials

#awk bellow sums up the percentages of each query on a blast hit. 
echo -e "\nNow we check which contigs have 97% or more of its length in the blast match and print ids to a file\n"
awk '{arr[$1]+=$15} END {for (i in arr) {print i,arr[i]}}' $fasta.blastn.cov.NoPartials > $fasta.blastn.cov.NoPartials.sum

#If more than 97% of the query is on the blast hit, its likely its the mitogenome.
cat $fasta.blastn.cov.NoPartials.sum | awk '{if ($2 >= 97) print $0}' > $fasta.blastn.cov.NoPartials.sum.97

#get the IDs for the queries that are more than 97% on the blast match
cat $fasta.blastn.cov.NoPartials.sum.97 | awk '{print $1}' > $fasta.blastn.cov.NoPartials.sum.97.id
echo -e "Those ones\n"
cat $fasta.blastn.cov.NoPartials.sum.97.id


echo -e "\nAs we have more than one, let's go with the largest contig and let's circularize it so we can be sure our mitogenone is complete\n"

#get the fasta sequences
python scripts/filterfasta.py -i $fasta.blastn.cov.NoPartials.sum.97.id $fasta > $fasta.blastn.cov.NoPartials.sum.97.fasta

#sort the fasta by the larger length and save it's id to a file
python scripts/get_Larger.py $fasta.blastn.cov.NoPartials.sum.97.fasta > $fasta.blastn.cov.NoPartials.sum.97.LargerContig.fasta.id

#get the ultimate fasta
python scripts/filterfasta.py -i $fasta.blastn.cov.NoPartials.sum.97.LargerContig.fasta.id $fasta.blastn.cov.NoPartials.sum.97.fasta > $fasta.blastn.cov.NoPartials.sum.97.LargerContig.fasta.fasta

#find the coordinates where mitogenome cirularises 
python scripts/circularizationCheck.original.py $fasta.blastn.cov.NoPartials.sum.97.LargerContig.fasta.fasta

#cut the fasta to get only one copy of the mitogenome
python scripts/cut_coords.py $fasta.blastn.cov.NoPartials.sum.97.LargerContig.fasta.fasta > mito.circu.fasta

#annotate the mitogenome with mitofinder
scripts/MitoFinder/mitofinder -j mito.circu.ANNOTATION -a mito.circu.fasta -r $mito_gb -o 5

echo "your mito genome is the file mito.circu.fasta. Please look inside the mitofinder Final Result folder to find your mitogenome annotated in genbank format."
