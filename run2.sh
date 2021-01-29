#!/bin/bash

set -e -o pipefail

#++++                  This script is part of:                    ++++
#++++                  CCS mitogenome                             ++++
#++++                  Darwin Tree of Life Assembly Pipeline      ++++
#++++                  Credit: Marcela Uliano-Silva               ++++


if [ -z $1 ]; then

	echo "use $0 -h for help"
	exit 0
elif [ $1 == "-h" ]; then

cat<<EOF
	 Usage: '$0 -c contigs.fasta -f close-related_mitogenome.fasta -g close-related_mitogenome.gb -t threads'
        -r: Pacbio HiFi reads from your species
	-f: Close-related species mitogenome in fasta format
	-g: Close-related species mitogenome in genbank format 
	-t: Number of threads for (i) hifiasm and (ii) the blast search 
        -l: Length of the close-related mitogenome (to filter out possible NUMTs)
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

while getopts ":r:c:f:g:t:l:o:" opt; do

	case $opt in
	
	r)
		reads="$OPTARG"
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
	l)
		length="$OPTARG"
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

echo -e "\nFirst we map your PacbioHiFi reads to the close-related mitogenome\n"

minimap2 -t ${threads} --secondary=no -ax map-pb ${genbank} ${reads} | samtools view -@10 -S -b -F4 -F 0x800 > ${genbank}.HiFiMapped.bam

echo -e "\nNow we filter out any mapped reads that are larger than the reference mitogenome to avoid NUMTS.\n"

#Get the mapped reads in fasta format
samtools fasta ${genbank}.HiFiMapped.bam > ${genbank}.HiFiMapped.bam.fasta

#Filter out any mapped reads that are larger than the reference mitogenome to avoid NUMTs
python scripts/filterfasta.py -l ${length} -n ${genbank}.HiFiMapped.bam.fasta > ${genbank}.HiFiMapped.bam.filtered.fasta

echo -e "\nNow let's run hifiasm to assemble the mapped and filtered reads!.\n"

hifiasm -t${threads} -o ${genbank}.HiFiMapped.bam.filtered.assembled ${genbank}.HiFiMapped.bam.filtered.fasta 2>hifiasm.log

scripts/gfa2fa ${genbank}.HiFiMapped.bam.filtered.assembled.p_ctg.gfa > ${genbank}.HiFiMapped.bam.filtered.assembled.p_ctg.fa
scripts/gfa2fa ${genbank}.HiFiMapped.bam.filtered.assembled.a_ctg.gfa > ${genbank}.HiFiMapped.bam.filtered.assembled.a_ctg.fa
cat ${genbank}.HiFiMapped.bam.filtered.assembled.p_ctg.fa ${genbank}.HiFiMapped.bam.filtered.assembled.a_ctg.fa > hifiasm.contigs.fasta

echo -e "\nNow let's run the blast of the assembled contigs with the close-related mitogenome\n"

makeblastdb -in ${fasta} -dbtype nucl
echo -e "\nmakeblastdb done. Running blast with the assembled hifiasm contigs\n"

blastn -query hifiasm.contigs.fasta -db ${fasta} -num_threads ${threads} -out contigs.blastn -outfmt '6 std qlen slen'
echo -e "Blast done!\n"

#the next script parses a series of conditions to exclude blast with NUMTs. 
python scripts/parse_blast.py 

#Next, we extract the mitogenome contig
python scripts/filterfasta.py -i contig.id ${contigs} > ${contigs}.mito.fa

#We check for circularisation
python scripts/circularizationCheck.modified.py ${contigs}.mito.fa

#If it circularises, we cut the fasta to get only one copy of the mitogenome
python scripts/cut_coords.py ${contigs}.mito.fa  > mitogenome.fasta

#annotate the mitogenome with mitofinder
mitofinder -j mitogenome.annotation -a mitogenome.fasta -r ${genbank} -o ${mitocode}
rm *.nsq *.nin *.nhr *.xml
echo -e "\nPipeline done!!!\n Your mito genome is the file mitogenome.fasta. \n Annotation: Please look inside 'mitogenome.annotation/mitogenome.annotation_Final_Results/' folder to find your mitogenome annotated in genbank format.\n"
printf "\n\nDone!" &&
printf "\n\nCompleted at: $(date "+%Y-%m-%d %H-%M-%S")\n\n"
