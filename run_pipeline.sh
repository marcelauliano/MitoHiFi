fasta=$1
mito_fa=$2
mito_gb=$3
threads=$4

echo -e "\nFirst let's run the blast with the close-related mitogenome"

blast/bin/makeblastdb -in $mito_fa -dbtype nucl
echo -e "\nmakeblastdb done. Running blast with CCS contigs\n"

blast/bin/blastn -query $fasta -db $mito_fa -num_threads $threads -out $fasta.blastn -outfmt '6 std qlen slen'
echo -e "blast done!\n"

cat $fasta.blastn | awk '{print$1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$11"\t"$12"\t"$13"\t"$14"\t"(($4*100)/$13)"\t"(($4*100)/$14)}' > $fasta.blastn.cov
#exclude contigs smaller than close-related species mitogenome size
cat $fasta.blastn.cov | awk '{ if ($13 >= $14) print $0}' > $fasta.blastn.cov.NoPartials
echo -e "\nNow we check which contigs have 97% or more of its length in the blast match and print ids to a file\n"
awk '{arr[$1]+=$15} END {for (i in arr) {print i,arr[i]}}' $fasta.blastn.cov.NoPartials > $fasta.blastn.cov.NoPartials.sum
cat $fasta.blastn.cov.NoPartials.sum | awk '{if ($2 >= 97) print $0}' > $fasta.blastn.cov.NoPartials.sum.97
cat $fasta.blastn.cov.NoPartials.sum.97 | awk '{print $1}' > $fasta.blastn.cov.NoPartials.sum.97.id
echo -e "Those ones\n"
cat $fasta.blastn.cov.NoPartials.sum.97.id
echo -e "\nAs we have more than one, let's go with the largest contig and let's circularize it so we can be sure our mitogenone is complete\n"
python scripts/filterfasta.py -i $fasta.blastn.cov.NoPartials.sum.97.id $fasta > $fasta.blastn.cov.NoPartials.sum.97.fasta
python scripts/get_Larger.py $fasta.blastn.cov.NoPartials.sum.97.fasta > $fasta.blastn.cov.NoPartials.sum.97.LargerContig.fasta.id
python scripts/filterfasta.py -i $fasta.blastn.cov.NoPartials.sum.97.LargerContig.fasta.id $fasta.blastn.cov.NoPartials.sum.97.fasta > $fasta.blastn.cov.NoPartials.sum.97.LargerContig.fasta.fasta
python scripts/circularizationCheck.original.py $fasta.blastn.cov.NoPartials.sum.97.LargerContig.fasta.fasta
python scripts/cut_coords.py $fasta.blastn.cov.NoPartials.sum.97.LargerContig.fasta.fasta > mito.circu.fasta
/software/team311/mu2/MitoFinder/mitofinder -j mito.circu.ANNOTATION -a mito.circu.fasta -r $mito_gb -o 5
