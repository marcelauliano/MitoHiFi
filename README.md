# Mito HiFi 

------ This is a beta version! -------


This repository contains software and scripts to circularise, cut and annotate mitogenomes assembled redudantly by assembly softwares such as HiCanu or Hifiasm.

The dissemination of PacBio HiFi reads makes the assemble of high-quality mitogenome straigh forward. Because of the circular nature of the molecule, however, the mitocontig is usually assembled rudandantly resulting in a contiguos multiple-copies mitogenome (IMPROVE THIS DESCRIPTION!!). This pipeline was developted to finalize the assembly and annotation of the mitogenome by (i) finding it among total assembled nuclear contigs, (ii) circularising and cutting it to represent only one copy of the circular molecule and (iii) annotation and presenting it in fasta and genbank format.


<b>Dependencies:</b>

- Blast (makeblastdb and blastn) have to be installed and export to your PATH
- mitoFinder: has to be installed [mitoFinder](https://github.com/RemiAllio/MitoFinder) and export to your PATH 
- Biopython
- Python Pandas

<b>Installation</b>

### Get and install MitoHiFi (Linux)

```

git clone https://github.com/marcelauliano/MitoHiFi.git

```

### Test run

```
cd MitoHiFi
cd exampleFiles
ln -s ../scripts
ln -s ../run_MitoHiFi.sh
sh run_MitoHiFi.sh -c test.fa -f NC_016067.1.fasta -g NC_016067.1.gb -t 1 -o 5

```
### Required arguments

```
Usage: 'sh run_mitoHiFi.sh -c contigs.fasta -f close-related_mitogenome.fasta -g close-related_mitogenome.gb -t threads'
Parameters:	
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
 ```
 
 ### Output
 
 ```
 mitogenome.fasta  - this is your final mitogenome in fasta format
 
Inside 'mitogenome.annotation/mitogenome.annotation_Final_Results/' you find mitogenome.annotation_mtDNA_contig.gb, which is the annotation of your mitogenome performed by mitofinder and outputed in genbank format. 
 
```
### Further

If you would like to rotate your mitogenome to start at tRNA-Phe, check the coordinates of it on your genbank file and run:
```
python scripts/rotate.py -i mitogenome.fasta -r <coordinate> > mitogenome.rotated.fa
```

 ### Description intermediate outputs
 
<b>parsed_blast.txt:</b>   - tab separated file with 4 columns as follows


 - qseqid - the ID of your input contigs
 - %q_in_match - a percentage of the length of your contig in a blast match with the close related species  
 - leng_query - length of your contig
 - s_length  lenght of the close-related mitogenome given
 
 <b>contigs.blastn:</b> - outfm 6 regular tab output plus 2 extra columns which contain length_of_query and length_of_subject 
 

 

 
 For more information mu2@sanger.ac.uk
