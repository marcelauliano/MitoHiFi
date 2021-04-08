# MitoHiFi 

------ This is v2 -------

MitoHiFi is a python pipeline distributed under the [license](https://github.com/marcelauliano/MitoHiFi/blob/master/scripts/LICENSE) 

--------------------------------------


<b>MitoHiFi</b> assembles a species mitogenome from Pacbio HiFi reads.

With Mitoifi.v2 you can start from raw Pacbio HiFi reads (flag -r) or from assembled contigs (flag -c). You will also going to need a close-related mitochondria sequence in fasta and gb format. We have an internal script that can download it for you from NCBI (findMitoReference.py).



*The dissemination of high-quality long reads - such as PacBio HiFi - makes the assembly of high-quality mitogenome straight forward. Because of the circular nature of the molecule, however, the mitocontig is usually assembled redundantly resulting in multiple-copy mitogenome-contigs. This pipeline was developed to finalise the assembly and annotation of the mitogenome. It will also dected different variants of the mitogenome present in your sample. At the end you are going to have all the variants assembled and annotated, and MitoHiFi.v2 is going to choose a final consensus sequence. In addtion, you will find an aligment of all the variants to facilitate your analysis of mitochondria heteroplasmy.*

MitoHifi v2 will:

(i) extract mito reads and assemble them with hifiasm (flag -r), or find the mito contigs among assembled contigs (flag -c)    
(ii) indetify and separate NUMTS from real mitochontigs  
(iii) generate a circularized, non-redudant and annotated version of all the mitochondria sequences present in your sample
(iv) choose a final consensus as the final mitochondria

-----
-----

### Installation

There are two ways to install MitoHifi.v2 at the moment; (i) mannually - and then you will need to have all the dependencies installed in your PATH, or with (ii) a singularity image.

### Manual installation - Dependencies

- BLAST+ (makeblastdb and blastn) have to be installed and export to your PATH
- MitoFinder: has to be installed [MitoFinder](https://github.com/RemiAllio/MitoFinder) and export to your PATH 
- Biopython
- Python Pandas
- Mafft
- cdhit
- hifiasm
- samtools
- minimap2

<b>Installation</b>

### Once dependencies are done, install MitoHiFi.v2 (Linux)

```

git clone https://github.com/marcelauliano/MitoHiFi.git

```

### Running MitoHiFi.v2 from a Singularity image

We recommned using singularity versions => 3.7, as lower versions do not support spaces in the arguments, and you would not be able to pass more than one set of reads to -r

MitoHiFi.v2 was wrapped up into singularity container which you can use as:

singularity exec --bind /path/on/disk/to/data/:/data/ /path/to/mitohifi-v3.sif  mitohifi-v3-fromCirc_v02.11.3.py -r "/data/f1.fasta /data/f2.fasta /data/f3.fasta" -f /data/reference.fasta -g /data/reference.gb  -t 20 -o 2

Singluarity versions lower than 3.7 do not support spaces in the arguments, so if you want to pass several read datasets as in the example above use singularity version 3.7 or higher. 

The script for creating the reference files is incorporated into singularity image and can be used as follows:

singularity exec --bind /path/on/disk/to/data/:/data/ /path/to/mitohifi-v3.sif  findMitoReference.py --species "Cryptosula pallasiana" --email your@email.for.ncbi.db.query --outfolder /data/ --min_length 16000 


### Running MitoHifi_v2 with test data

```
cd MitoHiFi
cd exampleFiles
ln -s ../scripts
ln -s ../run_MitoHiFi.sh
sh run_MitoHiFi.sh -c test.fa -f NC_016067.1.fasta -g NC_016067.1.gb -t 1 -o 5

```
### Required arguments

1-) To run this pipeline, first you need a close-related mitochondria in fasta and genbank format. We have a script that can help you find this input. Giving the name of the species you are assembling, the script is going to look for the closest mitochondria it can find on NCBI. You can give the parameter -s to the script if you would like to download a partial mitochondria, but only for a species of the same genus. Otherwise, without -s, the script is going to search for complete mitochondrias only and as close as possible to your species on interest.

Using the species in our test data as an example, you would do:

findMitoReference.py --species "Cryptosula pallasiana" --email your@email.for.ncbi.db.query --outfolder /data/ --min_length 16000

This will output you xxx.fasta and xxx.gb that you will use as flags -f and -g in the main pipeline.

2-) Now, you need to decide if you are running MitoHiFi.v2 from (i) raw reads, in which case the pipeline is going to mapp your reads to the close-related species and then assemble the with Hifiasm, or (ii) your already have a Pacbio HiFi assembly and you are going to give the contigs to the pipeline (flag -c).

2.1-) If you are starting from raw reads, your required commands are:

```
Usage: 'python mitohifi_v2.py -r "f1.fasta f2.fasta f3.fasta" -f reference.fasta -g reference.gb  -t <int> -o <int> '

Parameters descriptions:
	-r: PacBio HiFi reads
	-f: Close-related species mitogenome in fasta format
	-g: Close-related species mitogenome in genbank format 
	-t: Number of threads for minimap2, hifiasm and the blast search 
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

2.2-) If you are starting from assembled contigs, your required commands are:

```
Usage: 'python mitohifi_v2.py -c contigs.fasta -f reference.fasta -g reference.gb  -t <int> -o <int> '

Parameters descriptions:
	-c: contigs # from assemblers such as Hicanu or Hifiasm
	-f: Close-related species mitogenome in fasta format
	-g: Close-related species mitogenome in genbank format 
	-t: Number of threads for minimap2, hifiasm and the blast search 
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



- 
- To run this pipeline you need 3 inputs: (i) your multifasta contig files, and a close-related species mitochondrial genome in (ii) fasta and in (iii) genbank format.

```
Usage: 'sh run_mitoHiFi.sh -c contigs.fasta -f close-related_mitogenome.fasta -g close-related_mitogenome.gb -t threads -o <integer>'
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
If you would like to finalise another contig as the mitogenome - for studies of heteroplasmy for example - check the output <b>parsed_blast.txt</b> , choose the contig ID, use the scripts/filterfasta.py to extract only that contig and run the pipeline with it again.

```
python scripts/filterfasta.py -i contig.id <your_initial_input> > contig.id.fa

sh run_MitoHiFi.sh -c contig.id.fa -f <close-related-mito>.fasta -g <close-related-mito>.gb -t <int> -o <int>
```

 ### Description intermediate outputs
 
<b>parsed_blast.txt</b>   - tab separated file with 4 columns as follows


 - qseqid - the ID of your input contigs
 - %q_in_match - a percentage of the length of your contig in a blast match with the close related species  
 - leng_query - length of your contig
 - s_length  lenght of the close-related mitogenome given
 
 <b>contigs.blastn</b> - outfmt 6 blast output plus 2 extra columns containing respectively length_of_query and length_of_subject 
 
<b>circularisationCheck.txt</b>  - one liner separated by commas containing: the id of contig, if it circularises or not (True or False), start coordinate of circularisation, end coordinate of circularisation
 

### Citations ####



When using MitoHifi, please cite this github page and

Please cite MitoFinder:

- Allio, R, Schomaker‐Bastos, A, Romiguier, J, Prosdocimi, F, Nabholz, B, Delsuc, F. MitoFinder: Efficient automated large‐scale extraction of mitogenomic data in target enrichment phylogenomics. Mol Ecol Resour. 2020; 00: 1– 14. https://doi.org/10.1111/1755-0998.13160

And for tRNAs annotation:

- Laslett, D., & Canbäck, B. (2008). ARWEN: a program to detect tRNA genes in metazoan mitochondrial nucleotide sequences. Bioinformatics, 24(2), 172-175.

-------------------
 
 For more information mu2@sanger.ac.uk
