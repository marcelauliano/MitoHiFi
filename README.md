# MitoHiFi 

------ This is v2 -------

MitoHiFi is a python pipeline distributed under the [license](LICENSE)

--------------------------------------


## <b>MitoHiFi</b> is a python workflow that assembles a species mitogenome from Pacbio HiFi reads.

With Mitoifi.v2 you can start from raw Pacbio HiFi reads (flag -r) or from assembled contigs (flag -c). You will also going to need a close-related mitochondria sequence in fasta and gb format. We have an internal script that can download it for you from NCBI (findMitoReference.py).



*The dissemination of high-quality long reads - such as PacBio HiFi - makes the assembly of high-quality mitogenome straight forward. Because of the circular nature of the molecule, however, the mitocontig is usually assembled redundantly resulting in multiple-copy mitogenome-contigs. This pipeline was developed to finalise the assembly and annotation of the mitogenome. It will also dected different variants of the mitogenome present in your sample. At the end you are going to have all the variants assembled and annotated, and MitoHiFi.v2 is going to choose a final consensus sequence. In addtion, you will find an aligment of all the variants to facilitate your analysis of mitochondria heteroplasmy.*



MitoHifi v2 will:

(i) extract mito reads and assemble them with hifiasm (flag -r), or find the mito contigs among assembled contigs (flag -c)    
(ii) indentify and separate NUMTS from real mitochontigs  
(iii) generate a circularized, non-redudant and annotated version of all the mitochondria sequences present in your sample	
(iv) choose a final consensus as the final mitochondria


-----



### Installation

There are two ways to install MitoHifi.v2 at the moment; (i) mannually - and then you will need to have all the dependencies installed in your PATH, or with (ii) a singularity image.

### Manual installation - Dependencies

All the software listed bellow have to be installed and exported to your PATH. 
MitoHifi.v2 was developed and tested with the following software versions:

- Blast version 2.6.0 
- [MitoFinder version 1.4](https://github.com/RemiAllio/MitoFinder) 
- Biopython version 1.78 
- Pandas version 1.1.3 
- MAFFT version 7.475 
- HiFiasm 0.14-r312 
- CD-HIT version 4.8.1 
- samtools version 1.7 
- minimap version 2.17-r941  


### Once dependencies are done, install MitoHiFi.v2 (Linux)

```

git clone --branch mitohifi_v2 https://github.com/marcelauliano/MitoHiFi.git

```

### Running MitoHiFi.v2 from a Singularity image

We have wrapped up MitoHiFi.v2 code into a singularity container. We recommned using singularity versions => 3.7, as lower versions do not support spaces in the arguments, and you would not be able to pass more than one set of reads to the flag -r

MitoHiFi.v2 siungularity should be called as:

```
singularity exec --bind /path/on/disk/to/data/:/data/ /path/to/mitohifi-v2.sif  mitohifi_v2.py -r "/data/f1.fasta /data/f2.fasta /data/f3.fasta" -f /data/reference.fasta -g /data/reference.gb  -t 10 -o 2
```

Singluarity versions lower than 3.7 do not support spaces in the arguments, so if you want to pass several read datasets as in the example above use singularity version 3.7 or higher. 

The script for generating the reference files is incorporated into the singularity image and can be called as follows:

```
singularity exec --bind /path/on/disk/to/data/:/data/ /path/to/mitohifi-v3.sif  findMitoReference.py --species "Cryptosula pallasiana" --email your@email.for.ncbi.db.query --outfolder /data/ --min_length 16000 
```

### Required arguments

1-) To run this pipeline, first you need a close-related mitochondria in fasta and genbank format. We have a script that can help you find this input. Giving the name of the species you are assembling, the script is going to look for the closest mitochondria it can find on NCBI. You can give the parameter **-s** to the script if you would like to restrict your mitochondria search for species within your given genus, but this means the script can download partial mitochondrial sequences. Otherwise, without **-s**, the script is going to search for complete mitochondrias only and as close as possible to your species on interest.

Using the species in our test data as an example, you would do:

```
findMitoReference.py --species "Phalera bucephala" --email your@email.for.ncbi.db.query --outfolder /data/ --min_length 16000
```
This command will output you NC_016067.1.fasta and NC_016067.1.gb that you will use as flags **-f** and **-g** in the main pipeline.

2-) Now, you need to decide if you are running MitoHiFi.v2 from:
(i) raw reads, in which case the pipeline is going to map your reads to the close-related species (to pull out mito-reads and exclude possible NUMTS) and then assemble them using Hifiasm, or 
(ii) your already have a Pacbio HiFi assembly and you are going to give the contigs to the pipeline (flag **-c**).

2.1-) If you are starting from raw reads, your required arguments are:

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

2.2-) If you are starting from assembled contigs, your required arguments are:

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

### Running MitoHifi_v2 with test data

- Use your singularity container image or have all the dependencies in your path then,
- Download the data from the exampleFiles folder. The fasta and .gb file for NC_016067.1 will be your **-f** and **-g** inputs, respectively. Remember you could have gotten those files with the script findMitoReference.py.
- Now run the test with the assembled contigs called test.fa

```
'python mitohifi_v2.py -c test.fa -f NC_016067.1.fasta -g NC_016067.1.gb  -t <int> -o 5 '

 ```
 
 
 ### Outputs
 
MitoHifi will produce a series of folders with the results. The main result will be in your working folder and it constitutes of:
- final_mitogenome.fasta - the final mitochondria circularized and rotated to start at tRNA-Phe
- final_mitogenome.gb - the final mitochondria annotated in genbank format.

## Further outputs

 Folders:
 
 **contigs_filtering** will contain 3 outputs:
 
- parsed_blast.txt 
- parsed_blast_all.txt
- contigs_ids.txt
- contigs.blastn - outfmt 6 blast output plus 2 extra columns containing respectively length_of_query and length_of_subject 
 
 Columns descriptions of <b>parsed_blast.txt</b> and <b>parsed_blast.txt</b>: tab separated files with 4 columns as follows

 - qseqid - the ID of your input contigs
 - %q_in_match - a percentage of the length of your contig in a blast match with the close related species reference
 - leng_query - length of your contig
 - s_length - lenght of the close-related mitogenome given
 - perc - how much percente the leng_query is in relation to s_length
 
Folder **contigs_circularization** will contain files related to circularizing the mito contig. The most important file is 

- all_contigs.circularisationCheck.txt - describes the points of circularization for the chosen contigs. The script is ran iteractively until no more circularizatin is found, which will happen when you meet (False, -1, -1). The columns are: the contig id, if it circularises or not (True or False), start coordinate of circularisation, end coordinate of circularisation.

Folder **potential_contigs** will contain a folder for each contig present in parsed_blast.txt. Within each contig folder you will find the circularized and annotated mitosequence for that contig.


Folder **final_mitogenome_choice** will contain a few files, the most important one being

- all_mitogenomes.rotated.aligned.fa - this is an aligment of all the mithocondrial sequences assembled by the pipeline. Its possible you will find heteroplasmy in your sample, in which case you will have more than one version of the final mito represented. The pipeline chooses a final consensus by a majority rule, using cdhit-est to cluster sequenvces at a 80% identitty and chosing the largest one in that cluster as the final. If you want to study heteroplasmy of your sample, please investigate the *all_mitogenomes.rotated.aligned.fa* file further, and all the results in the **potential_contigs** folder.

## Important parameter to change and test (-p)

Mitohifi is going to pull possible mito contigs by blasting your contigs with the close-related mito. The Default parameter **-p** is going to chose any contig which has 50% or more of its length in the blast match. This is the default because with invertebrate taxa from the Darwin Tree of Life we have been seeing that the repetitive portion of the mitogenomes is not very conserved between some taxa. In these cases, a more stringent **-p** ends up excluding real mito sequences. Nevertheless, if you are working with more conserved taxa - such as mammals and other vertebrates - use higher -p (such as 80 or 90) for better results.


### Citations ####



When using MitoHifi, please cite this github page and

Please cite MitoFinder:

- Allio, R, Schomaker‐Bastos, A, Romiguier, J, Prosdocimi, F, Nabholz, B, Delsuc, F. MitoFinder: Efficient automated large‐scale extraction of mitogenomic data in target enrichment phylogenomics. Mol Ecol Resour. 2020; 00: 1– 14. https://doi.org/10.1111/1755-0998.13160

And for tRNAs annotation:

- Laslett, D., & Canbäck, B. (2008). ARWEN: a program to detect tRNA genes in metazoan mitochondrial nucleotide sequences. Bioinformatics, 24(2), 172-175.

-------------------
 
 For more information on python code and pipeline: mu2@sanger.ac.uk and jf18@sanger.ac.uk
 
 Questions on the Singularity: kk16@sanger.ac.uk
