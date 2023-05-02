# MitoHiFi 

------ This is v3.0.0 -------

MitoHiFi is a python pipeline distributed under the [MIT License](LICENSE)


MitoHiFi was first developed to assemble the mitogenomes for a wide range of species in the Darwin Tree of Life Project (DToL)  ![](dtol-logo-round-300x132.png)


Find out more [./docs/Darwin Tree of Life data portal](https://portal.darwintreeoflife.org/)

--------------------------------------

## 1. Background
**MitoHiFi is a python workflow that assembles mitogenomes from Pacbio HiFi reads.**

With MitoHiFi v3.0.0 you can start from the raw Pacbio HiFi reads (flag **-r**) or from the assembled contigs (flag **-c**). You also need a reference mitochondria sequence in FASTA and [GenBank format](https://www.ncbi.nlm.nih.gov/Sitemap/samplerecord.html). We provide an internal script (findMitoReference.py) that can be used to find and download the most closely-related reference genome for your species from NCBI.



*The dissemination of high-quality long reads - such as PacBio HiFi - makes the assembly of high-quality mitogenome straightforward. Because of the circular nature of the molecule, however, the mitogenome is usually assembled redundantly resulting in multiple-copy mitogenome contigs. This pipeline was developed to finalise the assembly and annotate the mitogenome. It also considers heteroplasmy and aims to assemble and annotate all the mtDNA variants in your sample, among which MitoHiFi v2.2 will choose a representative as the final genome assembly according to the similarity to the reference genome and gene completeness. In addtion, an multiple sequence alignment (MSA) of all the variants is provided to facilitate the analysis of mitochondria heteroplasmy.*



MitoHiFi v3.0.0 will:

(i) extract mito reads and assemble them with hifiasm (flag **-r**), or find the mito contigs among assembled contigs (flag **-c**)<br />    

(ii) indentify and separate NUMTS from real mitogenome contigs<br />  

(iii) generate a circularized, non-redudant and annotated version of all the mitochondria variants present in your sample<br />

(iv) choose a representative as the final mitochondria genome assembly<br />

(iv) plot annotation image and a reads coverage plot (if started with **-r**)

-----

## 2. Installation

There are two ways to install MitoHiFi v3.0.0:

(i) manually - and then you will need to have all the dependencies installed in your PATH.  

(ii) with a docker image so you can run it with a singularity (highly recommended).

### 2.1 Manual installation - Dependencies

All the software listed bellow have to be installed and exported to your PATH. 
MitoHifi v2.3 was developed and tested with the following software versions:

- BLAST version 2.6.0
- Bedtools version 2.30.0 
- [MitoFinder version 1.4](https://github.com/RemiAllio/MitoFinder) (default annotation tool) or [MITOS](https://gitlab.com/Bernt/MITOS) version 2.1.0 (alternative annotation tool)  
- Biopython version 1.78
- [DNA Features Viewer](https://github.com/Edinburgh-Genome-Foundry/DnaFeaturesViewer) version 3.1.0
- Matplotlib version 3.5.2 
- Pandas version 1.1.3
- [Pillow](https://pillow.readthedocs.io/en/stable/) version 6.2.1  
- MAFFT version 7.475 
- HiFiasm 0.16.1-r375 
- CD-HIT version 4.8.1 
- samtools version 1.7 
- minimap version 2.17-r941
- [MITOS](https://pypi.org/project/mitos/) 


All above dependencies (except MitoFinder) can be installed using [conda](https://conda.io/projects/conda/en/latest/index.html), and we share a YML file to facilitate their installation. The user can create a conda environment named `mitohifi_env` with all dependencies installed by running:  

```
conda env create -n mitohifi_env -f mitohifi_env.yml 
```

Due to incompatibility of python versions, MitoFinder needs to be manually installed following its official instructions. After that, the user must add MitoFinder to PATH:  

```
export PATH=</path/to/MitoFinder>:${PATH}
```

where `</path/to/MitoFinder>` needs to be replaced with the path where MitoFinder was installed. 

**Once all dependencies are installed, install MitoHiFi v2.3 (Linux)**

```
git clone https://github.com/marcelauliano/MitoHiFi.git
```

### 2.2 Running MitoHiFi v3.0.0 from a Singularity with a Docker image


We have wrapped up MitoHiFi v3.0.0 code into a singularity container. We recommend using singularity versions => 3.7, as lower versions do not support spaces in the arguments, and you would not be able to pass more than one set of reads to the flag **-r**

MitoHiFi.v3.0.0 siungularity image should be run as:

```
 singularity exec --bind /lustre/:/lustre/ docker://ghcr.io/marcelauliano/mitohifi:master mitohifi.py -r "/data/f1.fasta /data/f2.fasta /data/f3.fasta" -f /data/reference.fasta -g /data/reference.gb -t 10 -o 2 
```

The script for quering reference .fasta and .gb files from NCBI is incorporated into the singularity image and can be called as follows:

```
 singularity exec --bind /lustre/:/lustre/ docker://ghcr.io/marcelauliano/mitohifi:master findMitoReference.py --species "Cryptosula pallasiana" --email your@email.for.ncbi.db.query --outfolder /data/ --min_length 16000 
```

## 3. Parameter list

```
usage: MitoHiFi [-h] (-r <reads>.fasta | -c <contigs>.fasta) -f
                <relatedMito>.fasta -g <relatedMito>.gbk -t <THREADS> [-d]
                [-a {animal,plant,fungi}] [-p <PERC>] [-m <BLOOM FILTER>]
                [--max-read-len MAX_READ_LEN] [--mitos]
                [--circular-size CIRCULAR_SIZE]
                [--circular-offset CIRCULAR_OFFSET] [-winSize WINSIZE]
                [-covMap COVMAP] [-v] [-o <GENETIC CODE>]

required arguments:
  -r <reads>.fasta      -r: Pacbio Hifi Reads from your species
  -c <contigs>.fasta    -c: Assembled fasta contigs/scaffolds to be searched
                        to find mitogenome
  -f <relatedMito>.fasta
                        -f: Close-related Mitogenome is fasta format
  -g <relatedMito>.gbk  -k: Close-related species Mitogenome in genebank
                        format
  -t <THREADS>          -t: Number of threads for (i) hifiasm and (ii) the
                        blast search

optional arguments:
  -d                    -d: debug mode to output additional info on log
  -a {animal,plant,fungi}
                        -a: Choose between animal (default) or plant
  -p <PERC>             -p: Percentage of query in the blast match with close-
                        related mito
  -m <BLOOM FILTER>     -m: Number of bits for HiFiasm bloom filter [it maps
                        to -f in HiFiasm] (default = 0)
  --max-read-len MAX_READ_LEN
                        Maximum lenght of read relative to related mito
                        (default = 1.0x related mito length)
  --mitos               Use MITOS2 for annotation (opposed to default
                        MitoFinder
  --circular-size CIRCULAR_SIZE
                        Size to consider when checking for circularization
  --circular-offset CIRCULAR_OFFSET
                        Offset from start and finish to consider when looking
                        for circularization
  -winSize WINSIZE      Size of windows to calculate coverage over the
                        final_mitogenom
  -covMap COVMAP        Minimum mapping quality to filter reads when building
                        final coverage plot
  -v, --version         show program's version number and exit
  -o <GENETIC CODE>     -o: Organism genetic code following NCBI table (for
                        mitogenome annotation): 1. The Standard Code 2. The
                        Vertebrate MitochondrialCode 3. The Yeast
                        Mitochondrial Code 4. The Mold,Protozoan, and
                        Coelenterate Mitochondrial Code and the
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
                        Code
```

## 4. Running MitoHiFi

1-) To run MitoHiFi, first you need a close-related mitochondria in fasta and genbank format. We have a script that can help you fetch this data from NCBI database. Given the name of the species you are assembling, the script is going to look for the closest mitochondria it can find on NCBI. By default the script searches for an available mitochondria assembly of exactly same species. If it's not available the search goes on for a phylogenetically close candidate (based on NCBI taxonomy).

Using the species in our test data as an example, you would do:

```
findMitoReference.py --species "Phalera bucephala" --email your@email.for.ncbi.db.query --outfolder /data/ --min_length 16000
```
This command will give you NC_016067.1.fasta and NC_016067.1.gb that you can use for flags **-f** and **-g** in the main pipeline.

2-) Now, you need to decide if you want to run MitoHiFi v2.3 from:
(i) raw reads, in which case the pipeline will map the reads (flag **-r**) to the reference genome of the closely-related species to pull out mito-reads and exclude possible NUMTS and then assemble them using Hifiasm, or 
(ii) assembled contigs from the whole genome Pacbio HiFi sequence data, in which case the pipeline will take the contig file in FASTA format as input (flag **-c**) and map the contigs to the reference genome of the closely-related species to pull-out potential mitogenome contigs.

2.1-) If you are starting from raw reads, **-r** flag is required to provide the raw PacBio HiFi reads. Here is an example:

```
'singularity exec --bind /software/:/software/ docker://ghcr.io/marcelauliano/MitoHiFi:main mitohifi.py -r "f1.fasta f2.fasta f3.fasta" -f reference.fasta -g reference.gb -t <int> -o <int>'
 ```

2.2-) If you are starting from assembled contigs, **-c** flag is required to provide the FASTA file for assembled contigs. Here is an example:

```
'singularity exec --bind /software/:/software/ docker://ghcr.io/marcelauliano/MitoHiFi:main mitohifi.py -c contigs.fasta -f reference.fasta -g reference.gb -t <int> -o <int>'
```

### 4.1 Running MitoHiFi with test data

- Use your singularity container image or have all the dependencies in your path then,
- In the exampleFiles folder, the fasta and .gb file for NC_016067.1 will be your **-f** and **-g** inputs, respectively. Remember you could have gotten those files with the script findMitoReference.py.
- Now run the test with the assembled contigs called test.fa

```
'python mitohifi.py -c exampleFiles/test.fa -f exampleFiles/NC_016067.1.fasta -g exampleFiles/NC_016067.1.gb -t 1 -o 5'
 ```
 
## 5. Output files
### 5.1 Main Outputs  

MitoHifi will produce a series of folders with the results. The main results will be in your working folder and they are:
- final_mitogenome.fasta - the final mitochondria circularized and rotated to start at tRNA-Phe
- final_mitogenome.gb - the final mitochondria annotated in GenBank format.  
- final_mitogenome.coverage.png - the sequencing coverage throughout the final mitogenome  
- final_mitogenome.annotation.png - the predicted genes throughout the final mitogenome
- contigs_annotations.png - annotation plots for all potential contigs
- coverage_plot.png - reads coverage plot of filtered reads mapped to all potential contigs
- contigs_stats.tsv - it will show you the statistics of your assembled mitos such as the number of genes, size, whether it was circularized or not, if the sequence has frameshifts and etc...
- shared_genes.tsv - show comparison of annotation between close-related mitogenome and all potential contigs assembled

### 5.2 Further outputs

* Folder **contigs_filtering** will contain three outputs:
 
  - parsed_blast.txt 
  - parsed_blast_all.txt
  - contigs_ids.txt
  - contigs.blastn - outfmt 6 blast output plus 2 extra columns containing respectively length_of_query and length_of_subject 
 
Columns descriptions of <b>parsed_blast.txt</b> and <b>parsed_blast_all.txt</b>: tab separated files with 4 columns as follows

    - qseqid - the ID of your input contigs
    - %q_in_match - a percentage of the length of your contig in a blast match with the close related species reference
    - leng_query - length of your contig
    - s_length - lenght of the close-related mitogenome given
    - perc - how much percente the leng_query is in relation to s_length
 
* Folder **contigs_circularization** will contain files related to the circularization of the mito contigs. The most important file is 

  - all_contigs.circularisationCheck.txt - describes the points of circularization for all contigs. The script is run iteratively until no more circularizatin is found, which will happen when you meet (*False, -1, -1*). The columns are: the contig id, if it circularises or not (*True* or *False*), start coordinate of circularisation, end coordinate of circularisation.

* Folder **potential_contigs** will contain a folder for each contig present in parsed_blast.txt. Within each contig folder you will find the circularized and annotated mitosequence for that contig.

* Folder **final_mitogenome_choice** will contain a few files, the most important one is

  - all_mitogenomes.rotated.aligned.fa - this is an aligment of all the mithocondrial sequences assembled by the pipeline. It is possible you will find heteroplasmy in your sample, in which case you will have more than one version of the final mito presented. The pipeline chooses a final representative by a majority rule, using cdhit-est to cluster sequenvces at a 80% identitty and chosing the largest one in that cluster as the final. If you want to study heteroplasmy of your sample, please investigate the *all_mitogenomes.rotated.aligned.fa* file further, and all the results in the **potential_contigs** folder.


* Folder **coverage_mapping** will contain bam files for a quick inspection on IGV-like softwares, the most important are:  
  - HiFi-vs-final_mitogenome.sorted.bam - contains mapping information of filtered HiFi reads against the final mitogenome 
  - HiFi-vs-potential_contigs.sorted.bam - contains mapping information of HiFi reads against all potential contigs

* If you run the pipeline with flag **-r** you will have a further folder called **reads_mapping_and_assembly** which will contain

  - gbk.HiFiMapped.bam.fasta - all reads that mapped to the closely-related mito
  - gbk.HiFiMapped.bam.filtered.fasta - mapped reads filtered by size. We remove any reads that are larger than the size of the close-related mito (or as specified by **--max-read-len** ) as a rough way to filter out numpts
  - hifiasm.contigs.fasta - final hifiasm primary and alternate contigs concatenated. This is the file used to find your mitos.

## 6. New parameter for plants!!

MitoHiFi is still not optmized to assemble a plant mitochondria or chloroplast. But if you have a contig you consider to be one of those, you can use MitoHiFi with the flag **-c** to finalize your organelle! It will circularize (if that is the case) and annotate, and output statistics for you. To do this, you need to use the parameter **-a plant** when calling the main script mitohifi.py .

Also, the script findMitoReference.py can now search for a chloroplast instead of a mitochondria. Use flag **-t chloroplast** 

## 7. Important parameter to change and test 

### 7.1 (-p)

Mitohifi is going to pull possible mito contigs by blasting your contigs with the closely-related mito reference genome. The Default parameter **-p** is going to chose any contig which has 50% or more of its length in the blast match. This is the default because with invertebrate taxa from the Darwin Tree of Life project we have been seeing that the repetitive portion of the mitogenomes is not very conserved between some taxa. In these cases, a more stringent **-p** ends up excluding real mito sequences. Nevertheless, if you are working with more conserved taxa - such as mammals and other vertebrates - use higher **-p** (such as 80 or 90) for better results.

### 7.2. Annotation run with MITOS

The default annotator for MitoHiFi is MitoFinder, but the user can annotate with MITOS by flagging **--mitos** while starting a MitoHiFi run.

### 7.3 Plots

The user can change **-winSize** and **-covMap** parameters to tun the final coverage plots. 

## 8. Citations

When using MitoHiFi, please cite our preprint: 

MitoHiFi: a python pipeline for mitochondrial genome assembly from PacBio High Fidelity reads
Marcela Uliano-Silva, João Gabriel R. N. Ferreira, Ksenia Krasheninnikova, Darwin Tree of Life Consortium, Giulio Formenti, Linelle Abueg, James Torrance, Eugene W. Myers, Richard Durbin, Mark Blaxter, Shane A. McCarthy
bioRxiv 2022.12.23.521667; doi: https://doi.org/10.1101/2022.12.23.521667
 
and

Please cite MitoFinder:

- Allio, R, Schomaker‐Bastos, A, Romiguier, J, Prosdocimi, F, Nabholz, B, Delsuc, F. MitoFinder: Efficient automated large‐scale extraction of mitogenomic data in target enrichment phylogenomics. Mol Ecol Resour. 2020; 00: 1– 14. https://doi.org/10.1111/1755-0998.13160

or MITOS:
- M. Bernt, A. Donath, F. Jühling, F. Externbrink, C. Florentz, G. Fritzsch, J. Pütz, M. Middendorf, P. F. Stadler MITOS: Improved de novo Metazoan Mitochondrial Genome Annotation Molecular Phylogenetics and Evolution 2013, 69(2):313-319. 

And for tRNAs annotation:

- Laslett, D., & Canbäck, B. (2008). ARWEN: a program to detect tRNA genes in metazoan mitochondrial nucleotide sequences. Bioinformatics, 24(2), 172-175.

-------------------
 
 For more information on python code and pipeline: mu2@sanger.ac.uk and jf18@sanger.ac.uk
 
 Questions on the Singularity: Marcela Uliano mu2@sanger.ac.uk
