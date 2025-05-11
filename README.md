# MitoHiFi 

------ This is v3.2.2 ------

MitoHiFi v3.2.2 is a python pipeline distributed under [MIT License](LICENSE) !


MitoHiFi was first developed to assemble the mitogenomes for a wide range of species in the Darwin Tree of Life Project (DToL) 

![](./docs/dtol-logo-round-300x132.png)


Find out more [Darwin Tree of Life data portal](https://portal.darwintreeoflife.org/)

[![MitoHiFi Integration Test](https://github.com/marcelauliano/MitoHiFi/actions/workflows/github-actions-integration-test.yml/badge.svg)](https://github.com/marcelauliano/MitoHiFi/actions/workflows/github-actions-integration-test.yml)

--------------------------------------

## 1. Background
**MitoHiFi is a python workflow that assembles mitogenomes from Pacbio HiFi reads.**

With MitoHiFi v3.2 you can start from the raw Pacbio HiFi reads (flag **-r**) or from the assembled contigs (flag **-c**). You also need a reference mitochondria sequence in FASTA and [GenBank format](https://www.ncbi.nlm.nih.gov/Sitemap/samplerecord.html). We provide an internal script (findMitoReference.py) that can be used to find and download the most closely-related reference genome for your species from NCBI..



*The availability of high-quality long reads, such as PacBio HiFi, greatly simplifies the assembly of accurate mitogenomes. However, due to the circular nature of the molecule, mitogenomes are typically assembled redundantly, resulting in multicopy mitogenome contigs. To address this, we have developed a pipeline specifically designed to finalize the assembly and the annotation of the mitogenome. Our pipeline takes into account heteroplasmy, aiming to assemble and annotate all mtDNA variants present in your sample. Among these variants, MitoHiFi selects a representative as the final genome assembly based on several criteria including circularization and gene completeness. Additionally, MitoHiFi provides multiple intermediate outputs such as coverage and annotation plots and a multiple sequence alignment (MSA) of all the variants, which facilitates the analysis of mitochondrial heteroplasmy.*



MitoHiFi v3.2 will:

(i) Extract mitochondrial reads and assemble them with hifiasm (using the -r flag) or identify mitochondrial contigs among assembled contigs (using the -c flag).

(ii) Identify and separate NUMTS (Nuclear Mitochondrial DNA Sequences) from genuine mitogenome contigs.

(iii) Generate a circularized, non-redundant, and annotated version of all mitochondrial variants present in your sample.

(iv) Select a representative as the final mitochondrial genome assembly.

(v) Plot an annotation image and a coverage plot of the reads (if initiated with the -r flag).

You can view a diagram illustrating the general MitoHiFi workflow [here](./docs/Figure1.png). For documentation on the behaviour of each script visit [here](./docs/scripts_documentation.pdf). 

-----

## 2. Installation

Below, we describe the three different ways to install MitoHiFi.

### 2.1 Using Docker and Singularity

We provide a Docker container for MitoHiFi. The container is built using GitHub Actions and can be obtained by running the following command:

```
docker pull ghcr.io/marcelauliano/mitohifi:master
```

Once the container is pulled, you can execute MitoHiFi within [Docker](https://docs.sylabs.io/guides/2.6/user-guide/singularity_and_docker.html). If Docker is not available on your cluster computer, you can use Singularity instead. Different environments may require slightly different Singularity parameters, but a generic Singularity command would be as follows:

```
singularity exec --bind /path/to/container_directory:/path/to/container_directory docker://ghcr.io/marcelauliano/mitohifi:master mitohifi.py -h
```

### 2.2 Conda: Partially Installing Dependencies with Conda

We are unable to create a complete Conda recipe for MitoHiFi due to conflicting Python versions required by different dependencies, making it impossible to coexist within the same Conda environment. However, we provide a partial Conda recipe that installs most of the dependencies. To complete the installation, follow these steps:

1. Install MitoFinder and/or MITOS outside of Conda.
2. Ensure MitoFinder and/or MITOS are added to the PATH before starting the run.
Please note that MitoFinder and/or MITOS should be installed separately and made accessible via the PATH environment variable to ensure their proper integration with MitoHiFi. Once those are installed, do:
```
#Clone MitoHiFi git repo
git clone https://github.com/marcelauliano/MitoHiFi.git

#create a conda environment with our yml file that is inside MitoHiFi/environment
conda env create -n mitohifi_env -f MitoHiFi/environment/mitohifi_env.yml 
```

Add MitoFinder and/or MITOS to the PATH and then activate your mitohifi_env conda environment.

To run MitoHiFi, do:

```
(mitohifi_env) python MitoHiFi/src/mitohifi.py -h
```

### 2.3 Manually install all dependencies

This is the least recommended way to install MitoHiFi, but below you will find a list of software that needs to be installed and added to your PATH before installing MitoHiFi. The software versions provided are the latest ones we have tested and confirmed to be working:

  - python=3.7
  - samtools=1.11
  - cd-hit=4.8.1
  - minimap2=2.19
  - hifiasm=0.19.5
  - mafft=7.520
  - biopython=1.79
  - matplotlib=3.2.2
  - dna_features_viewer=3.1.2
  - pandas=1.3.5
  - bedtools=2.31.0
  - pillow=6.2.1
  - bcbio-gff=0.7.0
  - MitoFinder=v1.4.0
  - MITOS=2.1.0
  - ncbi-blast+

Once all dependencies are installed and in your PATH, git clone MitoHiFi and execute it:

```
git clone https://github.com/marcelauliano/MitoHiFi.git
python MitoHiFi/src/mitohifi.py -h
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

## 4. Running MitoHiFi with the test data

To run MitoHiFi, you first need a closely related mitochondrial sequence in both FASTA and GenBank formats. We provide a script that can assist you in retrieving this data from the NCBI database. When provided with the species name you are assembling, the script will search for the nearest available mitochondrial assembly on NCBI. By default, the script looks for an exact match of the species for an available mitochondrial assembly. If an exact match is not found, it continues the search for a phylogenetically close candidate based on NCBI taxonomy.

### 4.1 Running MitoHiFi with test dataset starting from reads (-r)

Using the file ilDeiPorc1.reads.100.fa in our test data as an input example, you would do:

```
findMitoReference.py --species "Deilephila porcellus" --outfolder /path/to/outputdir --min_length 14000
```
This command will download OQ694980.1.fasta and OQ694980.1.gb that you should use for flags **-f** and **-g** in the main pipeline. Attention: once the NCBI is updated with more mitogenomes, findMitoReference.py might download something else. 

Now run MitoHiFi with 4 CPUs (change -t to change CPU numbers):

```
python mitohifi.py -r MitoHiFi/tests/ilDeiPorc1.reads.100.fa -f OQ694980.1.fasta -g OQ694980.1.gb -t 4 -o 5 
```

### 4.1 Running MitoHiFi with test dataset starting from contigs (-c)

To test starting from contigs, you are going to use the file MitoHiFi/test/ilPhaBuce1_contig.fa . For that, first you need to download a reference close to that species. For that, do:

```
findMitoReference.py --species "Phalera bucephala" --outfolder /path/to/outputdir --min_length 14000
```

Once the close reference is downloaded, run MitoHiFi:

```
python mitohifi.py -c MitoHiFi/tests/ilPhaBuce1_contig.fa -f NC_072273.1.fasta -g NC_072273.1.fasta -t 4 -o 5 
```

 
## 5. Output files
### 5.1 Main Outputs  

MitoHifi will produce a series of folders with the results. The main results will be in your working folder and they are:
- **final_mitogenome.fasta** - the final mitochondria circularized and rotated to start at tRNA-Phe
- **final_mitogenome.gb** - the final mitochondria annotated in GenBank format.  
- **final_mitogenome.coverage.png** - the sequencing coverage throughout the final mitogenome  
- **final_mitogenome.annotation.png** - the predicted genes throughout the final mitogenome
- **contigs_annotations.png** - annotation plots for all potential contigs
- **coverage_plot.png** - reads coverage plot of filtered reads mapped to all potential contigs
- **contigs_stats.tsv** - containing the statistics of your assembled mitos such as the number of genes, size, whether it was circularized or not, if the sequence has frameshifts and etc...
- **shared_genes.tsv** - show comparison of annotation between close-related mitogenome and all potential contigs assembled

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

When using MitoHiFi, please cite our paper: 

Uliano-Silva, M., Ferreira, J.G.R.N., Krasheninnikova, K. et al. MitoHiFi: a python pipeline for mitochondrial genome assembly from PacBio high fidelity reads. BMC Bioinformatics 24, 288 (2023). https://doi.org/10.1186/s12859-023-05385-y 
 
and

Please cite MitoFinder if you use the default annotation tool:

- Allio, R, Schomaker‐Bastos, A, Romiguier, J, Prosdocimi, F, Nabholz, B, Delsuc, F. MitoFinder: Efficient automated large‐scale extraction of mitogenomic data in target enrichment phylogenomics. Mol Ecol Resour. 2020; 00: 1– 14. https://doi.org/10.1111/1755-0998.13160

or MITOS:
- M. Bernt, A. Donath, F. Jühling, F. Externbrink, C. Florentz, G. Fritzsch, J. Pütz, M. Middendorf, P. F. Stadler MITOS: Improved de novo Metazoan Mitochondrial Genome Annotation Molecular Phylogenetics and Evolution 2013, 69(2):313-319. 

And for tRNAs annotation:

- Laslett, D., & Canbäck, B. (2008). ARWEN: a program to detect tRNA genes in metazoan mitochondrial nucleotide sequences. Bioinformatics, 24(2), 172-175.

-------------------
 
 For more information on python code and pipeline: mu2@sanger.ac.uk and jf18@sanger.ac.uk
 
 Questions about the Docker container: Marcela Uliano-Silva mu2@sanger.ac.uk



## 9. Watch a lecture on MitoHiFi

Want to know everything about how to run MitoHiFi for animals, funghi and plants? Why the pipeline stops sometimes? Caveats and best practices? Have a watch here: https://youtube.com/watch?v=1NWHC2zkRmg 

