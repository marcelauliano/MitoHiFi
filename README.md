# Mito HiFi 

------ This is a beta version! -------


This repository contains software and scripts to circularise, cut and annotate mitogenomes assembled redudantly by assembly softwares such as HiCanu or Hifiasm.

The dissemination of PacBio HiFi reads makes the assemble of high-quality mitogenome straigh forward. Because of the circular nature of the molecule, however, the mitocontig is usually assembled rudandantly resulting in a contiguos multiple-copies mitogenome (IMPROVE THIS DESCRIPTION!!). This pipeline was developted to finalize the assembly and annotation of the mitogenome by (i) finding it among total assembled nuclear contigs, (ii) circularising and cutting it to represent only one copy of the circular molecule and (iii) annotation and presenting it in fasta and genbank format.


<b>Dependencies:</b>

- Blast (makeblastdb and blastn) have to be installed and export on your PATH
- mitoFinder: has to install [mitoFinder](https://github.com/RemiAllio/MitoFinder) and export on your PATH 
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
 * mitogenome.fasta *  - this is your final mitogenome in fasta format
 
Inside the folder mitogenome.annotation/mitogenome.annotation_Final_Results/ you find mitogenome.gb, which is the annotation of your mitogenome performed by mitofinder and outputed in genbank format. 
 
```
 
 ### Description of each intermediate output
 
 ``` 
 circularisationCheck.txt ...
 parsed_blast.txt ...
 !!!!! will be written soon !!!!!

 ```
 
 For more information mu2@sanger.ac.uk
