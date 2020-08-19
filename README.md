# Mito HiFi 

------ This is a beta version! -------


This repository contains software and scripts to circularise, cut and annotate mitogenomes assembled redudantly by assembly softwares such as HiCanu or Hifiasm.

The dissemination of PacBio HiFi reads made it possible for mitogenomes to be assembled straight away. Because of the circular nature of the molecule, however, the mitocontig is usually assembled rudandantly resulting in a contiguos many copies of the mitogenome (improve this description!!). This pipeline was developted to finalize the assembly and annotation of the mitogenome by (i) finding it among total assembled contigs, (ii) circularizing and cutting it to represent only one copy of the circular molecule and (iii) annotation it and presenting it in fasta and genbank format.


<b>Dependencies:</b>

- Blast (makeblastdb and blastn) have to be installed and export on your PATH
- mitoFinder: has to install [mitoFinder](https://github.com/RemiAllio/MitoFinder) and export on your PATH 
- Biopython

<b>Installation</b>

### Get and install MitoHiFi (Linux)

```

git clone https://github.com/marcelauliano/MitoHiFi.git

```

### Test run

```
cd MitoHiFi
cp ExampleFiles
ln -s ../scripts
sh run_mitoHiFi.sh -c test.fa -f NC_016067.1.fasta -g NC_016067.1.gb -t 1

```
### Required arguments

```
Usage: 'sh run_mitoHiFi.sh -c contigs.fasta -f close-related_mitogenome.fasta -g close-related_mitogenome.gb -t threads'
Parameters:	
	-c: assemnbled fasta contigs/scaffolds to be searched to find mitogenome
	-f: Close-related species mitogenome in fasta format
	-g: Close-related species mitogenome in genbank format 
	-t: Number of threads for the blast search 
 ```
 
 ### Description of each intermediate output
 
 ``` 
 will be written soon

 ```
 
 For more information mu2@sanger.ac.uk
