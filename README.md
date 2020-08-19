# Mito HiFi 

This repository contains software and scripts to circularise, cut and annotate mitogenomes assembled redudantly by assembly softwares such as HiCanu or Hifiasm.

The dissemination of PacBio HiFi reads made it possible for mitogenomes to be assembled straight away. Because of the circular nature of the molecule, however, the mitocontig is usually assembled rudandantly resulting in a contiguos many copies of the mitogenome (improve this description!!). This pipeline was developted to finalize the assembly and annotation of the mitogenome by (i) finding it among total assembled contigs, (ii) circularizing and cutting it to represent only one copy of the circular molecule and (iii) annotation it and presenting it in fasta and genbank format.


<b>Dependencies:</b>

- Blast (makeblastdb and blastn) have to be installed and export them to your PATH
- mitoFinder: You have to install [mitoFinder](https://github.com/RemiAllio/MitoFinder) and export it to your PATH 

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
sh MitoHiFi -c test.fa -f NC_016067.1.fasta -g NC_016067.1.gb -t 1

```
### Required arguments

```
Usage: -c contigs.fasta -f close-related_mitogenome.fasta -g close-related_mitogenome.gb -t threads'
	-c: assemnbled fasta contigs/scaffolds to be searched to find mitogenome
	-f: Close-related species mitogenome in fasta format
	-g: Close-related species mitogenome in genbank format 
	-t: Number of threads for the blast search 
 ```
