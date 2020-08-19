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
sh MitoHiFi -c test.fa -f Nxxx.fasta -g Nxxx.gb -t 1

```
