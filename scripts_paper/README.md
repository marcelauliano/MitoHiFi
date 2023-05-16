# 1. Reassembling mitogenomes with MitoHiFi
The scripts in this folder allow for the re-run of all mitogenomes assembled with MitoHiFi present in the manuscript _MitoHiFi: a python pipeline for mitochondrial genome assembly from PacBio High Fidelity reads_.

**Caution is advised as this could mean downloading several terabytes of data in certain cases.**

## 1.1. Reassembling Darwin Tree of Life mitogenomes with MitoHiFi

As stated in the paper, all mitogenomes were assembled with default MitoHiFi parameters, apart from the funghi. Here we present a bash script that would allow one to re-assemble all mitogenomes from the Darwin Tree of Life Program. Next you find another bash script to assemble (i) the funghi species. 

One will need the file AdditionalTable1-dtol.tsv and singularity to be installed int the PATH to run the following:


[./runMitoHiFi_default.sh](runMitoHiFi_default.sh)


! Note that one needs to change the BASE_DIR variable inside the script and INPUT_FILE in case one is running it only for a number of species. At the moment, the script is set to run with [AdditonalFile1-dtol.tsv](AdditonalFile1-dtol.tsv)

## 1.2. Reassemnbling Funghi mitogenomes with MitoHiFi

For the species _Mucor piriformis, Flammulina velutipes, Pleurotus ostreatus_ and _Agaricus bisporus_, reads were first assembled with MBG (parameters  -k 1001 -w 250 -a 5 -u 150) to obtain contigs that were then input to MitoHiFi with the -c flag. Bellow you find a bash script to reproduce these runs.

For this run you have to have singularity and [MBG](https://github.com/maickrau/MBG) installed and on your PATH.




