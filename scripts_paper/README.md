# 1. Reassembling mitogenomes with MitoHiFi
The scripts in this folder allow for the re-run of all mitogenomes assembled with MitoHiFi present in the manuscript _MitoHiFi: a python pipeline for mitochondrial genome assembly from PacBio High Fidelity reads_.

Caution is adivised as it would mean dowloading several Terabites of data in some cases.

## 1.1. Reassembling Darwin Tree of Life mitogenomes with MitoHiFi

As stated in the paper, all mitogenomes were assembled with default MitoHiFi parameters, apart from the funghi. Here we present a bash script that would allow one to re-assemble all mitogenomes from the Darwin Tree of Life Program. Next you find another bash script to assemble (i) the funghi species. 

One will need the file AdditionalTable1-dtol.tsv and singularity to be installed int the PATH to run the following:

```
[./runMitoHiFi_default.sh](runMitoHiFi_default.sh)
```

! Note that one needs to change the BASE_DIR variable inside the script and INPUT_FILE in case one is running it only for a number of species.

## 1.2. Reassemnbling Funghi mitogenomes with MitoHiFi


 
 
 Finalized mitogenomes IDs can be found in column two of AdditionalTable1.tsv.

But if one want's to re-run MitoHiFi for all or for a number of species, one can use the AdditionalTable1.tsv - or a portion of it - and use the base script presented as follows. To run the following one need to have singularity installed on the PATH.

In the folder where you want to 
