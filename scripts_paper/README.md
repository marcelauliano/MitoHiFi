# 1. Reassembling mitogenomes with MitoHiFi
The scripts in this folder allow for the re-run of all mitogenomes assembled with MitoHiFi and presented in the manuscript _MitoHiFi: a python pipeline for mitochondrial genome assembly from PacBio High Fidelity reads_.

 | :warning: WARNING          |
|:---------------------------|
| **Caution: Potential for downloading several terabytes of data.**      |
 


## 1.1. Reassembling mitogenomes with MitoHiFi

All the species assembled with MitoHiFi with default parameters can be reassembled as follows. One will need the file [AdditionalTable1-dtol.tsv](AdditionalTable1-dtol.tsv), [singularity](https://docs.sylabs.io/guides/2.6/user-guide/singularity_and_docker.html) and [sra-tools](https://github.com/ncbi/sra-tools) to be installed and on the PATH. 

The following bash script will create a folder for each species, and inside each it will download the PacBio data from [SRA](https://www.ncbi.nlm.nih.gov/sra) and run MitoHiFi.


[./runMitoHiFi_default.sh](runMitoHiFi_default.sh)


:white_check_mark: Please ensure that you modify the BASE_DIR variable within the script and the INPUT_FILE if you intend to run it for a specific number of species. The script is configured to run with 4 threads by default. To adjust the number of threads, modify the -t parameter on line 59.

## 1.2. Reassemnbling Funghi mitogenomes with MitoHiFi

For the species _Mucor piriformis, Flammulina velutipes, Pleurotus ostreatus_ and _Agaricus bisporus_, reads were first assembled with MBG (parameters  -k 1001 -w 250 -a 5 -u 150) to obtain contigs that were then input to MitoHiFi with the -c flag. Bellow you find a bash script to reproduce these runs.

For this run you have to have [singularity](https://docs.sylabs.io/guides/2.6/user-guide/singularity_and_docker.html), [MBG](https://github.com/maickrau/MBG) and [sra-tools](https://github.com/ncbi/sra-tools) installed and on your PATH.

[./run_MitoHiFi_funghi.sh](run_MitoHiFi_funghi.sh)

!Note that the above is running with file [AdditonalFile1-funghi.tsv](AdditonalFile1-funghi.tsv). One needs to modify paramter -t in lines 41 and 69 of [run_MitoHiFi_funghi.sh](run_MitoHiFi_funghi.sh) to run it with more threads.






