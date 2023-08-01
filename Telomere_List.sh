#!/bin/bash
##---------------------------------------------------------------------
## This script is for separating telomere sequences out on MinION data (.fast5) using
##  the GPU nodes on EASLEY.
##  Brynleigh Payne. July 2023. This was modified to call the Brown Anole chromosome telomere Data
##    Thanks to Auburn Uni HPC administration for feedback to get this to work.
##    Some of the config has come from information here:
##-------------------- Parameters for running on EASLEY ----------------
#SBATCH --job-name=Telomere_List           # job name
#SBATCH --nodes=2                         # node(s) requried for job
#SBATCH --ntasks=10                       # number of tasks across all nodes
#SBATCH --partition=general                  # name of partition to submit job
#SBATCH --time=10:00:00                    # Run time (D-HH:MM:SS)
#SBATCH --output=job-%jTelomere_List_Test8.out          # Output file. %j is replaced with job ID
#SBATCH --error=job-%jTelomere_List_Test8.err           # Error file. %j is replaced with job ID
#SBATCH --mail-type=ALL                   # will send email for begin,end,fail
#SBATCH --mail-user=bep0022@auburn.edu
##---------------------------------------------------------------------

##Linearize all of the input files
##Separate out a specific sequence motif (TTAGGG)
#awk '/^>/ {printf("%s%s\t",(N>0?"\n":""),$0);N++;next;} {printf("%s",$0);} END {printf("\n");}' *.fastq |\
#awk -F '\t' '{if(index($4,”TTAGGG”)!=0) printf("%s\n%s\n",$1,$4);}’ > Telomere_TTAGGG.fa

#For one barcode
grep -B1 TTAGGGTTAGGGTTAGGG BC22.fastq > BC22_TTAGGG.fastq
#For all of the barcodes
grep -B1 TTAGGGTTAGGGTTAGGG BC**.fastq > TTAGGG.fastq

##Separate out a specific sequence motif (CCCTAA)
#For one of the barcodes
grep -B1 CCCTAACCCTAACCCTAA BC22.fastq > BC22_CCCTAA.fastq
#For all of the barcodes
grep -B1 CCCTAACCCTAACCCTAA BC**.fastq > CCCTAA.fastq
