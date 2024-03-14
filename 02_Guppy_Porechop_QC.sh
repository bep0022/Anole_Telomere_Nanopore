#!/bin/bash
##---------------------------------------------------------------------
## This script is for doing High Accuracy Basecalling on MinION data (.fast5) using
##  the GPU nodes on EASLEY. Run either on gpu4 or gpu2 on Easley. It will not run on nova_gpu.
##  Tonia Schwartz. July 2023. This was written to call the Seedbeetle Mitochondrial Genome Data
##    Thanks to Gabrielle de la Silva, Dasia Simpson, and Auburn Uni HPC administration for feedback to get this to work.
##    Some of the config has come from information here: https://hackmd.io/@Miles/S12SKP115.
##-------------------- Parameters for running on EASLEY ----------------
#SBATCH --job-name=Guppy_Pore_QC_Test4           # job name
#SBATCH --nodes=2                         # node(s) requried for job
#SBATCH --ntasks=10                       # number of tasks across all nodes
#SBATCH --partition=general                  # name of partition to submit job
#SBATCH --time=10:00:00                    # Run time (D-HH:MM:SS)
#SBATCH --output=job-%j_G_P_QCTest4.out          # Output file. %j is replaced with job ID
#SBATCH --error=job-%j_G_P_QCTest4.err           # Error file. %j is replaced with job ID
#SBATCH --mail-type=ALL                   # will send email for begin,end,fail
#SBATCH --mail-user=bep0022@auburn.edu

##-------------------------------------------------------------------------------------
  #### /tools/guppy-6.4.6gpu/bin/guppy_basecaller
  #### source /opt/asn/etc/asn-bash-profiles-special/modules.sh

module load guppy/6.4.6gpu
module load fastqc
module load python/anaconda/3.8.6   #for NanoPlot, porechop, medaka, quast

##----- Define variable for directories.
  ## Source = This is where the directory of fast5 data files.
  ## DD = this is where the High Accuracy Basepair Called Data will go.
  ## WD = the working directory in scratch

SOURCE=/hosted/biosc/SchwartzLab/Nanopore/AnoleEmbryoOldiesTelomere/20230605_2000_MN39112_FAS96157_fb416b6a
DD=/hosted/biosc/SchwartzLab/Nanopore/AnoleEmbryoOldiesTelomere/BaseCalled_Guppy10_GPU_HAC_230622/pass
WD=/scratch/bep0022

##------ Make the directory for my results in my home folder
   ## the -p will  Create intermediate directories as required.

mkdir -p ${DD}
mkdir -p ${WD}/NanoPlot
mkdir -p ${WD}/FastQC           # this is for the bulk data after HAC, before Porechop
mkdir -p ${WD}/FastQC_postPorechop	# this is for the data after Porechop
mkdir -p ${WD}/Porechopped

##------ This is where Dr. Schwartz used guppy basecaller.

#####---------- Nanoplot to summarize the basecalled data
NanoPlot -t 10 --summary ${DD}/sequencing_summary.txt -o ${WD}/NanoPlot_HighAccurancy_Guppy_GPU_230704

#####--------- concatenate, move, and rename the high-accuracy base called .fastq.gz files that passed the filter to the working directory
cat ${DD}/*.fastq.gz  > ${WD}/GuppyGPU_HAC_All_230230704.fastq.gz
echo "finished concatenate"

######------- FASTQC to assess quality of the concatenated sequence data
## FastQC: run on each of the data files that have 'All' to check the quality of the data
## The output from this analysis is a folder of results and a zipped file of results

fastqc ${WD}/GuppyGPU_HAC_All_230230704.fastq.gz --outdir=${WD}/FastQC

#####------- Porechop to cut adapters and barcodes
### -i input file, your basecalled fastq file
### porechop -i input_reads.fastq.gz -b output_dir

porechop -i ${WD}/GuppyGPU_HAC_All_230230704.fastq.gz   -t 10 -b ${WD}/Porechopped
echo "finished porechop"

######------- FASTQC to assess quality of the concatenated sequence data
## FastQC: run on each of the data files that have 'All' to check the quality of the data
## The output from this analysis is a folder of results and a zipped file of results

######------- FASTQC to assess quality and number of each barcode files
fastqc ${WD}/Porechopped/*  --outdir=${WD}/FastQC_postPorechop
