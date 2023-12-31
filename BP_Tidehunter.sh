#!/bin/bash
##---------------------------------------------------------------------
## This script is for separating telomere sequences out on MinION data (.fast5) using
##  the GPU nodes on EASLEY.
##  Brynleigh Payne. July 2023. This was modified to call the Brown Anole chromosome telomere Data
##    Thanks to Auburn Uni HPC administration for feedback to get this to work.
##    Some of the config has come from information here:
##-------------------- Parameters for running on EASLEY ----------------
#SBATCH --job-name=BP_Tidehunter_All           # job name
#SBATCH --nodes=2                         # node(s) requried for job
#SBATCH --ntasks=10                       # number of tasks across all nodes
#SBATCH --partition=general                  # name of partition to submit job
#SBATCH --time=10:00:00                    # Run time (D-HH:MM:SS)
#SBATCH --output=job-%jBP_Tidehunter_All_Test20.out          # Output file. %j is replaced with job ID
#SBATCH --error=job-%jBP_Tidehunter_All_Test20.err           # Error file. %j is replaced with job ID
#SBATCH --mail-type=ALL                   # will send email for begin,end,fail
#SBATCH --mail-user=bep0022@auburn.edu
##---------------------------------------------------------------------

#### Load Python version 3.8.6 for use of Tidehunter and version 3.10.9 for use of Snakemake.
  #### Load Tidehunter (case sensitive). 
  ##The loading needs to be completed outside of the script, just included in here for accessibility.

  module load python/anaconda/3.8.6
  conda install -c bioconda tidehunter
  #module load python/anaconda/3.10.9


##Running Tidehunter without adding parameters
  #For one .fastq file-SUCCESSFUL
#TideHunter -f 2 /hosted/biosc/SchwartzLab/Nanopore/AnoleEmbryoOldiesTelomere/BaseCalled_Guppy10_GPU_HAC_230622/pass/fastq_runid_d236a8400e4b3d9b7f03775a8afc182c26ff215e_1_0.fastq > BP_Tidehunter
  #For all of the .fastq files-SUCCESSFUL-Test 12
#TideHunter -f 2 /hosted/biosc/SchwartzLab/Nanopore/AnoleEmbryoOldiesTelomere/BaseCalled_Guppy10_GPU_HAC_230622/pass/*.fastq > BP_Tidehunter_All.out

##Adding parameters-SUCCESSFUL-Test 20
TideHunter -u -t 10 -f 2 /hosted/biosc/SchwartzLab/Nanopore/AnoleEmbryoOldiesTelomere/BaseCalled_Guppy10_GPU_HAC_230622/pass/*.fastq > BP_Tidehunter_All_params.out
