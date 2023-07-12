#!/bin/bash
##---------------------------------------------------------------------
## This script is for doing High Accuracy Basecalling on MinION data (.fast5) using
##  the GPU nodes on EASLEY. Run either on gpu4 or gpu2 on Easley. It will not run on nova_gpu.
##  Tonia Schwartz. June 2023. This was written to call the Seedbeetle Mitochondrial Genome Data
##    Thanks to Gabrielle de la Silva, Dasia Simpson, and Auburn Uni HPC administration for feedback to get this to work.
##    Some of the config has come from information here: https://hackmd.io/@Miles/S12SKP115.
##-------------------- Parameters for running on EASLEY ----------------
#SBATCH --job-name=Guppy10_SB           # job name
#SBATCH --nodes=1                       # node(s) requried for job
#SBATCH --ntasks=10                     # number of tasks across all nodes
#SBATCH --partition=gpu4                  # name of partition to submit job
#SBATCH --time=5:00:00                   # Run time (D-HH:MM:SS)
#SBATCH --output=job-%j_Guppy10.out          # Output file. %j is replaced with job ID
#SBATCH --error=job-%j_Guppy10.err           # Error file. %j is replaced with job ID
#SBATCH --mail-type=ALL              # will send email for begin,end,fail
#SBATCH --mail-user=bep0022@auburn.edu
#SBATCH --gres=gpu:tesla:1
##-------------------------------------------------------------------------------------
  #### /tools/guppy-6.4.6gpu/bin/guppy_basecaller
  #### source /opt/asn/etc/asn-bash-profiles-special/modules.sh

module load guppy/6.4.6gpu


##----- Define variable for directories.
  ## This is where the directory of fast5 data files.

SOURCE=/hosted/biosc/SchwartzLab/Nanopore/AnoleEmbryoOldiesTelomere/20230605_2000_MN39112_FAS96157_fb416b6a
OUTDIR=/hosted/biosc/SchwartzLab/Nanopore/AnoleEmbryoOldiesTelomere/BaseCalled_Guppy10_GPU_HAC_230622

##------ Make the directory for my results in my home folder
   ## the -p will  Create intermediate directories as required.

mkdir -p $OUTDIR

##------- Basecall the data using Guppy GPU
  ### During basecalling, the samples can be demultiplexed,
  ## i.e. separated back into folders containing only reads from each sample, according to which barcode is detected.
  ## You can do this with Guppy, by telling it which barcode kit you used, via the “barcode_kits” variable
  ## If you don't use this, all the resulting .fastq.gz files will be put in one folder.
  ## guppy-6.4.6gpu/bin/guppy_basecaller

guppy_basecaller \
        --cpu_threads_per_caller 1    --num_callers 10 \
        --input_path ${SOURCE}        --save_path ${OUTDIR}  \
        --recursive \
              --gpu_runners_per_device 6    --chunks_per_runner 1024 -x auto \
        --flowcell FLO-MIN106         --kit SQK-NBD112-24 \
        --compress_fastq

## Other Options
        #--barcode_kits "SQK-NBD112-24" \
        #--config /tools/guppy-6.4.6gpu/data/dna_r10_450bps_hac.cfg   \
        #--align_ref /hosted/biosc/SchwartzLab/Nanopore/SeedBeetle_MitoVariation_1/SeedBeetleMtGenome_MF960125.fasta
        #--barcode_kits "SQK-NBD112-24"  \
        #--trim_barcodes \ This didn't work
        #--nested_output_folder \
        #--flowcell FLO-MIN112        --kit SQK-NBD112.24

