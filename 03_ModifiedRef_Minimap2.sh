#!/bin/bash
#####################
#---------------Mapping genome to reference
#######################
# Be sure to change these #SBATCH parameters as appropriate for YOUR job.
#############################################################
#SBATCH --job-name=ModifiedRef_Minimap2_Test5           # job name
#SBATCH --nodes=2                                                 #node(s) requried for job
#SBATCH --ntasks=8                          # number of tasks across all nodes
#SBATCH --partition=general                 # name of partition to submit job
#SBATCH --time=05:00:00                     # Run time (D-HH:MM:SS)
#SBATCH --output=job-%j_ModifiedRef_Minimap2_Test5.out          # Output file. %j is replaced with job ID
#SBATCH --error=job-%j_ModifiedRef_Minimap2_Test5.err           # Error file. %j is replaced with job ID
#SBATCH --mail-type=ALL              # will send email for begin,end,fail
#SBATCH --mail-user=bep0022@auburn.edu
#---------------

 # /tools/guppy-6.4.6gpu/bin/guppy_basecaller
#source /opt/asn/etc/asn-bash-profiles-special/modules.sh
#module load guppy/6.4.6gpu

########## Load Modules
module load fastqc
module load python/anaconda/3.8.6   #for NanoPlot, porechop, medaka, quast
module load minimap2/2.26
module load samtools/1.17
module load bcftools/1.17

##########  Define variables and make directories
#  DD is the Guppy GPU High Accuracy basecalls; WD are the working directories in scratch; REF is for the modified reference genome (added 2KB telomeri$
DD=/hosted/biosc/SchwartzLab/Nanopore/AnoleEmbryoOldiesTelomere/BaseCalled_Guppy10_GPU_HAC_230622/pass/*.fastq
WD=/scratch/bep0022
RESULTS=${WD}/Modified_Results_HighAccuracy_Guppy_GPU_20230708
REF=/home/bep0022/ReferenceGenomes/AnoSag2.1_Scaffold1to14_telo.fasta
i=ModifiedRef_Minimap2_Test

##  make the directories in SCRATCH for holding the raw data
## -p tells it to make any upper level directories that are not there. Notice how this will also make the WD.
mkdir ${RESULTS}
mkdir ${RESULTS}/ConsensusSequences_Modified

#####----- Align to the reference using MiniMap2
## minimap2 <reference genome file> <fastq files> > <alignment.sam> -ax map-ont
##	-a  Generate CIGAR and output alignments in the SAM format. Minimap2 outputs in PAF by default.
##	-x map-ont Align noisy long reads of ~10% error rate to a reference genome. This is the default mode.
minimap2 -a -x map-ont ${REF} ${WD}/Porechopped/*.fastq.gz > ${RESULTS}/${i}.sam

samtools view -bS ${RESULTS}/${i}.sam > ${RESULTS}/${i}.bam
samtools sort ${RESULTS}/${i}.bam -o ${RESULTS}/${i}_sorted.bam
samtools flagstat ${RESULTS}/${i}_sorted.bam > ${RESULTS}/${i}.stats.txt
samtools depth ${RESULTS}/${i}_sorted.bam |  awk '{sum+=$3} END { print "Average = ",sum/NR}' > ${RESULTS}/${i}.depth.txt
samtools consensus -f fasta ${RESULTS}/${i}_sorted.bam  -o ${RESULTS}/ConsensusSequences_Modified/ModifiedRef_SAM_cons.fa
