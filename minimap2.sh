#!/bin/bash
#####################
#---------------Mapping genome to reference
#######################
# Be sure to change these #SBATCH parameters as appropriate for YOUR job.
#############################################################
#SBATCH --job-name=Minimap2_Test           # job name
#SBATCH --nodes=2                                                 #node(s) requried for job
#SBATCH --ntasks=8                          # number of tasks across all nodes
#SBATCH --partition=general                 # name of partition to submit job
#SBATCH --time=05:00:00                     # Run time (D-HH:MM:SS)
#SBATCH --output=job-%j_Minimap2_Test2.out          # Output file. %j is replaced with job ID
#SBATCH --error=job-%j_Minimap2_Test2.err           # Error file. %j is replaced with job ID
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
#  DD is the Guppy GPU High Accuracy basecalls; WD are the working directories in scratch; REF and ANNO are for the reference genomes.
DD=/hosted/biosc/SchwartzLab/Nanopore/AnoleEmbryoOldiesTelomere/BaseCalled_Guppy10_GPU_HAC_230622/pass
WD=/scratch/bep0022
RESULTS=${WD}/Results_HighAccuracy_Guppy_GPU_20230708
REF=/hosted/biosc/SchwartzLab/ReferenceGenomes/Anolis/Anolis.sagrei.Bahama__Draft2021/AnoSag2.1_gene_annotation_032122_CDS.fa

##  make the directories in SCRATCH for holding the raw data
## -p tells it to make any upper level directories that are not there. Notice how this will also make the WD.
mkdir ${WD}/ConsensusSequences
mkdir ${RESULTS}

#####----- Align to the reference using MiniMap2
## minimap2 <reference genome file> <fastq files> > <alignment.sam> -ax map-ont
##	-a  Generate CIGAR and output alignments in the SAM format. Minimap2 outputs in PAF by default.
##	-x map-ont Align noisy long reads of ~10% error rate to a reference genome. This is the default mode.
minimap2 -a -x map-ont ${REF} ${WD}/Porechopped/*.fastq.gz  > ${i}.sam

samtools view -bS ${i}.sam > ${i}.bam
samtools sort ${i}.bam -o ${i}_sorted.bam
samtools flagstat ${i}_sorted.bam > ${i}.stats.txt
samtools depth ${i}_sorted.bam |  awk '{sum+=$3} END { print "Average = ",sum/NR}' > ${i}.depth.txt
samtools consensus -f fasta ${i}_sorted.bam  -o ${WD}/ConsensusSequences/${i}_SAM_cons.fa
