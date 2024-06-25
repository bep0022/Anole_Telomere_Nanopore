#!/bin/bash
#####################
#---------------This is a bash script for calling python scripts provided for finding the telomeres
#######################
# Be sure to change these #SBATCH parameters as appropriate for YOUR job.
#############################################################
#SBATCH --job-name=TelFinder           # job name
#SBATCH --nodes=2                                                 #node(s) requried for job
#SBATCH --ntasks=8                          # number of tasks across all nodes
#SBATCH --partition=general                 # name of partition to submit job
#SBATCH --time=08:00:00                     # Run time (D-HH:MM:SS)
#SBATCH --output=job-%j_TelFinder31.out          # Output file. %j is replaced with job ID
#SBATCH --error=job-%j_TelFinder31.err           # Error file. %j is replaced with job ID
#SBATCH --mail-type=ALL              # will send email for begin,end,fail
#SBATCH --mail-user=bep0022@auburn.edu
#---------------

##########  Define variables and make directories

WD=/scratch/bep0022/                                                                                            ##This is the working directory
DD=/hosted/biosc/SchwartzLab/Nanopore/AnoleEmbryoOldiesTelomere/BaseCalled_Guppy10_GPU_HAC_230622/pass/         ##This is where the .fastq files are located (input)
RESULTSD=/home/bep0022/TelFinder_AnoleTrial_BAM/                                                                ##This is where the results you want to be brought to
#REFD=/hosted/biosc/SchwartzLab/ReferenceGenomes/Anolis                                                         ##This is where the reference genome is located

#REF=GCA_037176765.1_rAnoSag1.mat_genomic.fna

## Make the directories and all subdirectories defined by the variables above
mkdir -p $RESULTSD
#mkdir -p /home/bep0022/TelFinder_fastq1/
#mkdir -p /home/bep0022/TelFinder_fastq2/
#mkdir -p /home/bep0022/TelFinder_fastq3/
#mkdir -p /home/bep0022/TelFinder_fastq4/
#mkdir -p /home/bep0022/TelFinder_fastq5/
#mkdir -p /home/bep0022/TelFinder_fastq6/

##########


########## Load modules for running the python scripts

module load python/anaconda/3.8.6
module load hisat2/2.2.1
module load stringtie/2.1.6
module load gcc/9.3.0
module load samtools/1.19
module load bcftools/1.17
module load gffreader/12.7


########## Run Python scripts

####### example; Q:format=fastq(paired)
##python3 TelFinder.py -f fq2 -1 example/fastq/DRR237583_1.fastq -2 example/fastq/DRR237583_2.fastaq -ref example/fastq/AFM_genomic.fna -o example/fastq/ -s AFMfq -t 80

###FindKmer
#python3 FindKmer.py -f fq2 -1 fastq_runid_d236a8400e4b3d9b7f03775a8afc182c26ff215e_0_0.fastq -2 fastq_runid_d236a8400e4b3d9b7f03775a8afc182c26ff215e_1_0.fastq -s Anolissagrei -ref GCA_037176765.1_rAnoSa$

###statMut
#python3 statMut.py -f fq2 -1 fastq_runid_d236a8400e4b3d9b7f03775a8afc182c26ff215e_0_0.fastq -2 fastq_runid_d236a8400e4b3d9b7f03775a8afc182c26ff215e_1_0.fastq -s Anolissagrei -ref GCA_037176765.1_rAnoSag$

###score
#python3 score.py -f fq2 -1 fastq_runid_d236a8400e4b3d9b7f03775a8afc182c26ff215e_0_0.fastq -2 fastq_runid_d236a8400e4b3d9b7f03775a8afc182c26ff215e_1_0.fastq -s Anolissagrei -ref GCA_037176765.1_rAnoSag1.$

###TelFinder
#python3 TelFinder.py -f fq2 -1 fastq_runid_d236a8400e4b3d9b7f03775a8afc182c26ff215e_0_0.fastq -2 fastq_runid_d236a8400e4b3d9b7f03775a8afc182c26ff215e_1_0.fastq -s Anolissagrei -ref GCA_037176765.1_rAnoS$



####### Q:format=fasta
##python3 TelFinder.py -f fa -inf example/fasta/NCR.fna -s NCR -o example/fasta -e left
##python3 TelFinder.py -f fa -inf example/fasta/NCR.fna -s NCR -o example/fasta -e right

####Anole Mat Genome
###FindKmer
#python3 FindKmer.py -f fa -inf GCA_037176765.1_rAnoSag1.mat_genomic.fna -s Anolissagrei -o /home/bep0022/TelFinder_AnoleMatGenome2/ -e left
#python3 FindKmer.py -f fa -inf GCA_037176765.1_rAnoSag1.mat_genomic.fna -s Anolissagrei -o /home/bep0022/TelFinder_AnoleMatGenome2/ -e right

###statMut
#python3 statMut.py -f fa -inf GCA_037176765.1_rAnoSag1.mat_genomic.fna -s Anolissagrei -o /home/bep0022/TelFinder_AnoleMatGenome2/ -e left
#python3 statMut.py -f fa -inf GCA_037176765.1_rAnoSag1.mat_genomic.fna -s Anolissagrei -o /home/bep0022/TelFinder_AnoleMatGenome2/ -e right

###score
#python3 score.py -f fa -inf GCA_037176765.1_rAnoSag1.mat_genomic.fna -s Anolissagrei -o /home/bep0022/TelFinder_AnoleMatGenome2/ -e left
#python3 score.py -f fa -inf GCA_037176765.1_rAnoSag1.mat_genomic.fna -s Anolissagrei -o /home/bep0022/TelFinder_AnoleMatGenome2/ -e right

###TelFinder
#python3 TelFinder.py -f fa -inf GCA_037176765.1_rAnoSag1.mat_genomic.fna -s Anolissagrei -o /home/bep0022/TelFinder_AnoleMatGenome2/ -e left
#python3 TelFinder.py -f fa -inf GCA_037176765.1_rAnoSag1.mat_genomic.fna -s Anolissagrei -o /home/bep0022/TelFinder_AnoleMatGenome2/ -e right


####Anole Trial Mapped Reads to Mat Genome
###FindKmer
#python3 FindKmer.py -f fa -inf newgenome_minimap2_SAM_cons.fa -s Anolissagrei -o /home/bep0022/TelFinder_AnoleTrialMappedReads/ $
#python3 FindKmer.py -f fa -inf newgenome_minimap2_SAM_cons.fa -s Anolissagrei -o /home/bep0022/TelFinder_AnoleTrialMappedReads/ $

###statMut
#python3 statMut.py -f fa -inf newgenome_minimap2_SAM_cons.fa -s Anolissagrei -o /home/bep0022/TelFinder_AnoleTrialMappedReads/ -$
#python3 statMut.py -f fa -inf newgenome_minimap2_SAM_cons.fa -s Anolissagrei -o /home/bep0022/TelFinder_AnoleTrialMappedReads/ -$

###score
#python3 score.py -f fa -inf newgenome_minimap2_SAM_cons.fa -s Anolissagrei -o /home/bep0022/TelFinder_AnoleTrialMappedReads/ -e $
#python3 score.py -f fa -inf newgenome_minimap2_SAM_cons.fa -s Anolissagrei -o /home/bep0022/TelFinder_AnoleTrialMappedReads/ -e $

###TelFinder
#python3 TelFinder.py -f fa -inf newgenome_minimap2_SAM_cons.fa -s Anolissagrei -o /home/bep0022/TelFinder_AnoleTrialMappedReads/$
#python3 TelFinder.py -f fa -inf newgenome_minimap2_SAM_cons.fa -s Anolissagrei -o /home/bep0022/TelFinder_AnoleTrialMappedReads/$


####Anole Trial Mapped Reads to Mat Genome (.bam)
###FindKmer
python3 FindKmer.py -f bam -inf newgenome_minimap2.bam -s Anolissagrei -o /home/bep0022/TelFinder_AnoleTrial_BAM/ -e left
python3 FindKmer.py -f bam -inf newgenome_minimap2.bam -s Anolissagrei -o /home/bep0022/TelFinder_AnoleTrial_BAM/ -e right

###statMut
python3 statMut.py -f bam -inf newgenome_minimap2.bam -s Anolissagrei -o /home/bep0022/TelFinder_AnoleTrial_BAM/ -e left
python3 statMut.py -f bam -inf newgenome_minimap2.bam -s Anolissagrei -o /home/bep0022/TelFinder_AnoleTrial_BAM/ -e right

###score
python3 score.py -f bam -inf newgenome_minimap2.bam -s Anolissagrei -o /home/bep0022/TelFinder_AnoleTrial_BAM/ -e left
python3 score.py -f bam -inf newgenome_minimap2.bam -s Anolissagrei -o /home/bep0022/TelFinder_AnoleTrial_BAM/ -e right

###TelFinder
python3 TelFinder.py -f bam -inf newgenome_minimap2.bam -s Anolissagrei -o /home/bep0022/TelFinder_AnoleTrial_BAM/ -e left
python3 TelFinder.py -f bam -inf newgenome_minimap2.bam -s Anolissagrei -o /home/bep0022/TelFinder_AnoleTrial_BAM/ -e right


####Human Genome
###FindKmer
#python3 FindKmer.py -f fa -inf GCA_000001405.29_GRCh38.p14_genomic.fna -s Homosapiens -o /home/bep0022/TelFinder_HumanGenomeRef/

###statMut
#python3 statMut.py -f fa -inf GCA_000001405.29_GRCh38.p14_genomic.fna -s Homosapiens -o /home/bep0022/TelFinder_HumanGenomeRef/

###score
#python3 score.py -f fa -inf GCA_000001405.29_GRCh38.p14_genomic.fna -s Homosapiens -o /home/bep0022/TelFinder_HumanGenomeRef/

###TelFinder
#python3 TelFinder.py -f fa -inf GCA_000001405.29_GRCh38.p14_genomic.fna -s Homosapiens -o /home/bep0022/TelFinder_HumanGenomeRef/


####Zebra Finch Genome
###FindKmer
#python3 FindKmer.py -f fa -inf GCA_003957565.4_bTaeGut1.4.pri_genomic.fna -s Taeniopygiaguttata -o /home/bep0022/TelFinder_Zebra$
#python3 FindKmer.py -f fa -inf GCA_003957565.4_bTaeGut1.4.pri_genomic.fna -s Taeniopygiaguttata -o /home/bep0022/TelFinder_Zebra$

###statMut
#python3 statMut.py -f fa -inf GCA_003957565.4_bTaeGut1.4.pri_genomic.fna -s Taeniopygiaguttata -o /home/bep0022/TelFinder_ZebraF$
#python3 statMut.py -f fa -inf GCA_003957565.4_bTaeGut1.4.pri_genomic.fna -s Taeniopygiaguttata -o /home/bep0022/TelFinder_ZebraF$

###score
#python3 score.py -f fa -inf GCA_003957565.4_bTaeGut1.4.pri_genomic.fna -s Taeniopygiaguttata -o /home/bep0022/TelFinder_ZebraFin$
#python3 score.py -f fa -inf GCA_003957565.4_bTaeGut1.4.pri_genomic.fna -s Taeniopygiaguttata -o /home/bep0022/TelFinder_ZebraFin$

###TelFinder
#python3 TelFinder.py -f fa -inf GCA_003957565.4_bTaeGut1.4.pri_genomic.fna -s Taeniopygiaguttata -o /home/bep0022/TelFinder_Zebr$
#python3 TelFinder.py -f fa -inf GCA_003957565.4_bTaeGut1.4.pri_genomic.fna -s Taeniopygiaguttata -o /home/bep0022/TelFinder_Zebr$



####### Q:format=fastq(single)
##python3 TelFinder.py -f fq1 -inf example/fastq/DD.fastq -ref example/fastq/AFM_genomic.fna -o example/fastq/ -s AFMfq -t 80

#### Anole trial raw reads
#### 0_0.fastq
###FindKmer
#python3 FindKmer.py -f fq1 -inf fastq_runid_d236a8400e4b3d9b7f03775a8afc182c26ff215e_0_0.fastq -ref GCA_037176765.1_rAnoSag1.mat$

###statMut
#python3 statMut.py -f fq1 -inf fastq_runid_d236a8400e4b3d9b7f03775a8afc182c26ff215e_0_0.fastq -ref GCA_037176765.1_rAnoSag1.mat_$

###score
#python3 score.py -f fq1 -inf fastq_runid_d236a8400e4b3d9b7f03775a8afc182c26ff215e_0_0.fastq -ref GCA_037176765.1_rAnoSag1.mat_ge$

###TelFinder
#python3 TelFinder.py -f fq1 -inf fastq_runid_d236a8400e4b3d9b7f03775a8afc182c26ff215e_0_0.fastq -ref GCA_037176765.1_rAnoSag1.ma$


#### 1_0.fastq
###FindKmer
#python3 FindKmer.py -f fq1 -inf fastq_runid_d236a8400e4b3d9b7f03775a8afc182c26ff215e_1_0.fastq -ref GCA_037176765.1_rAnoSag1.mat$

###statMut
#python3 statMut.py -f fq1 -inf fastq_runid_d236a8400e4b3d9b7f03775a8afc182c26ff215e_1_0.fastq -ref GCA_037176765.1_rAnoSag1.mat_$

###score
#python3 score.py -f fq1 -inf fastq_runid_d236a8400e4b3d9b7f03775a8afc182c26ff215e_1_0.fastq -ref GCA_037176765.1_rAnoSag1.mat_ge$

###TelFinder
#python3 TelFinder.py -f fq1 -inf fastq_runid_d236a8400e4b3d9b7f03775a8afc182c26ff215e_1_0.fastq -ref GCA_037176765.1_rAnoSag1.ma$


#### 2_0.fastq
###FindKmer
#python3 FindKmer.py -f fq1 -inf fastq_runid_d236a8400e4b3d9b7f03775a8afc182c26ff215e_2_0.fastq -ref GCA_037176765.1_rAnoSag1.mat$

###statMut
#python3 statMut.py -f fq1 -inf fastq_runid_d236a8400e4b3d9b7f03775a8afc182c26ff215e_2_0.fastq -ref GCA_037176765.1_rAnoSag1.mat_$

###score
#python3 score.py -f fq1 -inf fastq_runid_d236a8400e4b3d9b7f03775a8afc182c26ff215e_2_0.fastq -ref GCA_037176765.1_rAnoSag1.mat_ge$

###TelFinder
#python3 TelFinder.py -f fq1 -inf fastq_runid_d236a8400e4b3d9b7f03775a8afc182c26ff215e_2_0.fastq -ref GCA_037176765.1_rAnoSag1.ma$


#### 3_0.fastq
###FindKmer
#python3 FindKmer.py -f fq1 -inf fastq_runid_d236a8400e4b3d9b7f03775a8afc182c26ff215e_3_0.fastq -ref GCA_037176765.1_rAnoSag1.mat$

###statMut
#python3 statMut.py -f fq1 -inf fastq_runid_d236a8400e4b3d9b7f03775a8afc182c26ff215e_3_0.fastq -ref GCA_037176765.1_rAnoSag1.mat_$

###score
#python3 score.py -f fq1 -inf fastq_runid_d236a8400e4b3d9b7f03775a8afc182c26ff215e_3_0.fastq -ref GCA_037176765.1_rAnoSag1.mat_ge$

###TelFinder
#python3 TelFinder.py -f fq1 -inf fastq_runid_d236a8400e4b3d9b7f03775a8afc182c26ff215e_3_0.fastq -ref GCA_037176765.1_rAnoSag1.ma$


#### 4_0.fastq
###FindKmer
#python3 FindKmer.py -f fq1 -inf fastq_runid_d236a8400e4b3d9b7f03775a8afc182c26ff215e_4_0.fastq -ref GCA_037176765.1_rAnoSag1.mat$

###statMut
#python3 statMut.py -f fq1 -inf fastq_runid_d236a8400e4b3d9b7f03775a8afc182c26ff215e_4_0.fastq -ref GCA_037176765.1_rAnoSag1.mat_$

###score
#python3 score.py -f fq1 -inf fastq_runid_d236a8400e4b3d9b7f03775a8afc182c26ff215e_4_0.fastq -ref GCA_037176765.1_rAnoSag1.mat_ge$

###TelFinder
#python3 TelFinder.py -f fq1 -inf fastq_runid_d236a8400e4b3d9b7f03775a8afc182c26ff215e_4_0.fastq -ref GCA_037176765.1_rAnoSag1.ma$


#### 5_0.fastq
###FindKmer
#python3 FindKmer.py -f fq1 -inf fastq_runid_d236a8400e4b3d9b7f03775a8afc182c26ff215e_5_0.fastq -ref GCA_037176765.1_rAnoSag1.mat$

###statMut
#python3 statMut.py -f fq1 -inf fastq_runid_d236a8400e4b3d9b7f03775a8afc182c26ff215e_5_0.fastq -ref GCA_037176765.1_rAnoSag1.mat_$

###score
#python3 score.py -f fq1 -inf fastq_runid_d236a8400e4b3d9b7f03775a8afc182c26ff215e_5_0.fastq -ref GCA_037176765.1_rAnoSag1.mat_ge$

###TelFinder
#python3 TelFinder.py -f fq1 -inf fastq_runid_d236a8400e4b3d9b7f03775a8afc182c26ff215e_5_0.fastq -ref GCA_037176765.1_rAnoSag1.ma$


#### 6_0.fastq
###FindKmer
#python3 FindKmer.py -f fq1 -inf fastq_runid_d236a8400e4b3d9b7f03775a8afc182c26ff215e_6_0.fastq -ref GCA_037176765.1_rAnoSag1.mat$

###statMut
#python3 statMut.py -f fq1 -inf fastq_runid_d236a8400e4b3d9b7f03775a8afc182c26ff215e_6_0.fastq -ref GCA_037176765.1_rAnoSag1.mat_$

###score
#python3 score.py -f fq1 -inf fastq_runid_d236a8400e4b3d9b7f03775a8afc182c26ff215e_6_0.fastq -ref GCA_037176765.1_rAnoSag1.mat_ge$

###TelFinder
#python3 TelFinder.py -f fq1 -inf fastq_runid_d236a8400e4b3d9b7f03775a8afc182c26ff215e_6_0.fastq -ref GCA_037176765.1_rAnoSag1.ma$

