#!/bin/bash
##---------------------------------------------------------------------
## This script is for separating telomere sequences out on MinION data (.fast5) using
##  the GPU nodes on EASLEY.
##  Brynleigh Payne. July 2023. This was modified to call the Brown Anole chromosome telomere Data
##    Thanks to Auburn Uni HPC administration for feedback to get this to work.
##    Some of the config has come from information here:
##-------------------- Parameters for running on EASLEY ----------------
#SBATCH --job-name=Tidehunter_Test           # job name
#SBATCH --nodes=2                         # node(s) requried for job
#SBATCH --ntasks=10                       # number of tasks across all nodes
#SBATCH --partition=general                  # name of partition to submit job
#SBATCH --time=10:00:00                    # Run time (D-HH:MM:SS)
#SBATCH --output=job-%jTidehunter_Test3.out          # Output file. %j is replaced with job ID
#SBATCH --error=job-%jTidehunter_Test3.err           # Error file. %j is replaced with job ID
#SBATCH --mail-type=ALL                   # will send email for begin,end,fail
#SBATCH --mail-user=bep0022@auburn.edu
##---------------------------------------------------------------------

##Software to install
Minimap2: /home/bep0022/bin/minimap2
TideHunter: /home/bep0022/miniconda3/bin/TideHunter
Porechop: /home/bep0022/repos/Porechop
samtools: /home/bep0022/miniconda3/bin/samtools
bedtools: /home/bep0022/miniconda3/bin/bedtools
seqkit: /home/bep0022/bin/seqkit

##References
#Must download
starts: telo_start_position_WT_2kb.txt
coverage: 2
ref: /hosted/biosc/SchwartzLab/ReferenceGenomes/Anolis/Anolis.sagrei.Bahama__Draft2021/AnoSag2.1_gene_annotation_032122_CDS.fa
outdir: /hosted/biosc/SchwartzLab/Nanopore/AnoleEmbryoOldiesTelomere/snakemake_out
reads: /hosted/biosc/SchwartzLab/Nanopore/AnoleEmbryoOldiesTelomere/BaseCalled_Guppy10_GPU_HAC_230622/pass/fastq_runid_d236a8400e4b3d9b7f03775a8afc182c26ff215e_1_0.fastq
telo_seq: TTAGGG
threads: 48

while :
do
  case "$1" in
    -i | --in ) # input Tidehunter
      in=$2
      shift 2
      ;;
    -r | --reads ) # fastq
      reads=$2
      shift 2
      ;;
    -s | --sequence ) #telomere sequence (optional : defaults to TTAGGG)
      sequence=$2
      shift 2
      ;;
    -p | --seqkit_path ) #telomere sequence (optional : defaults to TTAGGG)
      seqkit_path=$2
      shift 2
      ;;
    -h | --help ) # help message
      helpmsg=1
      shift 1
      ;;
    *) break
      ;;
  esac
done

if [ -z $in ];then
  errormsg="Error : Tidehunter output must be given"
  helpmsg=1
fi
if [ -z $reads ];then
  errormsg="Error : fastq must be given"
  helpmsg=1
fi
if [ ! -z $helpmsg ];then
  echo "./filter.sh [-i /Tidehunter/directory] [-s telomere_sequence [string]] /path/to/basecalled/sample
  Filter Tidehunter fastq for reads with telomere and adapter seqence
  ;-i/--in;;Tidehunter input file
  ;-r/--reads;;Fastq
  ;-s/--sequence;;Telomere sequence (default TTAGGG)
  ;-p/--seqkit_path;;path to seqkit executable
  ;-h/--help;;show help message and exit"|\
    tr ";" "\t"
  echo $errormsg
  exit
fi
if [ -z $sequence ];then
  sequence=TTAGGG
fi

mkdir -p ./tmp
out=./tmp

grep $sequence $in | cut -f1 | sed 's/^/@/' > ${out}/ID.tmp
#num=$(grep $sequence $in | cut -f1 | wc -l)
#echo $num reads with $sequence

##Generate a filtered .fastq with only reads with telomere
awk 'NR==FNR{a[$0];next}$1 in a {x=NR+3}(NR<=x){print} ' ${out}/ID.tmp ${reads} > ${out}/filtered.fastq
##Filter out reads < 5kb
#awk 'BEGIN {FS = "\t" ; OFS = "\n"} {header = $0 ; getline seq ; getline qheader ; getline qseq ; if (length(seq) >= 5000) {print header, seq, qheader, qseq}}' < ${out}/telomere.fas$
##Filter for reads with more than 10 As in the last 100 bp
$seqkit_path grep -s -R -100:-1 -r -p AAAAAAAAAA ${out}/filtered.fastq > ${out}/As.fastq
##Filter for reads with more then 10 Ts in the first 100 bp
$seqkit_path grep -s -R 1:100 -r -p TTTTTTTTTT ${out}/filtered.fastq > ${out}/Ts.fastq
##Combine the reads that are tailed
cat ${out}/As.fastq ${out}/Ts.fastq
rm -r $out
