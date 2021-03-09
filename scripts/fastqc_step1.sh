#!/bin/bash
#SBATCH --job-name=checkreads    ## Name of the job
#SBATCH -A SEYEDAM_LAB            ## account to charge 
#SBATCH -p standard               ## partition/queue name
#SBATCH --nodes=1                 ## (-N) number of nodes to use
#SBATCH --array=1-5               ## number of tasks to launch (number of samples)
#SBATCH --cpus-per-task=1         ## number of cores the job needs
#SBATCH --output=checkreads-%J.out ## output log file
#SBATCH --error=checkreads-%J.err ## error log file

inpath="/share/crsp/lab/seyedam/share/newATAC/"
file=$inpath"prefixes_all.txt" # sample prefix file containing the names for the samples you want to process
prefix=`head -n $SLURM_ARRAY_TASK_ID  ${file} | tail -n 1`
 
module load fastqc/0.11.9

mkdir ${inpath}${prefix}
mkdir ${inpath}${prefix}/fastqc

fastqc ${inpath}data/${prefix}_R1.fastq.gz ${inpath}data/${prefix}_R2.fastq.gz -o ${inpath}${prefix}/fastqc/
