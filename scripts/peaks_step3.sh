#!/bin/bash
#SBATCH --job-name=callpeaks    ## Name of the job
#SBATCH -A SEYEDAM_LAB            ## account to charge 
#SBATCH -p standard               ## partition/queue name
#SBATCH --nodes=1                 ## (-N) number of nodes to use
#SBATCH --array=1-4               ## number of tasks to launch (number of samples)
#SBATCH --cpus-per-task=16         ## number of cores the job needs
#SBATCH --output=callpeaks-%J.out ## output log file
#SBATCH --error=callpeaks-%J.err ## error log file

inpath="/share/crsp/lab/seyedam/erebboah/bulk_atac_pipeline/"
file=$inpath"prefixes.txt" # sample prefix file containing the names for the samples you want to process
prefix=`head -n $SLURM_ARRAY_TASK_ID  ${file} | tail -n 1`

mkdir ${inpath}${prefix}/peaks
 
findPeaks ${inpath}${prefix}/mapped/homer-tags -LP .1 -poisson .1 -style factor -size 150 -minDist 50 -localSize 50000 -o ${inpath}${prefix}/peaks/${prefix}.150bp.peaks.txt

findPeaks ${inpath}${prefix}/mapped/homer-tags -LP .1 -poisson .1 -style factor -size 500 -minDist 50 -localSize 50000 -o ${inpath}${prefix}/peaks/${prefix}.500bp.peaks.txt
