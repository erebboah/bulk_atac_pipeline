#!/bin/bash
#SBATCH --job-name=finalpeaks    ## Name of the job
#SBATCH -A SEYEDAM_LAB            ## account to charge 
#SBATCH -p standard               ## partition/queue name
#SBATCH --nodes=1                 ## (-N) number of nodes to use
#SBATCH --array=1-2                 ## number of tasks to launch (number of samples in the study)
#SBATCH --cpus-per-task=1         ## number of cores the job needs
#SBATCH --output=finalpeaks-%J.out ## output log file
#SBATCH --error=finalpeaks-%J.err ## error log file

expt=$1 # e.g. C2C12_MB
inpath="/share/crsp/lab/seyedam/erebboah/bulk_atac_pipeline/"
file=$inpath"samples.txt" # e.g. C2C12_ER on one line, C2C12_IR on other line, etc. Easier to add more samples for big experiments.
prefix=`head -n $SLURM_ARRAY_TASK_ID  ${file} | tail -n 1`

mkdir ${inpath}${expt}

cp ${inpath}${prefix}/peaks150/idr-output/combined.peaks-top-set.txt ${inpath}${expt}/${prefix}_150bp.txt
cp ${inpath}${prefix}/peaks500/idr-output/combined.peaks-top-set.txt ${inpath}${expt}/${prefix}_500bp.txt

# make idr output into bed files
awk -v OFS='\t' '{print $2,$3,$4}' ${inpath}${expt}/${prefix}_150bp.txt | tail -n+2 | LC_COLLATE=C sort -k1,1 -k2,2n > ${inpath}${expt}/${prefix}_150bp.bed
awk -v OFS='\t' '{print $2,$3,$4}' ${inpath}${expt}/${prefix}_500bp.txt | tail -n+2 | LC_COLLATE=C sort -k1,1 -k2,2n > ${inpath}${expt}/${prefix}_500bp.bed

rm ${inpath}${expt}/${prefix}_150bp.txt
rm ${inpath}${expt}/${prefix}_500bp.txt
