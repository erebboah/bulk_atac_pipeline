#!/bin/bash
#SBATCH --job-name=doidr    ## Name of the job
#SBATCH -A SEYEDAM_LAB            ## account to charge 
#SBATCH -p standard               ## partition/queue name
#SBATCH --nodes=1                 ## (-N) number of nodes to use
#SBATCH --cpus-per-task=16         ## number of cores the job needs
#SBATCH --output=doidr-%J.out ## output log file
#SBATCH --error=doidr-%J.err ## error log file

rep1=$1 # rep name, e.g. C2C12_MB_ER1
rep2=$2 # e.g. C2C12_MB_ER2
sample=$3 #e.g. C2C12_MB_ER
inpath="/share/crsp/lab/seyedam/erebboah/bulk_atac_pipeline/"
prefix=${sample}

mkdir ${inpath}${prefix}
mkdir ${inpath}${prefix}/peaks150
mkdir ${inpath}${prefix}/peaks150/pseudoreps
mkdir ${inpath}${prefix}/peaks150/pseudoreps/individual
mkdir ${inpath}${prefix}/peaks150/pseudoreps/pooled
mkdir ${inpath}${prefix}/peaks150/replicates
mkdir ${inpath}${prefix}/peaks150/idr-output
mkdir ${inpath}${prefix}/peaks500
mkdir ${inpath}${prefix}/peaks500/pseudoreps
mkdir ${inpath}${prefix}/peaks500/pseudoreps/individual
mkdir ${inpath}${prefix}/peaks500/pseudoreps/pooled
mkdir ${inpath}${prefix}/peaks500/replicates
mkdir ${inpath}${prefix}/peaks500/idr-output

module load R/4.0.2

# Copy the peaks for each replicate into one folder
cp ${inpath}${rep1}/peaks/${rep1}.150bp.peaks.txt ${inpath}${prefix}/peaks150/${rep1}.150bp.peaks.txt
cp ${inpath}${rep2}/peaks/${rep2}.150bp.peaks.txt ${inpath}${prefix}/peaks150/${rep2}.150bp.peaks.txt

cp ${inpath}${rep1}/peaks/${rep1}.500bp.peaks.txt ${inpath}${prefix}/peaks500/${rep1}.500bp.peaks.txt
cp ${inpath}${rep2}/peaks/${rep2}.500bp.peaks.txt ${inpath}${prefix}/peaks500/${rep2}.500bp.peaks.txt

# Make a combined tag directory for all reps, find peaks for pooled combined
makeTagDirectory ${inpath}${prefix}/peaks150/combined -d ${inpath}${rep1}/mapped/${rep1}_homer-tags ${inpath}${rep2}/mapped/${rep2}_homer-tags
makeTagDirectory ${inpath}${prefix}/peaks500/combined -d ${inpath}${rep1}/mapped/${rep1}_homer-tags ${inpath}${rep2}/mapped/${rep2}_homer-tags

findPeaks ${inpath}${prefix}/peaks150/combined -LP .1 -poisson .1 -style factor -size 150 -minDist 50 -localSize 50000 -o ${inpath}${prefix}/peaks150/combined/combined.peaks.txt
findPeaks ${inpath}${prefix}/peaks500/combined -LP .1 -poisson .1 -style factor -size 500 -minDist 50 -localSize 50000 -o ${inpath}${prefix}/peaks500/combined/combined.peaks.txt

# Create pseudoreps for individual reps 
python /share/crsp/lab/seyedam/erebboah/bulk_atac_pipeline/scripts/run_idr.py pseudoreplicate -d ${inpath}${rep1}/mapped/${rep1}_homer-tags ${inpath}${rep2}/mapped/${rep2}_homer-tags -o ${inpath}${prefix}/peaks150/pseudoreps/individual
python /share/crsp/lab/seyedam/erebboah/bulk_atac_pipeline/scripts/run_idr.py pseudoreplicate -d ${inpath}${rep1}/mapped/${rep1}_homer-tags ${inpath}${rep2}/mapped/${rep2}_homer-tags -o ${inpath}${prefix}/peaks500/pseudoreps/individual

# Create pseudoreps for pooled
python /share/crsp/lab/seyedam/erebboah/bulk_atac_pipeline/scripts/run_idr.py pseudoreplicate -d ${inpath}${prefix}/peaks150/combined -o ${inpath}${prefix}/peaks150/pseudoreps/pooled
python /share/crsp/lab/seyedam/erebboah/bulk_atac_pipeline/scripts/run_idr.py pseudoreplicate -d ${inpath}${prefix}/peaks500/combined -o ${inpath}${prefix}/peaks500/pseudoreps/pooled

# Call peaks for individual pseudoreps
cd ${inpath}${prefix}/peaks150/pseudoreps/individual
for f in *
	do
	findPeaks $f -LP .1 -poisson .1 -style factor -size 150 -minDist 50 -localSize 50000 -o ${f}_peaks.txt
	done

cd ${inpath}${prefix}/peaks500/pseudoreps/individual
for f in *
        do
        findPeaks $f -LP .1 -poisson .1 -style factor -size 500 -minDist 50 -localSize 50000 -o ${f}_peaks.txt
        done

# Call peaks for combined pseudoreps
cd ${inpath}${prefix}/peaks150/pseudoreps/pooled
for f in *
	do
	findPeaks $f -LP .1 -poisson .1 -style factor -size 150 -minDist 50 -localSize 50000 -o ${f}_peaks.txt
	done

cd ${inpath}${prefix}/peaks500/pseudoreps/pooled
for f in *
        do
        findPeaks $f -LP .1 -poisson .1 -style factor -size 500 -minDist 50 -localSize 50000 -o ${f}_peaks.txt
        done

# Finally run IDR 
# Copy peaks called by homer for original reps to replicates folder
cp ${inpath}${prefix}/peaks150/*.txt ${inpath}${prefix}/peaks150/replicates
python /share/crsp/lab/seyedam/erebboah/bulk_atac_pipeline/scripts/run_idr.py idr \
	-p ${inpath}${prefix}/peaks150/replicates/*.txt \
	-pr ${inpath}${prefix}/peaks150/pseudoreps/individual/*.txt \
	-ppr ${inpath}${prefix}/peaks150/pseudoreps/pooled/*.txt \
	--pooled_peaks ${inpath}${prefix}/peaks150/combined/combined.peaks.txt \
	-o ${inpath}${prefix}/peaks150/idr-output \
	--threshold .01

cp ${inpath}${prefix}/peaks500/*.txt ${inpath}${prefix}/peaks500/replicates
python /share/crsp/lab/seyedam/erebboah/bulk_atac_pipeline/scripts/run_idr.py idr \
        -p ${inpath}${prefix}/peaks500/replicates/*.txt \
        -pr ${inpath}${prefix}/peaks500/pseudoreps/individual/*.txt \
        -ppr ${inpath}${prefix}/peaks500/pseudoreps/pooled/*.txt \
        --pooled_peaks ${inpath}${prefix}/peaks500/combined/combined.peaks.txt \
        -o ${inpath}${prefix}/peaks500/idr-output \
        --threshold .01
