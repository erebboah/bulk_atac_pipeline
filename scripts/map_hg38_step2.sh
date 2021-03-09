#!/bin/bash
#SBATCH --job-name=maptohg38    ## Name of the job
#SBATCH -A SEYEDAM_LAB            ## account to charge 
#SBATCH -p standard               ## partition/queue name
#SBATCH --nodes=1                 ## (-N) number of nodes to use
#SBATCH --array=1-2               ## number of tasks to launch (number of samples)
#SBATCH --cpus-per-task=16         ## number of cores the job needs
#SBATCH --output=maptohg38-%J.out ## output log file
#SBATCH --error=maptohg38-%J.err ## error log file

inpath="/share/crsp/lab/seyedam/share/newATAC/"
file=$inpath"prefixes_hg38.txt" # sample prefix file containing the names for the samples you want to process
ref="/share/crsp/lab/seyedam/share/newATAC/ref/hg38/hg38"
refM="/share/crsp/lab/seyedam/share/newATAC/ref/chrM/chrM"
prefix=`head -n $SLURM_ARRAY_TASK_ID  ${file} | tail -n 1`
 
module load bowtie2/2.4.1
module load samtools/1.10
module load picard-tools/1.87
module load java/1.8.0
module load python/3.8.0
module load ucsc-tools/v393

mkdir ${inpath}${prefix}
mkdir ${inpath}${prefix}/mapped/

# First map the reads to mitochondrial chrom only
bowtie2 --very-sensitive --no-discordant -X 3000 -k 6 -p 16 --un-conc ${inpath}${prefix}/mapped/${prefix}_unaligned.fastq.gz -x ${refM} -1 ${inpath}data/${prefix}_R1.fastq.gz -2 ${inpath}data/${prefix}_R2.fastq.gz | samtools view -Sb - > ${inpath}${prefix}/mapped/${prefix}_chrM.bam

# Next step, map unaligned reads to the hg38 genome
bowtie2 --very-sensitive -X 3000 --no-discordant -k 6 -p 16 -x ${ref} -1 ${inpath}${prefix}/mapped/${prefix}_unaligned.fastq.1.gz -2 ${inpath}${prefix}/mapped/${prefix}_unaligned.fastq.2.gz | samtools view -Sb - > ${inpath}${prefix}/mapped/${prefix}.bam

# Sort the merged BAM file
java -jar /opt/apps/picard-tools/1.87/SortSam.jar \
I=${inpath}${prefix}/mapped/${prefix}.bam \
O=${inpath}${prefix}/mapped/${prefix}.sort.bam \
SORT_ORDER=coordinate

# Remove PCR duplicates
java -jar  /opt/apps/picard-tools/1.87/MarkDuplicates.jar \
I=${inpath}${prefix}/mapped/${prefix}.sort.bam \
O=${inpath}${prefix}/mapped/${prefix}.sort.nodup.bam \
M=${inpath}${prefix}/mapped/${prefix}.duplicates_metric.txt \
REMOVE_DUPLICATES=true

# Adjust the read start sites to represent the center of the transposon binding event. 
# Transposon binds as a dimer and inserts two adapters separated by 9 bps. 
# Reads aligning to the + strand are offset by +4 bps, and reads aligning to the – strand are offset −5 bps.

samtools view ${inpath}${prefix}/mapped/${prefix}.sort.nodup.bam | python shift.reads.py ${inpath}${prefix}/mapped/${prefix}_shifted_reads.sam

# Make a header
samtools view -H ${inpath}${prefix}/mapped/${prefix}.sort.nodup.bam > ${inpath}${prefix}/mapped/${prefix}.sort.nodup.header.sam
cat ${inpath}${prefix}/mapped/${prefix}.sort.nodup.header.sam ${inpath}${prefix}/mapped/${prefix}_shifted_reads.sam | samtools view -Sb - > ${inpath}${prefix}/mapped/${prefix}_shifted_reads.bam
rm ${inpath}${prefix}/mapped/${prefix}_shifted_reads.sam

# Sort BAM file
samtools sort ${inpath}${prefix}/mapped/${prefix}_shifted_reads.bam -o ${inpath}${prefix}/mapped/${prefix}_shifted_reads_sorted.bam
rm ${inpath}${prefix}/mapped/${prefix}_shifted_reads.bam

# Build BAM index
samtools index ${inpath}${prefix}/mapped/${prefix}_shifted_reads_sorted.bam

# Make Homer tag directory
makeTagDirectory ${inpath}${prefix}/mapped/homer-tags ${inpath}${prefix}/mapped/${prefix}_shifted_reads_sorted.bam -format sam

# Generate bigWig file - for initial visualization of the data in the genome browser
bamCoverage --bam ${inpath}${prefix}/mapped/${prefix}_shifted_reads_sorted.bam -o ${inpath}${prefix}/mapped/${prefix}.bigWig --binSize 30 -e -of bigwig --scaleFactor 1.0 --normalizeUsing RPKM -p 16
