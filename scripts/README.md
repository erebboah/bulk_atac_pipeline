### fastqc_step1.sh 
Requires fastqc/0.11.9 and paired reads with `_R1.fastq.gz` and `_R2.fastq.gz` file extensions.

Takes in a file of sample names (`prefixes_all.txt`) and launches an array of jobs, one per sample. Adjust number of jobs with `--array=1-5`.

Makes a sample folder in the main directory (`/share/crsp/lab/seyedam/share/newATAC/`) that will store results.

### map_hg38_step2.sh and map_mm10_step2.sh
Requires bowtie2/2.4.1, samtools/1.10, picard-tools/1.87, java/1.8.0, python/3.8.0, and ucsc-tools/v393, and paired reads with `_R1.fastq.gz` and `_R2.fastq.gz` file extensions.

Does NOT require `fastqc_step1.sh` to be run first.

First maps to the mitochondrial chromosome to remove mitochondrial reads, resulting in `sample_chrM.bam`.
```
bowtie2 --very-sensitive --no-discordant -X 3000 -k 6 -p 16 --un-conc sample_unaligned.fastq.gz
```

The unaligned reads are used to map to the genome, resulting in `sample.bam`.
```
bowtie2 --very-sensitive -X 3000 --no-discordant -k 6 -p 16 -x genome -1 sample_unaligned.fastq.1.gz -2 sample_unaligned.fastq.2.gz
```

The resulting bam file is sorted (`sample.sort.bam`) and PCR duplicates removed with `picard-tools` (`sample.sort.nodup.bam`). 

A custom script is used to shift reads according to transposon binding. Reads aligning to the + strand are offset by +4 bps, and reads aligning to the – strand are offset −5 bps (`sample_shifted_reads.sam`). 

The script has to make a `sam` header (`sample.sort.nodup.header.sam`) in order to convert the `sam` to `bam` and finally sort (`sample_shifted_reads_sorted.bam`).

The `sam` file is removed to save space but intermediate `bam` files are kept.

A [tag directory](http://homer.ucsd.edu/homer/ngs/tagDir.html) is made with `homer` which will be used for peak calling. 

Finally, `bigwig` files are made using `bamCoverage` and the final shifted, sorted `bam` file with the following parameters: `--binSize 30 -e -of bigwig --scaleFactor 1.0 --normalizeUsing RPKM -p 16`
