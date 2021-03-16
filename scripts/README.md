### fastqc_step1.sh 
Requires fastqc/0.11.9 and paired reads with `_R1.fastq.gz` and `_R2.fastq.gz` file extensions.

Takes in a file of sample names (`prefixes_all.txt`) and launches an array of jobs, one per sample. Adjust number of jobs with `--array=1-N`.

Makes a sample folder in the main directory (`/share/crsp/lab/seyedam/share/newATAC/`) that will store results.

### map_step2.sh 
Requires bowtie2/2.4.1, samtools/1.10, picard-tools/1.87, java/1.8.0, python/3.8.0, and ucsc-tools/v393, and paired reads with `_R1.fastq.gz` and `_R2.fastq.gz` file extensions.

Two inputs, in order, are paths pointing to 1) a `bowtie2` reference and 2) a `bowtie2` chromosome M reference.

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

A custom script (`shift.reads.py`) is used to shift reads according to transposon binding. Reads aligning to the + strand are offset by +4 bps, and reads aligning to the – strand are offset −5 bps (`sample_shifted_reads.sam`). 

The script has to make a `sam` header (`sample.sort.nodup.header.sam`) in order to convert the `sam` to `bam` and finally sort (`sample_shifted_reads_sorted.bam`).

The `sam` file is removed to save space but intermediate `bam` files are kept.

A [tag directory](http://homer.ucsd.edu/homer/ngs/tagDir.html) is made with `homer` which will be used for peak calling. 

Finally, `bigwig` files are made using `bamCoverage` and the final shifted, sorted `bam` file with the following parameters: `--binSize 30 -e -of bigwig --scaleFactor 1.0 --normalizeUsing RPKM -p 16`

### peaks_step3.sh
Requires `homer` and the `homer-tags` folder created in the previous step.

Peaks with a p-value over local tag count of up to .1 and overall Poisson p-value of up to .1 pass filters. Peaks must also be greater than 50bp apart from one another and the tag density at peaks to be 4-fold greater than in the surrounding 50kb region. We call 150bp and 500bp peaks.
```
findPeaks homer-tags -LP .1 -poisson .1 -style factor -size 150 -minDist 50 -localSize 50000
findPeaks homer-tags -LP .1 -poisson .1 -style factor -size 500 -minDist 50 -localSize 50000
```

### idr_step4.sh
IDR is a process better explained [here](https://github.com/karmel/homer-idr) and is adapted from [ENCODE](https://projecteuclid.org/journals/annals-of-applied-statistics/volume-5/issue-3/Measuring-reproducibility-of-high-throughput/10.1214/11-AOAS466.full)/[Dr. Kundaje](https://sites.google.com/site/anshulkundaje/projects/idr/deprecated) but could stand to be updated. For now I think an old reproducibility pipeline is better than none at all...

Peaks are called with recommended settings, for `-size 150` and `size -500`, for combined peaks, individual pseudoreps, and pooled pseudoreps.
```
 -LP .1 -poisson .1 -style factor -size 150 -minDist 50 -localSize 50000
 -LP .1 -poisson .1 -style factor -size 500 -minDist 50 -localSize 50000
 ```
 
 IDR analysis is run with `-threshold 0.01`
 
 ### get_peaks_step5.sh
 This script simply converts the IDR output `combined.peaks-top-set.txt` for 150 and 500bp peaks for all samples into bed files using `awk`, `tail`, and `sort`. I have not figured out a way to do this and the subsequent steps of merging and filtering in one bash script.
