# Bulk ATAC-seq Analysis Pipeline
The goal of week 10 in Eco Evo 283 was to use what we have learned to do a small project relevant to our bioinformatics interests. 

I am adapting a bulk ATAC-seq pipeline written for old HPC using to one that will work using HPC3 in a more polished pipeline format.

Original scripts graciously provided by Dr. Rabi Murad and PhD candidate Klebea Carvalho (Mortazavi lab).

In 5 bash scripts, this pipeline does the following in summary:
1) Check read quality with `fastqc`
2) Map reads with `bowtie2`, remove duplicates with `Picard`, and shift reads due to transposase binding
3) Call 150bp and 500bp peaks with `Homer`
4) Run [IDR](https://github.com/karmel/homer-idr) (Irreproducibility Discovery Rate) on replicates
5) Merge peaks, remove regions in no-pass list, and make counts matrix with `Homer`

Please see [this](https://github.com/erebboah/bulk_atac_pipeline/tree/main/scripts) readme for details on each step.

## Check read quality
First, make a prefixes text file containing the names of samples you want to process.

I am working with 4 samples: 2 replicates of mouse C2C12 myoblasts built by me and my undergrad Isaryhia Rodriguez. The libraries were built using the [Active Motif kit](https://www.activemotif.com/documents/2182.pdf), cat. 53150, and sequenced on a NextSeq 500 with 50 million reads per paired-end 43x43 library. 

```
data/
    C2C12_MB_ER1_R1.fastq.gz
    C2C12_MB_ER1_R2.fastq.gz
    C2C12_MB_ER2_R1.fastq.gz
    C2C12_MB_ER2_R2.fastq.gz
    C2C12_MB_IR1_R1.fastq.gz
    C2C12_MB_IR1_R2.fastq.gz
    C2C12_MB_IR2_R1.fastq.gz
    C2C12_MB_IR2_R2.fastq.gz
```

```
cd data/
ls *_R1.fastq.gz | sed 's/_R1.fastq.gz//' > ../prefixes.txt
```

To run `fastqc` v. 0.11.9:
```
sbatch fastqc_step1.sh
```

All pipeline outputs are in a directory with the same sample names as in the prefixes file:
```
C2C12_MB_ER1/
    fastqc/
        C2C12_MB_ER1_R1_fastqc.html
        C2C12_MB_ER1_R1_fastqc.zip
        C2C12_MB_ER1_R2_fastqc.html
        C2C12_MB_ER1_R2_fastqc.zip
C2C12_MB_ER2/
    fastqc/
     ...     
```

The `fastqc` outputs are also hosted on our [lab website](/var/www/html/erebboah/bulk_atac/C2C12_MB). 

## Map, remove duplicates, and shift reads
If a `prefixes.txt` file has not been made, make one now of the samples you want to process.

Make sure there is a `bowtie2` indexed M chromosome and reference genome. If not, make one from an [M chromosome fasta](http://hgdownload.soe.ucsc.edu/goldenPath/mm10/chromosomes/) and [reference fasta](https://www.encodeproject.org/data-standards/reference-sequences/).

For mouse:
```
cd ref/chrM_mm10/
module load bowtie2/2.4.1
bowtie2-build mm10_chrM.fa chrM
```

```
cd ref/mm10/
module load bowtie2/2.4.1
bowtie2-build mm10_no_alt_analysis_set_ENCODE.fasta mm10
```

Make sure `Homer` is installed. If not, install v. 4.11 from [Anaconda](https://anaconda.org/bioconda/homer). 
```
conda install -c bioconda homer
```

Also install genomes:
```
perl /data/homezvol2/$USER/miniconda3/share/homer/.//configureHomer.pl -install mm10
perl /data/homezvol2/$USER/miniconda3/share/homer/.//configureHomer.pl -install hg38
```

To run:
```
sbatch map_step2.sh ../ref/mm10/mm10 ../ref/chrM_mm10/chrM
```

The outputs are in a subdirectory called `mapped`:
```
C2C12_MB_ER1/
    fastqc/
    mapped/
        C2C12_MB_ER1.bam
        C2C12_MB_ER1.bigWig
        C2C12_MB_ER1_chrM.bam
        C2C12_MB_ER1.duplicates_metric.txt
        C2C12_MB_ER1_shifted_reads_sorted.bam
        C2C12_MB_ER1_shifted_reads_sorted.bam.bai
        C2C12_MB_ER1.sort.bam
        C2C12_MB_ER1.sort.nodup.bam
        C2C12_MB_ER1.sort.nodup.header.sam
        C2C12_MB_ER1_unaligned.fastq.1.gz
        AC2C12_MB_ER1_unaligned.fastq.2.gz
        homer-tags/
            chr10.tags.tsv
            chr11.tags.tsv
            ...
            tagAutocorrelation.txt
            tagCountDistribution.txt
            tagInfo.txt
            tagLengthDistribution.txt
C2C12_MB_ER2/
    fastqc/
    mapped/
        ...
```

The directory `homer-tags` will be used during peak calling.

## Call peaks with Homer
This script also makes use of the `prefixes.txt` file and requires the `homer-tags` directory. It calls both 150bp and 500bp peaks using `Homer`.

To run:
```
sbatch peaks_step3.sh
```

The outputs are in a subdirectory called `peaks`:
```
C2C12_MB_ER1/
    fastqc/
    mapped/
    peaks/
        C2C12_MB_ER1.150bp.peaks.txt
        C2C12_MB_ER1.500bp.peaks.txt
C2C12_MB_ER2/
    fastqc/
    mapped/
    peaks/
        C2C12_MB_ER2.150bp.peaks.txt
        C2C12_MB_ER2.500bp.peaks.txt
```

## Run IDR
