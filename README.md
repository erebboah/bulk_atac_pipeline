# Bulk ATAC-seq Analysis Pipeline
Original scripts graciously provided by Dr. Murad and PhD candidate Klebea Carvalho.

In 6 bash scripts, this pipeline does the following:
1) Check read quality with fastqc
2) Map reads with bowtie, remove duplicates with Picard, and shift reads due to transposase binding
3) Call 150bp and 500bp peaks with Homer
4) Run IDR (Irreproducibility Discovery Rate) on replicates
5) Merge peaks
6) Make counts matrix with Homer


## Check read quality
First, make a prefixes file containing the names of samples you want to process:
```
ls AT_AC*_R1.fastq.gz | sed 's/_R1.fastq.gz//' > ../prefixes.txt
```

To run fastqc v. 0.11.9:
```
sbatch fastqc_step1.sh
```

All pipeline outputs are in a directory with the same sample names as in the prefixes file:
```
AT_AC_5_S3/
    fastqc/
        AT_AC_5_S3_R1_fastqc.html
        AT_AC_5_S3_R1_fastqc.zip
        AT_AC_5_S3_R2_fastqc.html
        AT_AC_5_S3_R2_fastqc.zip
AT_AC_6_S4/
    fastqc/
     ...     
```

## Map, remove duplicates, and shift reads
To run:
```
sbatch map_step2.sh
```
