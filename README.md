# Bulk ATAC-seq Analysis Pipeline
Original scripts graciously provided by Dr. Murad and PhD candidate Klebea Carvalho.

In 6 bash scripts, this pipeline does the following:
1) Check read quality with `fastqc`
2) Map reads with `bowtie2`, remove duplicates with `Picard`, and shift reads due to transposase binding
3) Call 150bp and 500bp peaks with `Homer`
4) Run [IDR](https://github.com/karmel/homer-idr) (Irreproducibility Discovery Rate) on replicates
5) Merge peaks
6) Make counts matrix with `Homer`


## Check read quality
First, make a prefixes text file containing the names of samples you want to process:
```
ls AT_AC*_R1.fastq.gz | sed 's/_R1.fastq.gz//' > ../prefixes.txt
```

To run `fastqc` v. 0.11.9:
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
First, make sure there is a `bowtie2` indexed M chromosome and reference genome. If not, make one from a reference fasta:
```
cd ref/chrM/
module load bowtie2/2.4.1
bowtie2-build Homo_sapiens.GRCh38.dna.chromosome.MT.fa chrM
```

```
cd ref/hg38/
module load bowtie2/2.4.1
bowtie2-build GRCh38_encode.fasta hg38
```

Next, make sure `Homer` is installed. If not, install v. 4.11 from [Anaconda](https://anaconda.org/bioconda/homer). 
```
conda install -c bioconda homer
```

To run:
```
sbatch map_step2.sh
```

The outputs are in a subdirectory called `mapped`:
```
AT_AC_5_S3/
    fastqc/
    mapped/
        stuff
AT_AC_6_S4/
    fastqc/
    mapped/
        ...
```

## Call peaks with Homer
Call both 150bp and 500bp peaks.

To run:
```
sbatch peaks_step3.sh
```
