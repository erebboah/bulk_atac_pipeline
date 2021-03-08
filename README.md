# Bulk ATAC-seq Analysis Pipeline
The goal of week 10 in Eco Evo 283 was to use what we have learned to do a small project relevant to our bioinformatics interests. 

I am adapting a bulk ATAC-seq pipeline written for old HPC using outdated software to one that will work using HPC3 with newer packages.

Original scripts graciously provided by Dr. Rabi Murad and PhD candidate Klebea Carvalho (Mortazavi lab).

In 6 bash scripts, this pipeline does the following in summary:
1) Check read quality with `fastqc`
2) Map reads with `bowtie2`, remove duplicates with `Picard`, and shift reads due to transposase binding
3) Call 150bp and 500bp peaks with `Homer`
4) Run [IDR](https://github.com/karmel/homer-idr) (Irreproducibility Discovery Rate) on replicates (UPDATE THIS)
5) Merge peaks
6) Make counts matrix with `Homer`

Please see this readme for details on each step.

## Check read quality
First, make a prefixes text file containing the names of samples you want to process (all samples can be run together regardless of genome):
```
cd data/
ls *_R1.fastq.gz | sed 's/_R1.fastq.gz//' > ../prefixes_all.txt
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
First, make a genome-specific prefixes text file containing the names of samples you want to map:
```
cd data/
ls AT_AC_*_S*_R1.fastq.gz | sed 's/_R1.fastq.gz//' > ../prefixes_hg38.txt
```

Next, make sure there is a `bowtie2` indexed M chromosome and reference genome. If not, make one from a reference fasta:
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

Finally, make sure `Homer` is installed. If not, install v. 4.11 from [Anaconda](https://anaconda.org/bioconda/homer). 
```
conda install -c bioconda homer
```

To run:
```
sbatch map_hg38_step2.sh
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

There is another mapping script called `map_mm10_step2.sh` that aligns to the mm10 genome.

## Call peaks with Homer
Call both 150bp and 500bp peaks. Makes use of `prefixes_all.txt` files since I am processing both human and mouse samples but either the file or script can be edited.

To run:
```
sbatch peaks_step3.sh
```
