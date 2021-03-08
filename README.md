# Bulk ATAC-seq Analysis Pipeline
Original scripts graciously provided by Dr. Murad and PhD candidate Klebea Carvalho.

In 7 bash scripts, this pipeline does the following:
1) Check read quality with fastqc
2) (Optional) trim reads with cutadapt
3) Map reads with bowtie, remove duplicates with Picard, and shift reads due to transposase binding.
4) Call 150bp and 500bp peaks with Homer
5) Run IDR (Irreproducibility Discovery Rate) on replicates
6) Merge peaks
7) Make counts matrix with Homer

## Check read quality
```
qsub fastqc_step1.sh
```
