# Bulk ATAC-seq Analysis Pipeline
The goal of week 10 in Eco Evo 283 was to use what we have learned to do a small project relevant to our bioinformatics interests. 

I am adapting a bulk ATAC-seq pipeline written for old HPC using to one that will work using HPC3 in a more polished pipeline format.

Original scripts graciously provided by Dr. Rabi Murad and PhD candidate Klebea Carvalho (Mortazavi lab, UCI).

In 5 bash scripts, this pipeline does the following:
1) Check read quality with `fastqc`
2) Map reads with `bowtie2`, remove duplicates with `Picard`, and shift reads due to transposase binding
3) Call 150bp and 500bp peaks with `Homer`
4) Run [IDR](https://github.com/karmel/homer-idr) (Irreproducibility Discovery Rate) on replicates
5) Convert IDR peaks to bed files

The user can then generate the final merged peak set, and use `Homer`'s `annotatePeaks.pl` function to make a counts matrix with the final merged peak set and `homer-tags` directories made in the mapping step. 

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

The `fastqc` outputs are also hosted on our [lab website](http://crick.bio.uci.edu/erebboah/bulk_atac/C2C12_MB/). 

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
        C2C12_MB_ER1_unaligned.fastq.2.gz
        C2C12_MB_ER1_homer-tags/
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
    ...
```

## Run IDR
This script makes use of the `.150bp.peaks.txt` and `.500bp.peaks.txt` files from the previous step.

To run per sample, with two replicates each:
```
sbatch idr_step4.sh C2C12_MB_ER1 C2C12_MB_ER2 C2C12_MB_ER
```

```
sbatch idr_step4.sh C2C12_MB_IR1 C2C12_MB_IR2 C2C12_MB_IR
```

The extensive outputs from `IDR` are in a new directory with the sample name. The directory structure is the same between `peaks150` and `peaks500`.

Tag directories (`combined`) for 150 and 500bp peaks are made using both replicates (e.g. C2C12_MB_ER1 and C2C12_MB_ER2) and peaks are called (`combined.peaks.txt`).

Using IDR python code, pseudorep tag directories are created for individual and pooled replicates (`pseudoreps/individual` and `pseudoreps/pooled`) and peaks are called (e.g. `pseudoreps/individual/C2C12_MB_ER1_homer-tags-Pseudorep1_peaks.txt` and `pseudoreps/pooled/combined-Pseudorep1_peaks.txt`).

Finally, IDR analysis is run to produce a final set of peaks passing a cutoff of 0.01 (`combined.peaks-top-set.txt`). These peaks are confident and more likely real.
```
C2C12_MB_ER/
    peaks150/
        C2C12_MB_ER1.150bp.peaks.txt
        C2C12_MB_ER2.150bp.peaks.txt
        combined/
            ...
            combined.peaks.txt
        idr-output/
            combined.peaks-top-set.txt
            narrowpeaks/
            plots/
            pooled_comparisons/
            pseudorep_comparisons/
            replicate_comparisons/
        pseudoreps/
            individual/
                C2C12_MB_ER1_homer-tags-Pseudorep1/
                C2C12_MB_ER1_homer-tags-Pseudorep1_peaks.txt
                C2C12_MB_ER1_homer-tags-Pseudorep2/
                C2C12_MB_ER1_homer-tags-Pseudorep2_peaks.txt
                C2C12_MB_ER2_homer-tags-Pseudorep1/
                C2C12_MB_ER2_homer-tags-Pseudorep1_peaks.txt
                C2C12_MB_ER2_homer-tags-Pseudorep2/
                C2C12_MB_ER2_homer-tags-Pseudorep2_peaks.txt
            pooled/
                combined-Pseudorep1/
                combined-Pseudorep1_peaks.txt
                combined-Pseudorep2/
                combined-Pseudorep2_peaks.txt
        replicates/
            C2C12_MB_ER1.150bp.peaks.txt
            C2C12_MB_ER2.150bp.peaks.txt
    peaks500/
        ...
        idr-output/
            combined.peaks-top-set.txt
C2C12_MB_IR/
    peaks150/
        ...
        idr-output/
            combined.peaks-top-set.txt
    peaks500/ 
        ...
        idr-output/
            combined.peaks-top-set.txt
```

## Generate final peak set
The `combined.peaks-top-set.txt` from each of the 150 and 500bp IDR analyses for all the samples in the experiment is the input into merging and filtering peaks.

The script also uses a `samples.txt` file with all of the samples you want to consider as part of the experiment.

To make bed files for all of the peaks:
```
sbatch get_peaks_step5.sh C2C12_MB
```

The output is in a new folder with the experiment name passed as the argument to the script:
```
C2C12_MB/
    C2C12_MB_ER_150bp.bed
    C2C12_MB_ER_500bp.bed
    C2C12_MB_IR_150bp.bed
    C2C12_MB_IR_500bp.bed
```

To merge peaks and filter out a no-pass list/exclusion list of repeat regions which can be found [here](https://github.com/Boyle-Lab/Blacklist/tree/master/lists), you can run the following:

```
cat C2C12_MB/*.bed > C2C12_MB/merged.peaks.bed
bedtools intersect -v -a C2C12_MB/merged.peaks.bed -b ref/mm10-nopasslist.v2.bed > C2C12_MB/merged.peaks.filt.bed
```

And finally a python script appends a Peak ID to each region from Peak1 to PeakN:
```
python scripts/addPeakName.py C2C12_MB/merged.peaks.filt.bed C2C12_MB/merged.peaks.filt.final.bed
```

The outputs are in the experiment folder:
```
C2C12_MB/
    C2C12_MB_ER_150bp.bed
    C2C12_MB_ER_500bp.bed
    C2C12_MB_IR_150bp.bed
    C2C12_MB_IR_500bp.bed
    merged.peaks.bed
    merged.peaks.filt.bed
    merged.peaks.filt.final.bed
```

## Generate counts matrix
Since the replicates for each experiment will change, I find it's easier to run the last step manually as well and use tab to complete to add in all the `homer-tags` directories that you used to make the merged peak set.

Make sure `Homer` is available and you have your genome of interest installed.

I like to output the final matrix in the same folder as the final peak set.

```
annotatePeaks.pl C2C12_MB/merged.peaks.filt.final.bed mm10 -raw -annStats C2C12_MB/annotationStats.txt -d C2C12_MB_ER1/mapped/C2C12_MB_ER1_homer-tags  C2C12_MB_ER2/mapped/C2C12_MB_ER2_homer-tags C2C12_MB_IR1/mapped/C2C12_MB_IR1_homer-tags C2C12_MB_IR2/mapped/C2C12_MB_IR2_homer-tags > C2C12_MB/matrix.raw.txt
```

The header is unwieldly because it includes the path of all the tag directories; the columns (1-19) contain: 
1. peakID
2. chromosome
3. start
4. stop
5. strand
6. peak score
7. focus ratio
8. annotation
9. detailed annotation
10. distance to TSS
11. nearest promoter
12. promoter ID
13. nearest Unigene ID
14. nearest Refseq ID
15. nearest Ensemble ID
16. gene name
17. gene alias
18. gene description
19. gene type,

followed by columns of counts for the input tag directories in the order they were listed in columns 20 onward.

I remove the header and keep columns of interest and the counts columns. Since I have 4 samples, I keep columns 20-23.
```
tail -n +2 C2C12_MB/matrix.raw.tsv | cut -f 1,2,3,4,15,16,20,21,22,23 > C2C12_MB/matrix.tsv
```
```
C2C12_MB/
    annotationStats.txt
    C2C12_MB_ER_150bp.bed
    C2C12_MB_ER_500bp.bed
    C2C12_MB_IR_150bp.bed
    C2C12_MB_IR_500bp.bed
    matrix.raw.tsv
    matrix.tsv
    merged.peaks.bed
    merged.peaks.filt.bed
    merged.peaks.filt.final.bed
    
```

# Results
## Genome browser
I uploaded the `bigWig` files from the mapping step to the genome browser using our lab server to [host](http://crick.bio.uci.edu/erebboah/bulk_atac/C2C12_MB/) them.

The track can be loaded by clicking [here](https://genome.ucsc.edu/s/erebboah/C2C12_MB_bulkATAC).

This screenshot shows the 4 C2C12 myoblast bulk ATAC samples at the Myod1 locus, a myogenic transcription factor we expect to be expressed in these cells. One replicate `C2C12_MB_IR1` looks like it is slightly less open at the promoter than the other 3 reps.

![myod1](https://github.com/erebboah/bulk_atac_pipeline/blob/main/figures/myod1_atac.png?raw=true)

## Top annotated genes
("Top" by C2C12_ER2, because it looks the best at first glance at the browser).

```
sort -rnk8  matrix.tsv | head -n 25 | awk -F "\t" '{print $6}' | sort | uniq
```
1. Adgrg3
2. Asf1b
3. Dnpep
4. Mir181c
5. Mir23a
6. Mob3a
7. Mtmr4
8. Myh9
9. Rfx1
10. Sp1

Since these cells are proliferative precursors to differentiated skeletal muscle, I don't expect a lot of myogenic genes to be the most accessible, but rather cell cycle genes.

Asf1b is a histone chaperone. MicroRNA Mir23a has evidence in [literature](https://www.sciencedirect.com/science/article/pii/S0014482712003060) as an inhibitor of myogenic differentiation in C2C12 (indicating that accessibility may go down in differentiated cells). 

MicroRNAs are known to be critical in skeletal muscle development and are known as ["myomiRs"](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4871711/). However, miRNAs such as miR-675 within a lncRNA such as H19 will also show up in chromatin data, maybe an artifact.

Interestingly, we also capture high accessibility of Myh9, which is a large ["non-muscle"](https://en.wikipedia.org/wiki/MYH9) myosin heavy chain gene, but plays a role in early vertebrate development.

Importantly, C2C12 is a flexible cell line that can be coaxed into differentiating to myotubes, adipocytes, and ostocytes.

## GREAT annotation
To see what genes are highly accessible, I uploaded the final merged peak set to [GREAT](http://great.stanford.edu/public/html/), which is another tool to annotate peaks, as well as does a ontology analysis of the genes....

## Future directions
This pipeline could still be more user-friendly, especially the last few steps.

As we add more samples, especially of differentiated cells, we can perform differential accessibility testing between samples using our counts matrix. 

With enough ATAC-seq samples, we can perform interesting analyses such as transcription factor footprinting and motif enrichment.
