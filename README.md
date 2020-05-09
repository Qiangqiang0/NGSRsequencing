# Bioinformatics


## 1. genome resequencing


### 1. Basics

preparations: genome

out: genome index

raw_data(fastq) --> trimmomatic(fastq) --> botwie(sam) -->samtools(bam,sorted) --> GATK(vcf)

### 3. code

0. fastdump

```bash
# sra fastq data, `fastq-dump` and `fasterq-dump`
# --split-3 PE to two files, single end to one file
fasterq-dump --split-3 -O outdir seq.sra
```
1. fastqc

```bash

fastqc -o outdir -t threads fastq1 fastq2
```

2. trimmomantic

```bash
# single end with SE
java -jar trimmomatic-0.38.jar SE -threads 4 -phred33 fastq output_fastq ILLUMINACLIP:/Trimmomatic-0.38/adapters/TruSeq2-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36

## paird end with PE
java -jar trimmomatic-0.38.jar PE -threads 4 -phred33 forward.fastq inverse.fastq output_forward.fastq output_inverse.fastq ILLUMINACLIP:/Trimmomatic-0.38/adapters/TruSeq2-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36

```





