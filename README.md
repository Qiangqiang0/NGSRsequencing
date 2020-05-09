# Bioinformatics


## 1. genome resequencing


### 1. Basics

preparations: genome

out: genome index

raw_data(fastq) --> trimmomatic(fastq) --> bwa(sam) -->samtools(bam,sorted) --> GATK(vcf)

### 2. code

0. fastdump

```bash
# sra fastq data, `fastq-dump` and `fasterq-dump`
# --split-3 PE to two files, single end to one file
fasterq-dump --split-3 -O outdir seq.sra
```

0. reference index `bwa` or `samtools`

```bash
bwa index ref.fa
#samtools fadix ref.fa
```

1. fastqc

```bash

fastqc -o outdir -t threads fastq1 fastq2
```

2. trimmomantic and another fastqc on the output

strongly suggested to page: http://www.usadellab.org/cms/?page=trimmomatic for details

```bash
# single end with SE
java -jar trimmomatic-0.38.jar SE -threads 4 -phred33 fastq output_fastq ILLUMINACLIP:/Trimmomatic-0.38/adapters/TruSeq2-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36

## paird end with PE
java -jar trimmomatic-0.38.jar PE -threads 4 -phred33 forward.fastq inverse.fastq output_forward.fastq output_inverse.fastq ILLUMINACLIP::TruSeq3-PE.fa:2:30:10:2:keepBothReads LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
```

3. bwa

https://www.cnblogs.com/Formulate0303/p/7826944.html

http://davetang.org/wiki/tiki-index.php?page=SAM 

bwa mem works better for 70-100bp

baw -M: mark shorter split hits as secondary, otherwise it will be makred as supplementary

```bash
bwa mem -t 2 -R "@RG\tID:A\tSM:A" -M ref.fa  fastq > fastq_se.sam
bwa mem -t 2 -R "@RG\tID:A\tSM:A" -M ref.fa  fastq1 fastq2 > fastq_pe.sam
```

4. samtools: sort 
```bash
samtools view -bS fastq.sam > fastq.bam
samtools sort fastq.bam > fastq.sorted.bam
# or 
samtolls view -bS fastq.sam | samtools sort >fastq.sorted.bam
# or 
java -jar picard.jar SortSam 
```

