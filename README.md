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

0. reference index `bwa` and `samtools`

```bash
bwa index ref.fa
samtools faidx ref.fa
```

0. install gatk
```bash
# conda

conda install gatk
# run gatk-register to get the download page
gatk-register /dowloadedgatk.jar|tar.bz2

```
1. fastqc

```bash

fastqc -o outdir -t threads fastq1 fastq2
```

2. trimmomantic and another fastqc on the output

strongly suggested to page: http://www.usadellab.org/cms/?page=trimmomatic for details

```bash
# single end with SE
java -jar trimmomatic-0.38.jar SE -threads 4 (-phred33|-phred64) fastq output_fastq ILLUMINACLIP:/Trimmomatic-0.38/adapters/TruSeq2-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36

## paird end with PE
java -jar trimmomatic-0.38.jar PE -threads 4 (-phred33|-phred64) forward.fastq inverse.fastq output_forward.fastq output_inverse.fastq ILLUMINACLIP::TruSeq3-PE.fa:2:30:10:2:keepBothReads LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
```

3. bwa

https://www.cnblogs.com/Formulate0303/p/7826944.html

http://davetang.org/wiki/tiki-index.php?page=SAM 

bwa mem works better for 70-100bp

baw -M: mark shorter split hits as secondary, otherwise it will be makred as supplementary

```bash
bwa mem -t threads -R "@RG\tID:A\tSM:A" -M ref.fa  fastq > fastq_se.sam
bwa mem -t threads -R "@RG\tID:A\tSM:A" -M ref.fa  fastq1 fastq2 > fastq_pe.sam
```

4. samtools: sort 
```bash
samtools view -bS fastq.sam > fastq.bam
samtools sort fastq.bam > fastq.sorted.bam
# or 
samtools view -bS fastq.sam | samtools sort >fastq.sorted.bam
# or 
java -jar picard.jar SortSam SORT_ORDER="coordinate"
```

5. marking duplicates

duplicates come fom PCR amplification

```bash
java -jar pycard.jar MarkDuplicates I= fastq.sorted.bam  O= fastq.sorted.dup.bam M= fastq.sorted.dup.metrics
```

6. indel detection
* realign: get the output of intervals to show possible indels

* variant calling:HaplotypeCaller: slow but efficient,UnifiedGenotyper: less efficient but faster

  * __multi-sample calling__: slow but more trustable
  * __single-sample calling__: generate gvcf and combined into vcf

```bash
picard CreateSequenceDictionary REFERENCE= ref.fa OUTPUT= ref.dict
gatk3 -T RealignerTargetCreator -R  ref.fa  -I fastq.sorted.dup.bam -o possible_indel.intervels
gatk3 -T IndelRealigner -R ref.fa -I fastq.sorted.dup.bam -o fastq.sorted.dup.bam.re --targetIntervals possible_index.intervels

# single-sample calling
gatk3 -T HaplotypeCaller -R ref.fa -o fastq.vcf -I fastq.sorted.dup.bam.re --emitRefConfidence GVCF -nct threads


```




