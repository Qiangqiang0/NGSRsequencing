# NGS Resequencing

# GWAS, Construct genetic map, population evolution


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

#### 1. fastqc

```bash

fastqc -o outdir -t threads fastq1 fastq2
```

#### 2. trimmomantic and another fastqc on the output

strongly suggested to page: http://www.usadellab.org/cms/?page=trimmomatic for details

```bash
# single end with SE
java -jar trimmomatic-0.38.jar SE -threads 4 (-phred33|-phred64) fastq output_fastq ILLUMINACLIP:/Trimmomatic-0.38/adapters/TruSeq2-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36

## paird end with PE
java -jar trimmomatic-0.38.jar PE -threads 4 (-phred33|-phred64) forward.fastq inverse.fastq output_forward.fastq output_inverse.fastq ILLUMINACLIP::TruSeq3-PE.fa:2:30:10:2:keepBothReads LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
```

#### 3. bwa

https://www.cnblogs.com/Formulate0303/p/7826944.html

http://davetang.org/wiki/tiki-index.php?page=SAM 

bwa mem works better for 70-100bp

baw -M: mark shorter split hits as secondary, otherwise it will be makred as supplementary

```bash
bwa mem -t threads -R "@RG\tID:A\tSM:A" -M ref.fa  fastq > fastq_se.sam
bwa mem -t threads -R "@RG\tID:A\tSM:A" -M ref.fa  fastq1 fastq2 > fastq_pe.sam
```

#### 4. samtools: sort
```bash
samtools view -bS fastq.sam > fastq.bam
samtools sort fastq.bam > fastq.sorted.bam
# or 
samtools view -bS fastq.sam | samtools sort >fastq.sorted.bam
# or 
java -jar picard.jar SortSam SORT_ORDER="coordinate"


```

#### 5. marking duplicates

duplicates come fom PCR amplification

```bash
java -jar picard.jar MarkDuplicates I= fastq.sorted.bam  O= fastq.sorted.dup.bam M= fastq.sorted.dup.metrics
```

#### 6. indel detection
* realign: get the output of intervals to show possible indels

* variant calling:HaplotypeCaller: slow but efficient,UnifiedGenotyper: less efficient but faster

  * __multi-sample calling__: slow but more trustable
  * __single-sample calling__: generate gvcf and combined into vcf

  * using __--dbsnp__ to have snp reference

```bash
picard CreateSequenceDictionary REFERENCE= ref.fa OUTPUT= ref.dict
gatk3 -T RealignerTargetCreator -R  ref.fa  -I fastq.sorted.dup.bam -o possible_indel.intervals
gatk3 -T IndelRealigner -R ref.fa -I fastq.sorted.dup.bam -o fastq.sorted.dup.re.bam --targetIntervals possible_indel.intervals

# single-sample calling
gatk3 -T HaplotypeCaller -R ref.fa -o fastq.gvcf -I fastq.sorted.dup.re.bam --emitRefConfidence GVCF -nct 24  -variant_index_type LINEAR -variant_index_parameter 128000
gatk3 -T   CombineGVCFs -R  ref.fa -o fastq.vcf -V fastq1.gvcf -V ...

# multi-sample calling
gatk3 -T HaplotypeCaller -R ref.fa -o fastq.vcf -I fastq.sorted.dup.re.bam  -nct 24

# 
gatk3 -T GenotypeGVCFs -R ref.fa -V fastq.vcf -O common.vcf.gz

# select SNP
gatk3 -T  SelectVariants -R ref.fa -O SNP.vcf --variant common.vcf.gz --select-type-to-include SNP 

# select indel
gatk3 -T  SelectVariants -R ref.fa -O indel.vcf --variant common.vcf.gz --select-type-to-include INDEL
```

#### S1. BQSR: base quality score recalibration

rectify the quality score of base,for humans download the data, for other creatures, creat it by yourself

if you DO BSQR,plsease go to 
ref: https://www.bioinfo-scrounger.com/archives/642/
ref: https://www.plob.org/article/7009.html
```bash
gatk3 -T BaseRecalibrator -R ref.fa -knownSites dbsnp_137.hg19.vcf -knownSites Mills_and_1000G_gold_standard.indels.hg19.vcf -knownSites 1000G_phase1.indels.hg19.vcf -I fastq.bam -o fastq.bam.grp
gatk3 -T ApplyBQSR -R ref.fa -I fastq.bam --bqsr-recal-file fastq.bam.grp -O $sample.sorted.marked.BQSR.bam

```




## RAD-seq:

### 1. basics

using enzyme to digest genome, collect sequence around enzyme clevage sites.

process:

1. digestion;

2. add pcr primer barcode to the collected sequence

3. random breaking the sequence

4. add another pcr primer adaptor2

5. pcr on barcode and adaptor

6. seq

7. if (genome_is_avaiable){
	align to get SNP
} else{
	using different samples to discover SNP.
	avaible tools: __RADtools__, __Stacks__
}


## GBS

ref for GBS

ref for difference between GBS and RAD-seq: http://www.360doc.com/content/18/0220/11/33459258_730961471.shtml


## BSA

ref: https://m.imooc.com/article/details?article_id=268903