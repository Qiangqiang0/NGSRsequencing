# NGS Resequencing

# GWAS, Construct genetic map, population evolution


## 1. GATK TOOLS PIPLINE

* __germline__: preprocess --> haplotyperCaller --> VQSR --> downstream

* __somatic__: preprocess --> mutect

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

#### 5. Preprocess

duplicates come fom PCR amplification

BQSR
```bash

picard MarkDuplicates I= $inputs  O= $dup M= $metrics

samtools index $dup

#2. BQSR

bqsr=${dup%.*}.bqsr
BQSRReport=${bqsr}.table

#bqsr report
gatk4 BaseRecalibrator -R $ref -I $dup -O $BQSRReport $knownSites

gatk4 GatherBQSRReports -I $BQSRReport -O $BQSRReport.report


#Apply bqsr

gatk4 ApplyBQSR -R $ref -I $dup -O $bqsr.bam -bqsr $BQSRReport \
      --static-quantized-quals 10 --static-quantized-quals 20 --static-quantized-quals 30 \
      --add-output-sam-program-record \
      --create-output-bam-md5 \
      --use-original-qualities


```

#### 6. Germline indel

ref: https://gatk.broadinstitute.org/hc/en-us/articles/360035535932-Germline-short-variant-discovery-SNPs-Indels-

ref: https://github.com/gatk-workflows

* variant calling:HaplotypeCaller: slow but efficient,UnifiedGenotyper: less efficient but faster

* using __--dbsnp__ to have snp reference

```bash
#1. HaplotyperCaller

gatk4 HaplotypeCaller \
        -R $ref \
        -I $ACC \
        -O $outdir\$ACC.g.vcf \
        -ERC GVCF \
        -contamination 0

#2. consolidate with genomicsDBImport

#Check dbimport directory

gatk4 GenomicsDBImport \
        --genomicsdb-workspace-path dbimport \
        -L $intervals \
        --sample-name-map vcf.txt \
        --reader-threads 5 \
        --merge-input-intervals \
        --consolidate


#3. GenotypeGVCFs: joint-call corhort
gatk4 GenotypeGVCFs \
        -R $ref \
        -O out.db.vcf \
        -V gendb://dbimport \
        -L $intervals \
        --merge-input-intervals

```

#### 7. VQSR

If having joint data: VQSR

If having single g.vcf: CNNScoreVariant

You should also know Hard filter 

ref: https://gatk.broadinstitute.org/hc/en-us/articles/360035531612-Variant-Quality-Score-Recalibration-VQSR-
ref: https://gatk.broadinstitute.org/hc/en-us/articles/360035531112--How-to-Filter-variants-either-with-VQSR-or-by-hard-filtering

```bash

## 4.1 VariantFiltration
gatk4 VariantFiltration \
        -V out.db.vcf \
        --filter-expression "ExcessHet > 54.69" \
        --filter-name ExcessHet \
        -O out.db.excess.vcf.gz

## 4.2 sites only
gatk4 MakeSitesOnlyVcf \
        -I out.db.excess.vcf.gz \
        -O out.db.excess.sitesonly.vcf.gz

## 4.3 Calculate VQSLOD
#variantRecalibrator for indel
gatk4 VariantRecalibrator \
        -V out.db.excess.sitesonly.vcf.gz \
        ---trust-all-polymorphic \
            -tranche 100.0 -tranche 99.95 -tranche 99.9 -tranche 99.5 -tranche 99.0 -tranche 97.0 -tranche 96.0 -tranche 95.0 -tranche 94.0 -tranche 93.5 -tranche 93.0 -tranche 92.0 -tranche 91.0 -tranche 90.0 \
    -an FS -an ReadPosRankSum -an MQRankSum -an QD -an SOR -an DP \
    -mode INDEL \
    --max-gaussians 4 \
    -resource $resource \
    -O out.db.excess.sitesonly.indel.recal \
    --tranches-file out.db.excess.sitesonly.indel.tranches

#variantRecalibrator for snp
gatk4 VariantRecalibrator \
        -V out.db.excess.sitesonly.vcf.gz \
        ---trust-all-polymorphic \
            -tranche 100.0 -tranche 99.95 -tranche 99.9 -tranche 99.5 -tranche 99.0 -tranche 97.0 -tranche 96.0 -tranche 95.0 -tranche 94.0 -tranche 93.5 -tranche 93.0 -tranche 92.0 -tranche 91.0 -tranche 90.0 \
    -an FS -an ReadPosRankSum -an MQRankSum -an QD -an SOR -an DP \
    -mode SNP \
    --max-gaussians 4 \
    -resource $resource \
    -O out.db.excess.sitesonly.snp.recal \
    --tranches-file out.db.excess.sitesonly.snp.tranches

## 4.4 applyVQSR
#ApplyVQSR for indel
gatk4 ApplyVQSR \
        -V out.excess.vcf.gz \
        --recal-file out.db.excess.sitesonly.indel.recal \
        --tranches-file out.db.excess.sitesonly.indel.tranches \
        --turth-sentitivity-filter-level 99.7 \
        --creat-output-variant-index true \
        -mode INDEL \
        -O out.db.excess.indel.recal.vcf.gz



#ApplyVQSR for snp

gatk4 ApplyVQSR \
        -V out.excess.vcf.gz \
        --recal-file out.db.excess.sitesonly.snp.recal \
        --tranches-file out.db.excess.sitesonly.snp.tranches \
        --turth-sentitivity-filter-level 99.7 \
        --creat-output-variant-index true \
        -mode INDEL \
        -O out.db.excess.snp.recal.vcf.gz

gatk CollectVariantCallingMetrics \
    -I out.db.excess.indel.recal.vcf.gz \
    --DBSNP $dbsnp \
    -SD $dbsnp.dict \
    -O indel.metrics

gatk CollectVariantCallingMetrics \
    -I out.db.excess.snp.recal.vcf.gz \
    --DBSNP $dbsnp \
    -SD $dbsnp.dict \
    -O snp.metrics

```

HERE SHOWS HARD FILTER
```bash
#select variant
gatk SelectVariants \
    -V cohort.vcf.gz \
    -select-type SNP \
    -O snps.vcf.gz



gatk SelectVariants \
    -V cohort.vcf.gz \
    -select-type INDEL \
    -O indels.vcf.gz


# variant filtration
gatk VariantFiltration \
    -V snps.vcf.gz \
    -filter "QD < 2.0" --filter-name "QD2" \
    -filter "QUAL < 30.0" --filter-name "QUAL30" \
    -filter "SOR > 3.0" --filter-name "SOR3" \
    -filter "FS > 60.0" --filter-name "FS60" \
    -filter "MQ < 40.0" --filter-name "MQ40" \
    -filter "MQRankSum < -12.5" --filter-name "MQRankSum-12.5" \
    -filter "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum-8" \
    -O snps_filtered.vcf.gz

gatk VariantFiltration \
    -V indels.vcf.gz \
    -filter "QD < 2.0" --filter-name "QD2" \
    -filter "QUAL < 30.0" --filter-name "QUAL30" \
    -filter "FS > 200.0" --filter-name "FS200" \
    -filter "ReadPosRankSum < -20.0" --filter-name "ReadPosRankSum-20" \
    -O indels_filtered.vcf.gz



# filter evaluation
gatk CollectVariantCallingMetrics \
    -I filtered.vcf.gz \
    --DBSNP Homo_sapiens_assembly38.dbsnp138.vcf \
    -SD Homo_sapiens_assembly38.dict \
    -O metrics
~

```

### somantic detection

ref: https://gatk.broadinstitute.org/hc/en-us/articles/360035890491?id=11127
ref: https://gatk.broadinstitute.org/hc/en-us/articles/360035890631-Panel-of-Normals-PON-

Panel of Normals(PON):A Panel of Normal or PON is a type of resource used in somatic variant analysis. Depending on the type of variant you're looking for, the PON will be generated differently. What all PONs have in common is that (1) they are made from normal samples (in this context, "normal" means derived from healthy tissue that is believed to not have any somatic alterations) and (2) their main purpose is to capture recurrent technical artifacts in order to improve the results of the variant calling analysis.

There is no definitive rule for how many samples should be used to make a PON (even a small PON is better than no PON) but in practice we recommend aiming for a minimum of 40.


```bash
# 1. create pon

gatk4 Mutect2  -R ${ref} -I \
	$input \
	-max-mnp-distance 0 --independent-mates  \
	-L $intervals -O $input.vcf.gz

gatk4 GenomicsDBImport --genomicsdb-workspace-path $pon_db -R $ref -V $input -L $intervals

gatk4 SelectVariants -R ${ref} -V gendb://$pon_db -o check.vcf

gatk4 CreateSomaticsPanelOfNormals -R $ref -V gendb://$pon_db -O pon_out.vcf 


# OR you could try this ref: https://gatk.broadinstitute.org/hc/en-us/articles/360036900132-CreateSomaticPanelOfNormals-BETA-

 gatk Mutect2 \
   -R ref_fasta.fa \
   -I normal1.bam \
   -tumor normal1_sample_name \
   --germline-resource af-only-gnomad.vcf.gz \
   -O normal1_for_pon.vcf.gz

 gatk CreateSomaticPanelOfNormals \
   -vcfs normal1_for_pon_vcf.gz \
   -vcfs normal2_for_pon_vcf.gz \
   -vcfs normal3_for_pon_vcf.gz \
   -O pon.vcf.gz
 
# 2.  mutect
#ref: https://gatk.broadinstitute.org/hc/en-us/articles/360035894731-Somatic-short-variant-discovery-SNVs-Indels-

tumor_command_line="-I tumor.bam -tumor tumor_name"
normal_command_line="-I normal.bam -normal normal_name"

gatk4 Mutect2 \
    -R {ref} \
    $tumor_command_line \
    $normal_command_line \
    --germline-resource $germline-resource \
    -pon pon_out.vcf \
    --f1r2-tar-gz flr2.tar.gz
    -L intervals \
    -O somatic.vcf.gz

 # if single tumor sample then without normal_command_line


# 3. contamination assessment 

gatk4 GetPileupSummaries \
     -R {ref} \
	 -I $tumorBam \
     -V $small_exac_common_3.hg38.vcf.gz \
     -O tumor-pileups.table

# normal tissue is optional
gatk4 GetPileupSummaries \
     -R {ref} \
     -I $NormalBam \
     -V $small_exac_common_3.hg38.vcf.gz \
     -O Normal-pileups.table

# read orientation model
gatk LearnReadOrientationModel -I f1r2.tar.gz -O read-orientation-model.tar.gz

# cauculate contamition
gatk4 CalculateContamination \
    -I tumor-pileups.table \
    -O tumor-cauculation.table \
    -matched Normal-pileups.table

# Filter mutect calls
gatk4 FilterMutectCalls \
    -V somatic.vcf.gz \
    --contamination-table tumor-cauculation.table \
    --ob-priors read-orientation-model.tar.gz \
    -O somatic_match_m2_oncefiltered.vcf.gz

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


### 2. Stacks


## GBS

ref for GBS

ref for difference between GBS and RAD-seq: http://www.360doc.com/content/18/0220/11/33459258_730961471.shtml


### 1. TASSEL-GBS




## BSA

ref: https://m.imooc.com/article/details?article_id=268903


## vcf data analysis

### 1. annovar based documentation

1. what is annovar?

ref: https://doc-openbio.readthedocs.io/projects/annovar/en/latest/

	* Gene-based annotation: identify whether SNPs or CNVs cause protein coding changes and the amino acids that are affected

	* Region-based annotation: identify variants in specific genomic regions, for example, conserved regions among 44 species, predicted transcription factor binding sites, segmental duplication regions...

	* Filter-based annotation: identify variants that are documented in specific databases

	* Other functionalities: Retrieve the nucleotide sequence in any user-specific genomic positions in batch, identify a candidate gene list for Mendelian diseases from exome data, and other utilities.



### 2. snpEff

download snpEff from: 

sneff build:

```bash
cd snpEff
mkdir genomes
mkdir arath

```

snpEff 


### gff3 to gtf tools

gffread in cufflinks


gffread gff3 -T -o gtf

gffread gtf -o- > gff3
1. ncbi - genome
2. ensemble - http://asia.ensembl.org/info/data/ftp/index.html

	* http://plants.ensembl.org/index.html