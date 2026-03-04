#!/bin/bash

# Directory
# /users/vydang/variant-calling

#########################################
# Call Indels & SNPs
#########################################

# ---------- 1) Download Files ----------
# Read files
wget -P ~/variant-calling/reads \
ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/phase3/data/HG00096/sequence_read/SRR062634_1.filt.fastq.gz
wget -P ~/variant-calling/reads \
ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/phase3/data/HG00096/sequence_read/SRR062634_2.filt.fastq.gz

# Reference files
wget -P ~/variant-calling/data \
https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz
gunzip ~/variant-calling/data/hg38.fa.gz

# Index reference - .fai file
samtools faidx ~/variant-calling/data/hg38.fa

# reference dictionary
gatk CreateSequenceDictionary \
R=~/variant-calling/data/hg38.fa \
O=~/variant-calling/data/hg38.dict

# Known variant sites for BQSR from GATK
wget -P ~/variant-calling/data/ \
https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.dbsnp138.vcf.gz

# ---------- 2) Variant Calling ----------
ref="$HOME/variant-calling/data/hg38.fa"
known_sites="$HOME/variant-calling/data/Homo_sapiens_assembly38.dbsnp138.vcf.gz"
aligned_reads="$HOME/variant-calling/aligned-reads"
reads="$HOME/variant-calling/reads"
results="$HOME/variant-calling/results"
data="$HOME/variant-calling/data"

# --------------------
# 1) Quality Control
# --------------------
fastqc ${reads}/SRR062634_1.filt.fastq.gz -o ${reads}/
fastqc ${reads}/SRR062634_2.filt.fastq.gz -o ${reads}/
# No trimming required

# --------------------
# 2) Map to reference (BWA-MEM)
# --------------------
# BWA index reference
bwa index ${ref}

# BWA Alignment
bwa mem -t 4 -R "@RG\tID:SRR062634\tPL:ILLUMINA\tSM:SRR062634" ${ref} ${reads}/SRR062634_1.filt.fastq.gz ${reads}/SRR062634_2.filt.fastq.gz > ${aligned_reads}/SRR062634.paired.sam

# Flagstat file
samtools flagstat ${aligned_reads}/SRR062634.paired.sam > ${data}/SRR062634-alignment.paired.sam.txt

# --------------------
# 3) Mark duplicates and sort (GATK4)
# --------------------
gatk MarkDuplicatesSpark -I ${aligned_reads}/SRR062634.paired.sam -O ${aligned_reads}/SRR062634_sorted_dedup_reads.bam

samtools flagstat ${aligned_reads}/SRR062634_sorted_dedup_reads.bam > ${data}/SRR062634_sorted_dedup_reads.bam.txt

# --------------------
# 4) Base quality recalibration
# --------------------

# Build model
gatk BaseRecalibrator -I ${aligned_reads}/SRR062634_sorted_dedup_reads.bam -R ${ref} --known-sites ${known_sites} -O ${data}/recal_data.table

# Apply model to adjust the base quality scores
gatk ApplyBQSR -I ${aligned_reads}/SRR062634_sorted_dedup_reads.bam -R ${ref} --bqsr-recal-file ${data}/recal_data.table -O ${aligned_reads}/SRR062634_sorted_dedup_bqsr_reads.bam

# --------------------
# 5) Collect Alignment & Insert Size Metrics
# --------------------
gatk CollectAlignmentSummaryMetrics R=${ref} I=${aligned_reads}/SRR062634_sorted_dedup_bqsr_reads.bam O=${aligned_reads}/alignment_metrics.txt
gatk CollectInsertSizeMetrics INPUT=${aligned_reads}/SRR062634_sorted_dedup_bqsr_reads.bam OUTPUT=${aligned_reads}/insert_size_metrics.txt HISTOGRAM_FILE=${aligned_reads}/insert_size_histogram.pdf #provide distribution of insert sizes across reads

# Quality Control
multiqc ${aligned_reads}/ --outdir ${data}/

# --------------------
# 6) Call Variants - gatk haplotype caller
# --------------------
gatk HaplotypeCaller -R ${ref} -I ${aligned_reads}/SRR062634_sorted_dedup_bqsr_reads.bam -O ${results}/raw_variants.vcf

# Extract SNPs & INDELS
gatk SelectVariants -R ${ref} -V ${results}/raw_variants.vcf --select-type SNP -O ${results}/raw_snps.vcf
gatk SelectVariants -R ${ref} -V ${results}/raw_variants.vcf --select-type INDEL -O ${results}/raw_indels.vcf



