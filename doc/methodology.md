# Whole-Genome Germline Variant Calling - Methodology

## Dataset

Sample: HG00096  
Source: 1000 Genomes Project

Files:
- SRR062634_1.filt.fastq.gz
- SRR062634_2.filt.fastq.gz

Sequencing type:
Paired-end reads, 100 bp.

## 1) Raw Reads
Reads for sample HG00096 were downloaded from the 1000 Genomes Project with `wget`
Raw reads were not trimmed using Trimmomatic because quality scores passed initial inspection.

## 2) Alignment
Align reads to hg38 reference genome using **BWA-MEM** due to paired-end reads longer than ~70bp
Read group metadata was included for downstream GATK4 analysis
```bash 
bwa mem -t 4 -R "@RG\tID:SRR062634\tPL:ILLUMINA\tSM:SRR062634" hg38.fa SRR062634_1.filt.fastq SRR062634_2.filt.fastq > SRR062634.paired.sam
```

## 3) Convert to BAM + Mark Duplicates
Using **GATK** to sort the coordinates, identify and mark duplicates to prevent duplicate reads affecting variant detection.
```bash 
gatk MarkDuplicatesSpark -I SRR06234.paired.sam -O SRR062634_sorted_dedup_reads.bam
```
## 4) Recalibrate Base Quality
Using **GATK** to recalibrate base quality score and apply the model
```bash 
gatk BaseRecalibrator -I SRR062634_sorted_dedup_reads.bam -R hg38.fa --known-sites Homo_sapiens_assembly38.dbsnp138.vcf.gz -O recal_data.table
gatk ApplyBQSR -I SRR062634_sorted_dedup_reads.bam -R hg38.fa --bqsr-recal-file recal_data.table -O SRR062634_sorted_dedup_bqsr_reads.bam
```





















