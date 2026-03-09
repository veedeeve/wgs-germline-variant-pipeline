# Whole-Genome Germline Variant Calling - Methodology

## Dataset

Sample: HG00096  
Source: 1000 Genomes Project

Files:
- SRR062634_1.filt.fastq.gz
- SRR062634_2.filt.fastq.gz

Sequencing type:
Paired-end reads, 100 bp.

## Variant Calling
### 1) Raw Reads
Reads for sample HG00096 were downloaded from the 1000 Genomes Project with `wget`
Raw reads were not trimmed using Trimmomatic because quality scores passed initial inspection.

### 2) Alignment
Align reads to hg38 reference genome using **BWA-MEM** due to paired-end reads longer than ~70bp
Read group metadata was included for downstream GATK4 analysis
```bash 
bwa mem -t 4 -R "@RG\tID:SRR062634\tPL:ILLUMINA\tSM:SRR062634" hg38.fa SRR062634_1.filt.fastq SRR062634_2.filt.fastq > SRR062634.paired.sam
```

### 3) Convert to BAM + Mark Duplicates
Using **GATK** to sort the coordinates, identify and mark duplicates to prevent duplicate reads affecting variant detection.
```bash 
gatk MarkDuplicatesSpark -I SRR06234.paired.sam -O SRR062634_sorted_dedup_reads.bam
```
### 4) Recalibrate Base Quality
Using **GATK** to recalibrate base quality score and apply the model
```bash 
gatk BaseRecalibrator -I SRR062634_sorted_dedup_reads.bam -R hg38.fa --known-sites Homo_sapiens_assembly38.dbsnp138.vcf.gz -O recal_data.table
gatk ApplyBQSR -I SRR062634_sorted_dedup_reads.bam -R hg38.fa --bqsr-recal-file recal_data.table -O SRR062634_sorted_dedup_bqsr_reads.bam
```
### 5) Collect Alignment and Insert Size Metrics
Obtain alignment statistics and insert size distributions using **GATK** tools.
```bash
gatk CollectAlignmentSummaryMetrics R=hg38.fa I=SRR062634_sorted_dedup_bqsr_reads.bam O=alignment_metrics.txt
gatk CollectInsertSizeMetrics INPUT=SRR062634_sorted_dedup_bqsr_reads.bam OUTPUT=insert_size_metrics.txt HISTOGRAM_FILE=insert_size_histogram.pdf
```
MultiQC was used to summarize quality metrics (located in results folder)

### 6) Call Variants
Variants were identified using **GATK HaplotypeCaller**, reconstructing haplotypes and detects SNPs and indels
```bash
gatk HaplotypeCaller -R hg38.fa -I SRR062634_sorted_dedup_bqsr_reads.bam -O raw_variants.vcf
```
### 7) Extract SNPS and INDELS
SNPs and indels were separated from raw VCF file using **GATK SelectVariants**
```bash
gatk SelectVariants -R hg38.fa -V raw_variants.vcf --select-type SNP -O raw_snps.vcf
gatk SelectVariants -R hg38.fa -V raw_variants.vcf --select-type INDEL -O raw_indels.vcf
```
## Variant Filtering and Annotation
### 8) Variant Filtering
SNPs and indels were filtered using **GATK Variant Filtration**, removing low-confidence variants based on standard hard-filtering thresholds. Indels require stricter filtering due to difficulty to detect, while SNPs are easier to detect.

For **SNPs**:
- Quality by depth (QD) < 2.0
- Strand bias (FS) > 60.0
- Mapping quality (MQ) < 40.0
- Strand odd ratio (SOR) > 4.0
- Mapping quality difference (MQRankSum) < 12.5
- Read position bias (ReadPosRankSum) < -8.0

For **indels**:
- QD < 2.0
- FS > 200.0
- SOR > 10.0

In addition, genotype-level filters were applied to remove low-confidence genotype calls:
- DP < 10
- GQ < 10

Only variants that passed all site-level and genotype-level filters were retained using **GATK SelectVariants** and ```cat``` to generate analysis ready SNP and indel VCF files

### 9) Variant Annotation
High-confidence SNP and indel sets were annotated with hg38 data source bundle. Adding functional annotation to each variant (gene-level and transcript-level).
```bash
gatk Funcotator \
	--variant analysis-ready-snps-filteredGT.vcf \
	--reference hg38.fa \
	--ref-version hg38 \
	--data-sources-path funcotator_dataSources.v1.7.20200521g \
	--output analysis-ready-snps-filteredGT-funcocated.vcf \
	--output-file-format VCF
```

Annotated VCF files were converted into tab-delimited tables using **GATK VariantsToTable**.
Fields extracted:
- AC
- AN
- DP
- AF
- Funcotation

```bash
gatk VariantsToTable \
	-V analysis-ready-snps-filteredGT-funcocated.vcf \
	-F AC \
	-F AN \
	-F DP \
	-F AF \
	-F FUNCOTATION \
	-O output_snps.table
```
















