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
	- `ID:SRR062634` — unique identifier for the read group
	- `PL:ILLUMINA` — sequencing platform, used by GATK for platform-specific error modeling
  	- `SM:SRR062634` — sample name, used by GATK to associate reads with the correct sample

### 3) Convert to BAM + Mark Duplicates
Using **GATK** to sort the coordinates, identify and mark duplicates to prevent duplicate reads affecting variant detection.
```bash 
gatk MarkDuplicatesSpark -I SRR06234.paired.sam -O SRR062634_sorted_dedup_reads.bam
```
Marking (not removing) duplicates prevents PCR artifacts from inflating variant allele frequencies while preserving the original read data

### 4) Recalibrate Base Quality
Using **GATK** to recalibrate base quality score and apply the model
```bash 
gatk BaseRecalibrator -I SRR062634_sorted_dedup_reads.bam -R hg38.fa --known-sites Homo_sapiens_assembly38.dbsnp138.vcf.gz -O recal_data.table
gatk ApplyBQSR -I SRR062634_sorted_dedup_reads.bam -R hg38.fa --bqsr-recal-file recal_data.table -O SRR062634_sorted_dedup_bqsr_reads.bam
```
1. Build recalibration model
	- dbSNP138 as known sites file, sourced from GATK resource bundle
2. Apply model

### 5) Collect Alignment and Insert Size Metrics
Obtain alignment statistics and insert size distributions using **GATK** tools.
```bash
gatk CollectAlignmentSummaryMetrics R=hg38.fa I=SRR062634_sorted_dedup_bqsr_reads.bam O=alignment_metrics.txt
gatk CollectInsertSizeMetrics INPUT=SRR062634_sorted_dedup_bqsr_reads.bam OUTPUT=insert_size_metrics.txt HISTOGRAM_FILE=insert_size_histogram.pdf
```
- MultiQC was used to summarize quality metrics (located in results folder)
- Insert size is the distance between pair-end reads

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
SNPs and indels were filtered using **GATK Variant Filtration**, removing low-confidence variants based on standard hard-filtering thresholds. 
- Indels require stricter filtering due to difficulty to detect, while SNPs are easier to detect.

**SNP filters:**
 
| Filter | Threshold | Rationale |
| :----- | :-------: | :-------- |
| Quality by Depth (QD) | < 2.0 | Flags low-confidence variant calls relative to depth |
| Strand Bias (FS) | > 60.0 | Flags variants supported predominantly by one strand |
| Mapping Quality (MQ) | < 40.0 | Removes poorly mapped reads |
| Strand Odds Ratio (SOR) | > 4.0 | Secondary strand bias metric |
| MQ Rank Sum (MQRankSum) | < -12.5 | Flags mapping quality differences between ref/alt reads |
| Read Position Bias (ReadPosRankSum) | < -8.0 | Flags variants near ends of reads |
 
**Indel filters:**
 
| Filter | Threshold | Rationale |
| :----- | :-------: | :-------- |
| Quality by Depth (QD) | < 2.0 | Flags low-confidence calls relative to depth |
| Strand Bias (FS) | > 200.0 | Stricter threshold due to indel complexity |
| Strand Odds Ratio (SOR) | > 10.0 | Stricter threshold due to indel complexity |
 
**Genotype-level filters applied to both:**
 
| Filter | Threshold | Rationale |
| :----- | :-------: | :-------- |
| Depth (DP) | < 10 | Removes low-coverage genotype calls |
| Genotype Quality (GQ) | < 10 | Removes low-confidence genotype assignments |

Only variants that passed all site-level and genotype-level filters were retained using **GATK SelectVariants** and `cat` to generate analysis ready SNP and indel VCF files

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
**Fields extracted:**
 
| Field | Description |
| :---- | :---------- |
| AC | Allele count |
| AN | Total allele number |
| DP | Depth of coverage |
| AF | Allele frequency |
| FUNCOTATION | Functional annotation (gene, transcript, consequence) |

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

---
 
## Parameter Decisions
 
| Step | Key Parameter | Value | Rationale |
| :--- | :------------ | :---- | :-------- |
| BWA-MEM | -R (ID) | SRR062634 | Unique read group identifier required by GATK |
| BWA-MEM | -R (PL) | ILLUMINA | Platform tag used by GATK for error modeling |
| BWA-MEM | -R (SM) | SRR062634 | Sample name required for GATK sample association |
| BQSR | --known-sites | dbSNP138 | Most comprehensive curated SNP resource for hg38 |
| HaplotypeCaller | output | raw_variants.vcf | Unfiltered VCF retains all candidates for type-specific filtering |
| Funcotator | --output-file-format | VCF | Preserves all variant fields for downstream flexibility |
| SNP filter | QD | < 2.0 | Flags low-confidence calls relative to depth |
| SNP filter | FS | > 60.0 | Flags variants supported predominantly by one strand |
| SNP filter | MQ | < 40.0 | Removes poorly mapped reads |
| SNP filter | SOR | > 4.0 | Secondary strand bias metric |
| SNP filter | MQRankSum | < -12.5 | Flags mapping quality differences between ref/alt reads |
| SNP filter | ReadPosRankSum | < -8.0 | Flags variants near ends of reads |
| Indel filter | FS | > 200.0 | Stricter threshold due to indel complexity |
| Indel filter | SOR | > 10.0 | Stricter threshold due to indel complexity |
| Genotype filter | DP | < 10 | Removes low-coverage genotype calls |
| Genotype filter | GQ | < 10 | Removes low-confidence genotype assignments |
