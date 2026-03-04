# Human Germline Variant Calling & Functional Annotation Pipeline

## Project Overview
This project implements an end-to-end germline variant discovery workflow following GATK Best Practices. The pipeline processes human whole genome sequencing (WGS) data from raw FASTQ files through alignment, base quality score recalibration (BQSR), variant calling, hard filtering, and functional annotation.

## Objective
- Perform high-confidence germline SNP and INDEL discovery from human WGS data
- Apply site-level and genotype-level quality filtering
- Annotate variants with gene-level functional information
- Generate structured outputs for interpretation and evaluation
- Demonstrate implementation of GATK Best Practices in a reproducible pipeline

## Key Findings
- Successfully implemented a complete germline variant calling workflow using GATK4
- Generated high-confidence SNP and INDEL callsets from aligned WGS data
- Applied separate hard filtering criteria for SNPs and INDELs
- Performed genotype-level quality filtering based on depth (DP) and genotype quality (GQ)
- Annotated variants using Funcotator to extract gene-level information (e.g., NBPF1)
- Produced tab-delimited variant tables suitable for downstream biological analysis

## File Structure
```
├── scripts/
│   ├── variant-calling.sh
│   ├── variant-filtering-annotation.sh
├── results/
│   ├── raw_variants.vcf
│   ├── filtered_snps.vcf
│   ├── filtered_indels.vcf
│   ├── analysis-ready-snps.vcf
│   ├── analysis-ready-indels.vcf
│   └── output_snps.table
├── docs/
│   └── methodology.md
└── README.md
```

## Results
The pipeline produced:
- Raw variant calls (SNPs and INDELs) from HaplotypeCaller
- Site-level filtered variant sets using GATK hard-filter thresholds
- Genotype-level filtered callsets
- Functionally annotated VCF files using Funcotator
- Tab-delimited variant tables for targeted gene analysis  
Key evaluation metrics:
```
Mapped Reads: 99.66%
Total Variants: [1,067,473]
SNP Count: [946,560]
INDEL Count: [120,913]
```

## Discussions
This project demonstrates implementation of a complete germline variant discovery workflow following GATK Best Practices, from raw FASTQ files to annotated and filtered VCF outputs. Quality control indicated strong mapping performance with 99.66% mapped. Because the input FASTQ files were pre-filtered, the overall coverage rate is lower than standards.  
Future work could include benchmarking against high-confidence truth sets and re-running the pipeline on higher-coverage datasets to compare performance characteristics.
