#!/bin/bash

# Directory
# /users/vydang/variant-calling

#########################################
# Filter & Annotate Variants
#########################################

ref="$HOME/variant-calling/data/hg38.fa"
results="$HOME/variant-calling/results"

# ---------- 1) Filter Variants ----------

# Filter SNPs
gatk VariantFiltration -R ${ref} \
	-V ${results}/raw_snps.vcf \
	-O ${results}/filtered_snps.vcf \
	-filter-name "QD_filter" -filter "QD < 2.0" \
	-filter-name "FS_filter" -filter "FS > 60.0" \
	-filter-name "MQ_filter" -filter "MQ < 40.0" \
	-filter-name "SOR_filter" -filter "SOR > 4.0" \
	-filter-name "MQRankSum_filter" -filter "MQRankSum < -12.5" \
	-filter-name "ReadPosRankSum_filter" -filter "ReadPosRankSum < -8.0" \
	-genotype-filter-expression "DP < 10" \
	-genotype-filter-name "DP_filter" \
	-genotype-filter-expression "GQ < 10" \
	-genotype-filter-name "GQ_filter" \

# Filter Indels
gatk VariantFiltration -R ${ref} \
	-V ${results}/raw_indels.vcf \
	-O ${results}/filtered_indels.vcf \
	-filter-name "QD_filter" -filter "QD < 2.0" \
	-filter-name "FS_filter" -filter "FS > 200.0" \
	-filter-name "SOR_filter" -filter "SOR > 10.0" \
	-genotype-filter-expression "DP < 10" \
	-genotype-filter-name "DP_filter" \
	-genotype-filter-expression "GQ < 10" \
	-genotype-filter-name "GQ_filter" \

# Select variants that PASS filters
gatk SelectVariants \
	--exclude-filtered \
	-V ${results}/filtered_snps.vcf \
	-O ${results}/analysis-ready-snps.vcf

gatk SelectVariants \
	--exclude-filtered \
	-V ${results}/filtered_indels.vcf \
	-O ${results}/analysis-ready-indels.vcf

cat ${results}/analysis-ready-snps.vcf | grep -v -E "DP_filter|GQ_filter" > ${results}/analysis-ready-snps-filteredGT.vcf
cat ${results}/analysis-ready-indels.vcf | grep -v -E "DP_filter|GQ_filter" > ${results}/analysis-ready-indels-filteredGT.vcf


# ---------- 2) Annotated Variants ----------
# Obtained GATK pre-packaged data (https://gatk.broadinstitute.org/hc/en-us/articles/30332018805787-Funcotator)
gatk_source="$HOME/variant-calling/funcotator_dataSources.v1.7.20200521g"

gatk Funcotator \
	--variant ${results}/analysis-ready-snps-filteredGT.vcf \
	--reference ${ref} \
	--ref-version hg38 \
	--data-sources-path ${gatk_source} \
	--output ${results}/analysis-ready-snps-filteredGT-funcocated.vcf \
	--output-file-format VCF

gatk Funcotator \
	--variant ${results}/analysis-ready-indels-filteredGT.vcf \
	--reference ${ref} \
	--ref-version hg38 \
	--data-sources-path ${gatk_source} \
	--output ${results}/analysis-ready-indels-filteredGT-funcocated.vcf \
	--output-file-format VCF


# Extract fields from VCF to tab-delimited table
gatk VariantsToTable \
	-V ${results}/analysis-ready-snps-filteredGT-funcocated.vcf \
	-F AC \
	-F AN \
	-F DP \
	-F AF \
	-F FUNCOTATION \
	-O ${results}/output_snps.table

gatk VariantsToTable \
	-V ${results}/analysis-ready-indels-filteredGT-funcocated.vcf \
	-F AC \
	-F AN \
	-F DP \
	-F AF \
	-F FUNCOTATION \
	-O ${results}/output_indels.table

# ---------- OPTIONAL) Obtain details for gene NBF1 ----------
cat ${results}/analysis-ready-snps-filteredGT-funcocated.vcf | grep "Funcotation fields are:" | sed 's/|/\t/g' > ${results}/output_curated_variants.txt

# Obtain lines for gene: NBF1
cat ${results}/output_snps.table | cut -f 5 | grep "NBPF1" | sed 's/|/\t/g' >> ${results}/output_curated_variants.txt



