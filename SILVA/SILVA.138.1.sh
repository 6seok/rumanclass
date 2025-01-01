#!/bin/bash

##################################################################################
qiime rescript get-silva-data \
--p-version '138.1' \
--p-target 'SSURef_NR99' \
--p-include-species-labels \
--o-silva-sequences silva_138.1_SSU_NR99_RNA-seqs.qza \
--o-silva-taxonomy silva_138.1_SSU_NR99_tax.qza

qiime rescript cull-seqs \
--i-sequences silva_138.1_SSU_NR99_RNA-seqs.qza \
--o-clean-sequences silva_138.1_SSU_NR99_RNA-seqs_cleaned.qza

qiime rescript filter-seqs-length-by-taxon \
--i-sequences silva_138.1_SSU_NR99_RNA-seqs_cleaned.qza \
--i-taxonomy silva_138.1_SSU_NR99_tax.qza \
--p-labels Archaea Bacteria \
--p-min-lens 900 1200 \
--o-filtered-seqs silva_138.1_SSU_NR99_RNA-seqs-filtered.qza \
--o-discarded-seqs silva_138.1_SSU_NR99_RNA-seqs-discard.qza

qiime rescript dereplicate \
--i-sequences silva_138.1_SSU_NR99_RNA-seqs-filtered.qza \
--i-taxa silva_138.1_SSU_NR99_tax.qza \
--p-mode 'uniq' \
--o-dereplicated-sequences silva_138.1_SSU_NR99_RNA-seqs_dereplicated_uniq.qza \
--o-dereplicated-taxa silva_138.1_SSU_NR99_tax_dereplicated_uniq.qza

##################################################################################

qiime feature-classifier fit-classifier-naive-bayes \
--i-reference-reads silva_138.1_SSU_NR99_RNA-seqs_dereplicated_uniq.qza \
--i-reference-taxonomy silva_138.1_SSU_NR99_tax_dereplicated_uniq.qza \
--o-classifier silva_138.1_SSU_NR99_classifier.qza

##################################################################################

qiime feature-classifier classify-sklearn \
--i-classifier silva_138.1_SSU_NR99_classifier.qza \
--i-reads merged_sv.qza \
--p-confidence disable --p-n-jobs 10 \
--o-classification SILVA_138.1_merged_classification.qza

##################################################################################

qiime tools export \
--input-path SILVA_138.1_merged_classification.qza \
--output-path ./
mv taxonomy.tsv SILVA_138.1_merged_classification.tsv
grep -v -e "Mitochondria" -e "Chloroplast" SILVA_138.1_merged_classification.tsv > SILVA_138.1_merged_classification_filtered.tsv

qiime tools import \
--input-path SILVA_138.1_merged_classification_filtered.tsv \
--output-path SILVA_138.1_merged_classification_filtered.qza \
--type 'FeatureData[Taxonomy]'

qiime taxa filter-table --i-table merged_processed.qza \
--i-taxonomy SILVA_138.1_merged_classification.qza \
--p-exclude Chloroplast,Mitochondria \
--o-filtered-table SILVA_138.1_merged_processed_filtered.qza

##################################################################################

qiime clawback generate-class-weights \
--i-reference-sequences silva_138.1_SSU_NR99_RNA-seqs_dereplicated_uniq.qza \
--i-reference-taxonomy silva_138.1_SSU_NR99_tax_dereplicated_uniq.qza \
--i-samples SILVA_138.1_merged_processed_filtered.qza \
--i-taxonomy-classification SILVA_138.1_merged_classification_filtered.qza \
--o-class-weight SILVA_138.1_merged_weighted.qza

##################################################################################

qiime feature-classifier fit-classifier-naive-bayes \
--i-reference-reads silva_138.1_SSU_NR99_RNA-seqs_dereplicated_uniq.qza \
--i-reference-taxonomy silva_138.1_SSU_NR99_tax_dereplicated_uniq.qza \
--i-class-weight SILVA_138.1_merged_weighted.qza \
--o-classifier SILVA_138.1_weighted_classifier.qza
