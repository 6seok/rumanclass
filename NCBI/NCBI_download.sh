#!/bin/bash

##################################################################################
qiime rescript get-ncbi-data \
--p-query '33175[BioProject] OR 33317[BioProject]' \
--o-sequences NCBI_RefSeqs_RNA.qza \
--o-taxonomy NCBI_RefSeqs_tax.qza

qiime rescript cull-seqs \
--i-sequences NCBI_RefSeqs_RNA.qza \
--o-clean-sequences NCBI_RefSeqs_RNA_cleaned.qza

qiime rescript filter-seqs-length-by-taxon \
--i-sequences NCBI_RefSeqs_RNA_cleaned.qza \
--i-taxonomy NCBI_RefSeqs_tax.qza \
--p-labels Archaea Bacteria \
--p-min-lens 900 1200 \
--o-filtered-seqs NCBI_RefSeqs_RNA_filtered.qza \
--o-discarded-seqs NCBI_RefSeqs_RNA_discard.qza

qiime rescript dereplicate \
--i-sequences NCBI_RefSeqs_RNA_filtered.qza \
--i-taxa NCBI_RefSeqs_tax.qza \
--p-mode 'uniq' \
--o-dereplicated-sequences NCBI_RefSeqs_RNA_dereplicated_uniq.qza \
--o-dereplicated-taxa NCBI_RefSeqs_tax_dereplicated_uniq.qza

##################################################################################

qiime feature-classifier fit-classifier-naive-bayes \
--i-reference-reads NCBI_RefSeqs_RNA_dereplicated_uniq.qza \
--i-reference-taxonomy NCBI_RefSeqs_tax_dereplicated_uniq.qza \
--o-classifier NCBI_RefSeqs_classifier.qza

##################################################################################

qiime feature-classifier classify-sklearn \
--i-classifier NCBI_RefSeqs_classifier.qza \
--i-reads merged_sv.qza \
--p-confidence disable --p-n-jobs 10 \
--o-classification NCBI_merged_classification.qza

##################################################################################

qiime tools export \
--input-path NCBI_merged_classification.qza \
--output-path 
mv taxonomy.tsv NCBI_merged_classification.tsv
grep -v -e "Mitochondria" -e "Chloroplast" NCBI_merged_classification.tsv > NCBI_merged_classification_filtered.tsv

qiime tools import \
--input-path NCBI_merged_classification_filtered.tsv \
--output-path NCBI_merged_classification_filtered.qza \
--type 'FeatureData[Taxonomy]'

qiime taxa filter-table --i-table merged_processed.qza \
--i-taxonomy NCBI_merged_classification.qza \
--p-exclude Chloroplast,Mitochondria \
--o-filtered-table NCBI_merged_processed_filtered.qza

##################################################################################

qiime clawback generate-class-weights \
--i-reference-sequences NCBI_RefSeqs_RNA_dereplicated_uniq.qza \
--i-reference-taxonomy NCBI_RefSeqs_tax_dereplicated_uniq.qza \
--i-samples NCBI_merged_processed_filtered.qza \
--i-taxonomy-classification NCBI_merged_classification_filtered.qza \
--o-class-weight NCBI_merged_weighted.qza

##################################################################################

qiime feature-classifier fit-classifier-naive-bayes \
--i-reference-reads NCBI_RefSeqs_RNA_dereplicated_uniq.qza \
--i-reference-taxonomy NCBI_RefSeqs_tax_dereplicated_uniq.qza \
--i-class-weight NCBI_merged_weighted.qza \
--o-classifier NCBI_weighted_classifier.qza
