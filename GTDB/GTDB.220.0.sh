#!/bin/bash/

##################################################################################
qiime rescript get-gtdb-data \
    --p-version '220.0' \
    --p-db-type 'SpeciesReps' \
    --p-domain 'Both' \
    --o-gtdb-sequences GTDB-220.0_RNA.qza \
    --o-gtdb-taxonomy GTDB_220.0_tax.qza \
    --verbose

qiime rescript cull-seqs \
    --i-sequences GTDB_220.0_RNA.qza \
    --o-clean-sequences GTDB_220.0_RNA_cleaned.qza

qiime rescript filter-seqs-length-by-taxon \
    --i-sequences GTDB_220.0_RNA_cleaned.qza \
    --i-taxonomy GTDB_220.0_tax.qza \
    --p-labels Archaea Bacteria \
    --p-min-lens 900 1200 \
    --o-filtered-seqs GTDB_220.0_RNA_filtered.qza \
    --o-discarded-seqs GTDB_220.0_RNA_discard.qza 

qiime rescript dereplicate \
    --i-sequences GTDB_220.0_RNA_filtered.qza \
    --i-taxa GTDB_220.0_tax.qza \
    --p-mode 'uniq' \
    --o-dereplicated-sequences GTDB_220.0_RNA_dereplicated_uniq.qza \
    --o-dereplicated-taxa GTDB_220.0_tax_dereplicated_uniq.qza

##################################################################################

qiime feature-classifier fit-classifier-naive-bayes \
    --i-reference-reads GTDB_220.0_RNA_dereplicated_uniq.qza \
    --i-reference-taxonomy GTDB_220.0_tax_dereplicated_uniq.qza \
    --o-classifier GTDB_220.0_classifier.qza

##################################################################################

qiime feature-classifier classify-sklearn \
--i-classifier GTDB_220.0_classifier.qza \
--i-reads merged_sv.qza \
--p-confidence disable --p-n-jobs 10 \
--o-classification GTDB_220.0_merged_classification.qza

##################################################################################

qiime tools export \
--input-path GTDB_220.0_merged_classification.qza \
--output-path ./
mv taxonomy.tsv GTDB_220.0_merged_classification.tsv
grep -v -e "Mitochondria" -e "Chloroplast" GTDB_220.0_merged_classification.tsv > GTDB_220.0_merged_classification_filtered.tsv

qiime tools import \
--input-path GTDB_220.0_merged_classification_filtered.tsv \
--output-path GTDB_220.0_merged_classification_filtered.qza \
--type 'FeatureData[Taxonomy]'

qiime taxa filter-table --i-table merged_processed.qza \
--i-taxonomy GTDB_220.0_merged_classification.qza \
--p-exclude Chloroplast,Mitochondria \
--o-filtered-table GTDB_220.0_merged_processed_filtered.qza

##################################################################################

qiime clawback generate-class-weights \
--i-reference-sequences GTDB_220.0_SSU_NR99_RNA-seqs_dereplicated_uniq.qza \
--i-reference-taxonomy GTDB_220.0_SSU_NR99_tax_dereplicated_uniq.qza \
--i-samples GTDB_220.0_merged_processed_filtered.qza \
--i-taxonomy-classification GTDB_220.0_merged_classification_filtered.qza \
--o-class-weight GTDB_220.0_merged_weighted.qza

##################################################################################

qiime feature-classifier fit-classifier-naive-bayes \
--i-reference-reads GTDB_220.0_SSU_NR99_RNA-seqs_dereplicated_uniq.qza \
--i-reference-taxonomy GTDB_220.0_SSU_NR99_tax_dereplicated_uniq.qza \
--i-class-weight GTDB_220.0_merged_weighted.qza \
--o-classifier GTDB_220.0_weighted_classifier.qza
