#!/bin/bash

# Get greengenes data
## https://ftp.microbio.me/greengenes_release/current/
wget https://ftp.microbio.me/greengenes_release/current/2022.10.backbone.full-length.fna.qza  --no-check-certificate
wget https://ftp.microbio.me/greengenes_release/current/2022.10.backbone.tax.qza  --no-check-certificate
mv 2022.10.backbone.full-length.fna.qza greengenes_2022.10.backbone.full-length.RNA.qza
mv 2022.10.backbone.tax.qza greengenes_2022.10.backbone.tax.qza

##################################################################################

qiime feature-classifier fit-classifier-naive-bayes \
--i-reference-reads greengenes_2022.10.backbone.full-length.RNA.qza \
--i-reference-taxonomy greengenes_2022.10.backbone.tax.qza \
--o-classifier greengenes_2022.10_classifier.qza

##################################################################################

qiime feature-classifier classify-sklearn \
--i-classifier greengenes_2022.10_classifier.qza \
--i-reads merged_sv.qza \
--p-confidence disable --p-n-jobs 10 \
--o-classification greengenes_2022.10_merged_classification.qza

##################################################################################

qiime tools export \
--input-path greengenes_2022.10_merged_classification.qza \
--output-path ./
mv taxonomy.tsv greengenes_2022.10_merged_classification.tsv
grep -v -e "Mitochondria" -e "Chloroplast" greengenes_2022.10_merged_classification.tsv > greengenes_2022.10_merged_classification_filtered.tsv

qiime tools import \
--input-path greengenes_2022.10_merged_classification_filtered.tsv \
--output-path greengenes_2022.10_merged_classification_filtered.qza \
--type 'FeatureData[Taxonomy]'

qiime taxa filter-table --i-table merged_processed.qza \
--i-taxonomy greengenes_2022.10_merged_classification.qza \
--p-exclude Chloroplast,Mitochondria \
--o-filtered-table greengenes_2022.10_merged_processed_filtered.qza

##################################################################################

qiime clawback generate-class-weights \
--i-reference-sequences greengenes_2022.10.backbone.full-length.RNA.qza \
--i-reference-taxonomy greengenes_2022.10.backbone.tax.qza \
--i-samples greengenes_2022.10_merged_processed_filtered.qza \
--i-taxonomy-classification greengenes_2022.10_merged_classification_filtered.qza \
--o-class-weight greengenes_2022.10_merged_weighted.qza

##################################################################################

qiime feature-classifier fit-classifier-naive-bayes \
--i-reference-reads greengenes_2022.10.backbone.full-length.RNA.qza \
--i-reference-taxonomy greengenes_2022.10.backbone.tax.qza \
--i-class-weight greengenes_2022.10_merged_weighted.qza \
--o-classifier greengenes_2022.10_weighted_classifier.qza
