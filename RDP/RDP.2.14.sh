#!/bin/bash

# Please download Seq and taxonomy file at 
#### https://sourceforge.net/projects/rdp-classifier/files/rdp-classifier/

##################################################################################
qiime tools import \
	--input-path RefOTUs.fa \
	--output-path rdp_ref_seqs.qza \
	--type 'FeatureData[Sequence]' \
	--input-format 'MixedCaseDNAFASTAFormat'

qiime tools import \
	--input-path Ref_ref_taxonomy.txt \
	--output-path rdp_ref_tax.qza \
	--type 'FeatureData[Taxonomy]' \
	--input-format 'HeaderlessTSVTaxonomyFormat'

 qiime rescript cull-seqs \
 --i-sequences RDP_ref_seqs.qza \
 --o-clean-sequences RDP_ref_seqs_cleaned.qza

qiime rescript filter-seqs-length-by-taxon \
--i-sequences RDP_ref_seqs_cleaned.qza \
--i-taxonomy RDP_ref_tax.qza \
--p-labels Archaea Bacteria \
--p-min-lens 900 1200 \
--o-filtered-seqs RDP_ref_seqs_filtered.qza \
--o-discarded-seqs RDP_ref_seqs_discard.qza 

qiime rescript dereplicate \
--i-sequences RDP_ref_seqs_filtered.qza \
--i-taxa RDP_ref_tax.qza \
--p-mode 'uniq' \
--o-dereplicated-sequences RDP_ref_seqs_dereplicated_uniq.qza \
--o-dereplicated-taxa RDP_ref_tax_dereplicated_uniq.qza

##################################################################################

qiime feature-classifier fit-classifier-naive-bayes \
--i-reference-reads RDP_ref_seqs_dereplicated_uniq.qza \
--i-reference-taxonomy RDP_ref_tax_dereplicated_uniq.qza \
--o-classifier RDP_classifier.qza

##################################################################################

qiime feature-classifier classify-sklearn \
--i-classifier RDP_classifier.qza \
--i-reads merged_sv.qza \
--p-confidence disable --p-n-jobs 10 \
--o-classification RDP_merged_classification.qza

##################################################################################

qiime tools export \
--input-path RDP_merged_classification.qza \
--output-path ./
mv taxonomy.tsv RDP_merged_classification.tsv
grep -v -e "Mitochondria" -e "Chloroplast" RDP_merged_classification.tsv > RDP_merged_classification_filtered.tsv

qiime tools import \
--input-path RDP_merged_classification_filtered.tsv \
--output-path RDP_merged_classification_filtered.qza \
--type 'FeatureData[Taxonomy]'

qiime taxa filter-table --i-table merged_processed.qza \
--i-taxonomy RDP_merged_classification.qza \
--p-exclude Chloroplast,Mitochondria \
--o-filtered-table RDP_merged_processed_filtered.qza

##################################################################################

qiime clawback generate-class-weights \
--i-reference-sequences RDP_ref_seqs_dereplicated_uniq.qza \
--i-reference-taxonomy RDP_ref_tax_dereplicated_uniq.qza \
--i-samples RDP_merged_processed_filtered.qza \
--i-taxonomy-classification RDP_merged_classification_filtered.qza \
--o-class-weight RDP_merged_weighted.qza

##################################################################################

qiime feature-classifier fit-classifier-naive-bayes \
--i-reference-reads RDP_ref_seqs_dereplicated_uniq.qza \
--i-reference-taxonomy RDP_ref_tax_dereplicated_uniq.qza \
--i-class-weight RDP_merged_weighted.qza \
--o-classifier RDP_weighted_classifier.qza
