qiime rescript get-gtdb-data \
    --p-version '220.0' \
    --p-db-type 'SpeciesReps' \
    --p-domain 'Both' \
    --o-gtdb-sequences tempGTDB-220.0_RNA.qza \
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

qiime feature-classifier fit-classifier-naive-bayes \
--i-reference-reads GTDB_220.0_RNA_dereplicated_uniq.qza \
--i-reference-taxonomy GTDB_220.0_tax_dereplicated_uniq.qza \
--o-classifier GTDB_220.0_classifier.qza
