#  Tutorial

## For shotgun sequencing
### 1. Pre-processing
```
fastp fastp \--in1 Sample_1.fastq.gz \
    --in2 Sample_2.fastq.gz \
    --out1 Sample_fastp_1.fastq.gz \
    --out2 Sample_fastp_1.fastq.gz \
    --unpaired1  Sample_fastp_unpaired_1.fastq.gz \
    --unpaired2  Sample__fastp_unpaired_2.fastq.gz \
    --html  Sample_fastp.html \
    --thread 16 --verbose
```
### 2. Filter host or feed ingredient genome
```
# Build bowtie2 index
    bowtie2-build \
    -f SNU_Hanwoo_genome.fna \
    ./SNU_Hanwoo_bowtie2_db/SNU_Hanwoo_genome.btindex \
    --threads 10
```
```
# Mapping using bowtie2
bowtie2 \
    --local -p 40 \
    -x SNU_Hanwoo_bowtie2_db/SNU_Hanwoo_genome.btindex \
    -1 Sample_fastp_unpaired_1.fastq.gz \
    -2 Sample_fastp_unpaired_2.fastq.gz \
    -S Sample_bowtie2.sam 2> ${main_out}_bowtie2_hanwoo.log
samtools view \
    --threads 40 -Sb Sample_bowtie2.sam > Sample_bowtie2.bam
samtools view \
    --threads 40 -bf 0x04 Sample_bowtie2.bam > Sample_Hanwoo_filtered.bam && \
    samtools bam2fq Sample_Hanwoo_filtered.bam > Sample_Hanwoo_filtered.fastq
```
### 3. SortMeRNA for extract rRNA sequence
```
sortmerna \
    --ref /home/ryukseok/_DATABASEs/sortmerna_db/pr2_version_5.0.0_SSU_UTAX.fasta \
    --reads Sample_filtered.fastq \
    --workdir Sample_SortMeRNA \
    --threads 30 \
    --fastx
```
### 4. QIIME2 import and dereplication
```
# Import samples into QIIME2 artifact
qiime tools import \
	--type 'SampleData[SequencesWithQuality]' \
	--input-format SingleEndFastqManifestPhred33V2 \
	--input-path manifest.tsv \
	--output-path shotgun_data.qza

# Initial quality filtering process based on quality scores
qiime quality-filter q-score \
--i-demux shotgun_data.qza \
--o-filtered-sequences shotgun_data_qfiltered.qza \
--o-filter-stats shotgun_data_qfilter-stats.qza

# Perform dereplicate sequences
qiime vsearch dereplicate-sequences \
	--i-sequences shotgun_data_qfiltered.qza \
	--o-dereplicated-table shotgun_data_table.qza \
	--o-dereplicated-sequences shotgun_data_rep-seqs.qza
```

## For amplicon sequencing
### 1. Pre-processing
If your paired-end sequences are already merged, you can directly import them into the workflow.
For paired-end sequences, primers are first removed using **Cutadapt**, and the sequences are then merged using **FLASH2** before being used in the workflow.

### 2. QIIME2 import and denoising
```
# Import samples into QIIME2 artifact
qiime tools import \
	--type 'SampleData[SequencesWithQuality]' \
	--input-format SingleEndFastqManifestPhred33V2 \
	--input-path manifest_flash2.tsv \
	--output-path BioProject.qza

# Initial quality filtering process based on quality scores
qiime quality-filter q-score \
--i-demux BioProject.qza \
--o-filtered-sequences BioProject_qfiltered.qza \
--o-filter-stats BioProject_qfilter-stats.qza

# Perform denoise sequences
qiime deblur denoise-16S \
	--i-demultiplexed-seqs BioProject_qfiltered.qza \
	--p-left-trim-len 0 \
	--p-trim-length 300 \
	--p-sample-stats \
	--p-jobs-to-start 25 \
	--o-table BioProject_table.qza \
	--o-representative-sequences BioProject_rep-seqs.qza \
	--o-stats BioProject_deblur-stats.qza
```

## Merge datasets
```
# Merge tables
qiime feature-table merge --i-tables \
    shotgun_data_A_rep-seqs.qza \
    shotgun_data_B_rep-seqs.qza \
    ...
    BioProject_A_table.qza \
    BioProject_B_table.qza \
    BioProject_C_table.qza \
    ...
--o-merged-table merged_table.qza
```
```
# Merge representative sequences
qiime feature-table merge-seqs --i-data \
    shotgun_data_A_rep-seqs.qza \
    shotgun_data_B_rep-seqs.qza \
    ...
    BioProject_A_rep-seqs.qza \
    BioProject_B_rep-seqs.qza \
    BioProject_C_rep-seqs.qza \
    ...
--o-merged-rep-seqs merged_rep-seqs.qza
```
## Export datasets
```
# Export tables
qiime tools export \
--input-path merged_table.qza
--output-path ./
mv feature-table.biom merged_table.biom
biom convert -i merged_table.biom -o merged_table.tsv --to-tsv 
```
```
# Export representative sequences
qiime tools export \
--input-path merged_rep-seqs.qza
--output-path ./
mv dna-sequences.fasta merged_rep-seqs.fasta
```

## Combined representative sequences and table
**IMPORTANT** | Use the `combined-multiprocessing.py` script available in this repository to merge the table and representative sequences based on feature IDs. While you can use unmerged tables and representative sequences, additional processing steps will be required.
```
python3 combined-multiprocessing.py \
    --input-seq merged_rep-seqs.fasta \
    --input-table merged_table.tsv \
    --output merged_processed.tsv \
    --threads 20
```

## Re-import as QIIME2 artifact and process input format for weighting
```
biom convert -i merged_processed.tsv -o merged_processed.biom
qiime tools import \
    --input-path merged_processed.biom \
    --output-path merged_processed.qza \
    --type 'FeatureTable[Frequency]'
```
```

```
