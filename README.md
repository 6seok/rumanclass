# RUMANCLASS (RUminal microbiome MANually CLASSifier)

## QIIME2 manually weighted classifier for species-specific ruminal e micribome analysis

This method was developed to provide enhanced resolution for ruminant microbiome analysis. It is an advanced version of the approach proposed by Kaehler et al. (2019). The method utilizes open-source sequences available from the NCBI BioProject respiratory database and employs preprocessing steps alongside the QIIME2 default workflow.

For detailed instructions, please refer to the **[Tutorial](https://github.com/6seok/rumanclass/tree/main#tutorial)**.

As of January 2025, a species-specific weighted classifier for Hanwoo cattle has been created. 
Future plans include developing classifiers for Holstein, Jersey, Charolais, and Aberdeen Angus breeds.

***

**If you have referenced this workflow or used the pre-built weighted classifier, please include the following citation:**
+ R Kang and T Park (2025), Manually weighted classifier achieves higher resolution in species-specific ruminal microbi-ome analysis compared to standard or average weighted classifiers ,  [UNDECISION], https://doi.org/UNDECISION
***

##  Tutorial

### Workflow
![image](https://github.com/user-attachments/assets/ee4b9abe-bac2-44ff-b10f-a3ed6defe993)

### For shotgun sequencing
#### 1. Pre-processing
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
#### 2. Filter host or feed ingredient genome
```
## Build bowtie2 index
    bowtie2-build \
    -f SNU_Hanwoo_genome.fna \
    ./SNU_Hanwoo_bowtie2_db/SNU_Hanwoo_genome.btindex \
    --threads 10
```
```
## Mapping using bowtie2
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
#### 3. SortMeRNA for extract rRNA sequence
```
sortmerna \
    --ref /home/ryukseok/_DATABASEs/sortmerna_db/pr2_version_5.0.0_SSU_UTAX.fasta \
    --reads Sample_filtered.fastq \
    --workdir Sample_SortMeRNA \
    --threads 30 \
    --fastx
```
#### 4. QIIME2 import and dereplication
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

### For amplicon sequencing
#### 1. Pre-processing
If your paired-end sequences are already merged, you can directly import them into the workflow.
For paired-end sequences, primers are first removed using **Cutadapt**, and the sequences are then merged using **FLASH2** before being used in the workflow.

#### 2. QIIME2 import and denoising
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

### Merge datasets
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

