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
1. Pre-processing
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
2. Filter host and feed ingredient genome
```
## Build bowtie2 index
bowtie2-build -f SNU_Hanwoo_genome.fna ./SNU_Hanwoo_bowtie2_db/SNU_Hanwoo_genome.btindex --threads 10
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
3. SortMeRNA for extract rRNA sequence
```

```
4. QIIME2 import
```

```



