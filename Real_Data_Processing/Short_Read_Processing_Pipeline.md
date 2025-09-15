# Short Read RNA-Seq Processing Pipeline

This document describes the complete pipeline for processing short-read RNA-seq data from FASTQ files to transcript quantification using multiple quantification tools.

## Overview

The pipeline processes human genome data starting from FASTQ files, using 4 different quantification software tools:
- Cufflinks
- StringTie
- featureCounts
- RSEM

## Prerequisites

### Required Software
- HISAT2
- SAMtools
- Cufflinks
- StringTie
- featureCounts (subread package)
- RSEM
- Bowtie2

### Required Data Files
- Human reference genome (GRCh38)
- GTF annotation file (Homo_sapiens.GRCh38.112.gtf)
- Paired-end FASTQ files

## Pipeline Steps

### Step 1: Genome Preparation

```bash
# Prepare merged chromosome file
# This step combines all human chromosomes (1-22, X, Y) from GRCh38 into a single FASTA file
cd index
zcat Homo_sapiens.GRCh38.dna.chromosome.{1..22}.fa.gz \
     Homo_sapiens.GRCh38.dna.chromosome.X.fa.gz \
     Homo_sapiens.GRCh38.dna.chromosome.Y.fa.gz > merged_chromosomes.fa
```

**Purpose**: Creates a unified reference genome file containing all chromosomes for alignment.

### Step 2: Extract Exon Information

```bash
cd human_gtf
# Extract exon information from GTF file
# This script processes the GTF file to extract exon coordinates and information
python extract_exons.py Homo_sapiens.GRCh38.112.gtf > human.exon
```

**Purpose**: Extracts exon annotations from the GTF file, which will be used for splice-aware alignment.

### Step 3: Build HISAT2 Index

```bash
# Build HISAT2 index with splice site and exon information
# This creates an index optimized for RNA-seq alignment with known splice junctions
hisat2-build --ss human_gtf/human.ss --exon human_gtf/human.exon index/merged_chromosomes.fa human_tran
```

**Parameters**:
- `--ss`: Splice sites file 
- `--exon`: Exon information file
- `human_tran`: Output index prefix

**Purpose**: Creates a splice-aware index that enables accurate alignment across exon-exon junctions.

### Step 4: Read Alignment and Processing

```bash
# Align paired-end reads using HISAT2
hisat2 -p 8 --dta -x ...path/to/human_tran -1 ENCFFXXXXX.fastq.gz -2 ENCFFXXXXX.fastq.gz -S human_short.sam

# Sort and index the alignment file
samtools sort -@ 8 -o human_short.bam human_short.sam
samtools index human_short.bam
```

**HISAT2 Parameters**:
- `-p 8`: Use 8 threads
- `--dta`: Downstream transcriptome analysis mode (optimized for StringTie/Cufflinks)
- `-x`: Path to index
- `-1/-2`: Paired-end FASTQ files
- `-S`: Output SAM file

**SAMtools Parameters**:
- `-@ 8`: Use 8 threads for sorting
- `-o`: Output BAM file

### Step 5: Quantification with Multiple Tools

#### 5.1 Cufflinks Quantification

```bash
cufflinks -p 8 -g ...path/to/human_gtf/Homo_sapiens.GRCh38.112.gtf -o ans/ human_short.bam
```

**Parameters**:
- `-p 8`: Use 8 threads
- `-g`: GTF annotation file (for guided assembly)
- `-o`: Output directory

**Output**: Transcript abundance estimates in FPKM values.

#### 5.2 StringTie Quantification

```bash
stringtie -e -p 8 -G ...path/to/human_gtf/Homo_sapiens.GRCh38.112.gtf -o ans/human_short_stringtie.gtf human_short.bam
```

**Parameters**:
- `-e`: Estimate abundances of given reference transcripts only
- `-p 8`: Use 8 threads
- `-G`: Reference annotation GTF file
- `-o`: Output GTF file with abundance estimates

**Output**: GTF file with transcript abundance estimates.

#### 5.3 featureCounts Quantification

```bash
featureCounts -a ...path/to/human_gtf/Homo_sapiens.GRCh38.112.gtf \
              -o ans/my_featureCounts_output.txt \
              -T 8 \
              -p \
              -t exon \
              -g transcript_id \
              --countReadPairs \
              human_short.bam
```

**Parameters**:
- `-a`: GTF annotation file
- `-o`: Output count file
- `-T 8`: Use 8 threads
- `-p`: Count paired-end reads
- `-t exon`: Feature type to count (exons)
- `-g transcript_id`: Attribute type for grouping (transcript level)
- `--countReadPairs`: Count read pairs instead of individual reads

**Output**: Raw read counts per transcript.

#### 5.4 RSEM Quantification

```bash
# Prepare RSEM reference
rsem-prepare-reference --gtf ...path/to/human_gtf/Homo_sapiens.GRCh38.112.gtf \
                       --bowtie2 \
                       -p 8 \
                       ...path/to/index/merged_chromosomes.fa \
                       human_ensembl

# Calculate expression
rsem-calculate-expression -p 8 \
                         --bowtie2 \
                         --paired-end \
                         ...path/to/short_fq/human_short_1.fastq.gz \
                         ...path/to/short_fq/human_short_2.fastq.gz \
                         human_ensembl \
                         rsem_ans
```

**RSEM-prepare-reference Parameters**:
- `--gtf`: GTF annotation file
- `--bowtie2`: Use Bowtie2 for alignment
- `-p 8`: Use 8 threads
- `human_ensembl`: Reference name prefix

**RSEM-calculate-expression Parameters**:
- `-p 8`: Use 8 threads
- `--bowtie2`: Use Bowtie2 aligner
- `--paired-end`: Input is paired-end reads
- `rsem_ans`: Output prefix

**Output**: Transcript and gene-level abundance estimates in TPM and expected counts.

## Output Files Summary

After running the complete pipeline, you will have:

1. **Cufflinks**: `ans/transcripts.gtf`, `ans/isoforms.fpkm_tracking`
2. **StringTie**: `ans/human_short_stringtie.gtf`
3. **featureCounts**: `ans/my_featureCounts_output.txt`
4. **RSEM**: `rsem_ans.genes.results`, `rsem_ans.isoforms.results`
