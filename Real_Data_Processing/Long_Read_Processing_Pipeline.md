# Long Read RNA-Seq Processing Pipeline

This document describes the complete pipeline for processing long-read RNA-seq data from FASTQ files to transcript quantification using multiple quantification tools.

## Overview

The pipeline processes human genome data starting from long-read FASTQ files (Oxford Nanopore or PacBio), using 4 different quantification software tools:
- Cufflinks
- StringTie
- featureCounts
- IsoQuant

## Prerequisites

### Required Software
- Minimap2
- SAMtools
- Cufflinks
- StringTie
- featureCounts (subread package)
- IsoQuant

### Required Data Files
- Human reference genome (GRCh38)
- GTF annotation file (Homo_sapiens.GRCh38.112.gtf)
- Long-read FASTQ/FASTA files (Oxford Nanopore or PacBio)

## Pipeline Steps

### Step 1: Genome Preparation

```bash
# Use the same merged chromosome file from short-read pipeline
# This step combines all human chromosomes (1-22, X, Y) from GRCh38 into a single FASTA file
cd index
zcat Homo_sapiens.GRCh38.dna.chromosome.{1..22}.fa.gz \
     Homo_sapiens.GRCh38.dna.chromosome.X.fa.gz \
     Homo_sapiens.GRCh38.dna.chromosome.Y.fa.gz > merged_chromosomes.fa
```

**Purpose**: Creates a unified reference genome file containing all chromosomes for alignment.

**Note**: If you have already prepared this file for short-read processing, you can skip this step.

### Step 2: Long Read Alignment

```bash
# Align long reads using minimap2 with splice-aware mode
minimap2 -ax splice ...path/to/index/merged_chromosomes.fa ENCFFXXXXXX.fasta > human_long.sam

# Sort and index the alignment file
samtools sort -@ 8 -o human_long.bam human_long.sam
samtools index human_long.bam
```

**Minimap2 Parameters**:
- `-a`: Output in SAM format
- `-x splice`: Preset for spliced alignment (long cDNA reads)
- Input: Reference genome and long-read FASTA file
- Output: SAM alignment file

**SAMtools Parameters**:
- `-@ 8`: Use 8 threads for sorting
- `-o`: Output BAM file

**Purpose**: Aligns long reads to the reference genome with splice-aware mapping optimized for RNA-seq data.

### Step 3: Quantification with Multiple Tools

#### 3.1 Cufflinks Quantification

```bash
cufflinks -p 8 -g ...path/to/human_gtf/Homo_sapiens.GRCh38.112.gtf -o cufflinks/ human_long.bam
```

**Parameters**:
- `-p 8`: Use 8 threads
- `-g`: GTF annotation file (for guided assembly)
- `-o`: Output directory

**Output**: Transcript abundance estimates in FPKM values.

#### 3.2 featureCounts Quantification

```bash
featureCounts -a ...path/to/human_gtf/Homo_sapiens.GRCh38.112.gtf \
              -o ans_fea/my_featureCounts_output.txt \
              -L \
              -t exon \
              -g transcript_id \
              human_long.bam
```

**Parameters**:
- `-a`: GTF annotation file
- `-o`: Output count file
- `-L`: Count long reads (for reads longer than 1kb)
- `-t exon`: Feature type to count (exons)
- `-g transcript_id`: Attribute type for grouping (transcript level)

**Output**: Raw read counts per transcript, optimized for long reads.

#### 3.3 IsoQuant Quantification

```bash
isoquant.py --reference ...path/to/index/merged_chromosomes.fa \
            --genedb ...path/to/human_gtf/Homo_sapiens.GRCh38.112.gtf \
            --fastq ENCFFXXXXXX.fasta \
            --data_type nanopore \
            -o is/isoquant_output
```

**Parameters**:
- `--reference`: Reference genome FASTA file
- `--genedb`: GTF annotation file
- `--fastq`: Input long-read FASTA file
- `--data_type nanopore`: Specify data type (use 'pacbio' for PacBio data)
- `-o`: Output directory

**Output**: Comprehensive transcript quantification and isoform analysis specifically designed for long reads.

#### 3.4 StringTie Quantification

```bash
stringtie -L -e -p 8 \
          -G ...path/to/human_gtf/Homo_sapiens.GRCh38.112.gtf \
          -o stringtie/human_long_stringtie.gtf \
          human_long.bam
```

**Parameters**:
- `-L`: Disable trimming of predicted transcripts based on coverage (important for long reads)
- `-e`: Estimate abundances of given reference transcripts only
- `-p 8`: Use 8 threads
- `-G`: Reference annotation GTF file
- `-o`: Output GTF file with abundance estimates

**Output**: GTF file with transcript abundance estimates optimized for long reads.

## Output Files Summary

After running the complete pipeline, you will have:

1. **Cufflinks**: `cufflinks/transcripts.gtf`, `cufflinks/isoforms.fpkm_tracking`
2. **featureCounts**: `ans_fea/my_featureCounts_output.txt`
3. **IsoQuant**: `is/isoquant_output/` (multiple output files including transcript counts)
4. **StringTie**: `stringtie/human_long_stringtie.gtf`
