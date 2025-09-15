# Real Human Genome Data Processing Pipeline

This directory contains detailed documentation for processing real human genome RNA-seq data using standard bioinformatics tools and multiple quantification methods.

## Contents

- `Short_Read_Processing_Pipeline.md` - Complete pipeline for short-read RNA-seq processing
- `Long_Read_Processing_Pipeline.md` - Complete pipeline for long-read RNA-seq processing

## Overview

These pipelines are designed for processing real human genome RNA-seq data and demonstrate how to:

1. **Short-read sequencing** - Process Illumina paired-end reads
2. **Long-read sequencing** - Process Oxford Nanopore or PacBio reads

Each pipeline uses 4 different quantification software tools to provide comprehensive transcript abundance estimation and comparison.

## Pipeline Flow

```
FASTQ Files
     ↓
Genome Preparation & Indexing
     ↓
Read Alignment
     ↓
Multiple Quantification Tools
```

**Note**: These processing pipelines are independent of the main TCCN condition number calculation method. 

## Quick Start

1. Follow the appropriate pipeline documentation (short-read or long-read)
2. Ensure all required software is installed
3. Prepare reference genome and annotation files
4. Run the complete pipeline

## Software Requirements

### For Short-read Processing
- HISAT2 (v2.2.1)
- SAMtools (v1.15)
- Cufflinks (v2.2.1)
- StringTie (v2.2.1)
- featureCounts (subread v2.0.2)
- RSEM (v1.3.3)

### For Long-read Processing
- Minimap2 (v2.24)
- SAMtools (v1.15)
- Cufflinks (v2.2.1)
- StringTie (v2.2.1)
- featureCounts (subread v2.0.2)
- IsoQuant (v3.3.0)

## Data Requirements

- Human reference genome (GRCh38)
- GTF annotation file (Ensembl release 112 recommended)
- Quality-controlled FASTQ files

## Support

For questions regarding the real data processing pipelines, please refer to the specific pipeline documentation or contact the authors.
