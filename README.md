# ONT_FLYE Workflow

A Nextflow workflow for assembling integrant-containing yeast genomes with Oxford Nanopore Technologies (ONT) long-read sequence data and evaluation of number of integrent copies in assembly.

## Overview

This workflow has 3 different modes. The input data are the critical difference between the modes

1. **Scan mode** uses Chopper to filter input fastq files to q10+ reads that exceed 3 different length thresholds:
    a) > 30 Kb
    b) > 40 Kb
    c) > 50 Kb

Once sequence subsets are created the assembly routine is executed.

2. **Downsample mode** uses the python script downsample_fastq.py to select random subsets of an input fastq that are 25%, 50%, and 75% of the original. Once subsets are collected
2. **Quality control plots** using NanoPlot - generates comprehensive QC reports
3. **Read length binning** for optimized downstream analysis

## Requirements

### Software

- Nextflow (≥22.10.0)
- Conda/Mamba
- Singularity (for HPC environments)

### Conda Environment

This workflow uses a pre-existing conda environment containing:

- `chopper` - for read filtering
- `NanoPlot` - for quality control visualization

## Quick Start

```bash
# Clone the repository
git clone https://github.com/KochInstitute-Bioinformatics/ont-flye-workflow.git
cd ont-flye-workflow

# Run the workflow with your data
nextflow run main.nf --input_fastq your_reads.fastq --outdir results
