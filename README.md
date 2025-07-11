# ONT_FLYE Workflow

A Nextflow workflow for Oxford Nanopore Technologies (ONT) long-read sequencing quality control and assembly preparation.

## Overview

This workflow performs size-based filtering and quality control on ONT long-read sequencing data:

1. **Size-based filtering** using Chopper - splits reads into length bins
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
