# ONT Flye Workflow

A Nextflow workflow for assembling integrant-containing yeast genomes with Oxford Nanopore Technologies (ONT) long-read sequence data and evaluation of transgene copy numbers in assemblies.

## Overview

This workflow provides three distinct analysis modes for processing ONT long-read sequencing data for different phases of analysis:

### 🔍 **Scan Mode**

Performs quality-based filtering and applies read-length thresholds to identify minimal read length requirements for assembly.

### 📊 **Downsample Mode**

Creates random subsets of reads at different coverage levels to assess assembly quality vs. sequencing depth.

### 🔄 **Bootstrap Mode**

Generates multiple replicate assemblies from random subsamples to evaluate assembly consistency and robustness.

## Workflow Modes

### Scan Mode

**Purpose**: Identify optimal read length thresholds for your dataset

**Process**:

1. Filters reads using Chopper to quality ≥10
2. Creates three size-filtered datasets:
   - Reads ≥30kb
   - Reads ≥40kb  
   - Reads ≥50kb
3. Assembles each filtered dataset with Flye
4. Performs transgene BLAST analysis
5. Generates comprehensive QC reports with NanoPlot

**Best for**: Initial dataset exploration and parameter optimization

### Downsample Mode

**Purpose**: Evaluate the relationship between sequencing depth and assembly quality

**Process**:

1. Creates random subsets at 25%, 50%, and 75% of original reads
2. Assembles each subset with Flye
3. Performs transgene BLAST analysis
4. Compares assembly metrics across coverage levels

**Best for**: Determining minimum sequencing depth requirements

### Bootstrap Mode

**Purpose**: Assess assembly reproducibility

**Process**:

1. Creates multiple (default: 10) random subsamples at specified fraction (default: 75%)
2. Assembles each replicate independently
3. Performs transgene BLAST analysis on all replicates
4. Enables statistical analysis of assembly variation

**Best for**: Production of confidence intervals

## Requirements

### Software Dependencies

- **Nextflow** (≥22.10.0)
- **Conda/Mamba** (for environment management)
- **Singularity** (recommended for HPC environments)

### Computational Resources

- **CPU**: 24+ cores recommended for Flye assembly
- **Memory**: 32+ GB RAM for large genomes
- **Storage**: ~10x input file size for intermediate files

## Installation

### Clone the repository

```bash
git clone https://github.com/KochInstitute-Bioinformatics/ont-flye-workflow.git
cd ont-flye-workflow
```

### Verify Nextflow installation

```bash
nextflow -version
```

## Usage

### Single Sample Analysis

```bash
# Scan mode - explore optimal read lengths
nextflow run main.nf \
  --input_fastq sample.fastq.gz \
  --name sample_001 \
  --mode scan \
  --outdir results

# Downsample mode - test coverage requirements  
nextflow run main.nf \
  --input_fastq sample.fastq.gz \
  --name sample_001 \
  --mode downsample \
  --outdir results

# Bootstrap mode - assess reproducibility
nextflow run main.nf \
  --input_fastq sample.fastq.gz \
  --name sample_001 \
  --mode bootstrap \
  --outdir results
```

### Batch Processing with CSV Input

```bash
# Process multiple samples with different modes
nextflow run main.nf --samples batch_samples.csv --outdir results
```

### Input Formats

#### CSV File Format for batch processing

Create a CSV file with the following columns:

```csv
name,fastq,transgene,mode
sample_001,/path/to/sample_001.fastq.gz,A-vector_herceptin_pEY345,scan
sample_002,/path/to/sample_002.fastq.gz,A-vector_herceptin_pEY345,downsample
sample_003,/path/to/sample_003.fastq.gz,A-vector_herceptin_pEY345,bootstrap
```

**Column descriptions**:

- `name`: Unique sample identifier
- `fastq`: Full path to FASTQ file
- `transgene`: Transgene name (must match files in `transgenes/` directory)
- `mode`: Analysis mode (`scan`, `downsample`, or `bootstrap`)

### Single sample processing parameters

| Parameter | Description | Example |
|-----------|-------------|---------|
| `--input_fastq` | Path to input FASTQ file (single sample mode) | `sample.fastq.gz` |
| `--name` | Sample name (single sample mode) | `sample_001` |
| `--samples` | Path to CSV file (batch mode) | `samples.csv` |

### Optional Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--outdir` | `results` | Output directory |
| `--mode` | `scan` | Analysis mode (single sample) |
| `--genome_size` | `10m` | Estimated genome size for Flye |
| `--min_quality` | `10` | Minimum read quality score |
| `--transgene` | `A-vector_herceptin_pEY345` | Transgene name |

### Mode-Specific Parameters

#### Scan Mode Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--size_ranges` | `[30k, 40k, 50k]` | Read length thresholds |

#### Downsample Mode Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--downsample_rates` | `[0.25, 0.5, 0.75]` | Fraction of reads to retain |

#### Bootstrap Mode Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--fraction` | `0.75` | Fraction of reads per replicate |
| `--replicates` | `10` | Number of bootstrap replicates |

## Output Structure

```text
results/
├── pipeline_info/           # Execution reports and logs
├── chopper/                 # Filtered reads (scan mode)
├── nanoplot/               # Quality control reports
├── flye/                   # Assembly outputs
├── blast/                  # Transgene analysis results
├── downsample/             # Downsampled reads (downsample mode)
├── bootstrap/              # Bootstrap replicates (bootstrap mode)
└── summary/                # Aggregated statistics
```

### Key Output Files

- **Assembly FASTA**: `flye/{sample}/assembly.fasta`
- **Assembly statistics**: `summary/assembly_stats_summary.csv`
- **QC reports**: `nanoplot/{sample}/NanoPlot-report.html`
- **Transgene analysis**: `blast/{sample}_blast_results.txt`

## Transgene Configuration

The workflow includes a transgene library in the `transgenes/` directory. To add new transgenes:

1. Add FASTA file to `transgenes/` directory
2. Update `transgenes/transgene_library.csv`
3. Specify transgene name in your input CSV or command line

## Performance Optimization

### HPC Configuration

### Example SLURM execution

In process - this needs to be verified and tested!

```bash
nextflow run main.nf --samples batch.csv -profile slurm -c custom.config
```

### Resource Tuning

Adjust resources in `nextflow.config` based on your system:

```text
process {
    withName: 'FLYE' {
        cpus = 48        // Increase for faster assembly
        memory = '64 GB' // Increase for large genomes
        time = '24.h'    // Adjust based on genome size
    }
}
```

## Troubleshooting

### Common Issues

**Assembly failures**:

- Check read quality and quantity
- Verify genome size estimate
- Increase memory allocation

**Missing transgene files**:

- Ensure transgene FASTA exists in `transgenes/` directory
- Check transgene name spelling

**Resource errors**:

- Adjust CPU/memory in `nextflow.config`
- Use appropriate execution profile

### Getting Help

Contact: [Charlie Whittaker](mailto:charliew@mit.edu)

### Display help message

nextflow run main.nf --help

### Check workflow version

nextflow run main.nf --version

## Citation

If you use this workflow in your research, please cite:

[citation]

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## Contributing

Contributions are welcome! Please feel free to submit a Pull Request.

## Support

For questions and support, please open an issue on the [GitHub repository](https://github.com/KochInstitute-Bioinformatics/ont-flye-workflow/issues).
