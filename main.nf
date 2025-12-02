#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

// Include the main workflow
include { ONT_FLYE } from './workflows/ont_flye'

// Parameters with defaults
params.help = false
params.input_fastq = null
params.name = null
params.samples = null
params.outdir = "results"

// Help message function
def helpMessage() {
    log.info"""
    Usage:
    
    Option 1 - Single sample with runtime naming:
    nextflow run main.nf --input_fastq <path_to_fastq> --name <sample_name>
    
    Option 2 - Multiple samples from CSV:
    nextflow run main.nf --samples <path_to_csv>
    
    Required arguments (Option 1):
    --input_fastq    Path to input FASTQ file
    --name           Sample name for output files
    
    Required arguments (Option 2):
    --samples        Path to CSV file with columns: name,fastq
    
    Optional arguments:
    --outdir         Output directory (default: results)
    --min_quality    Minimum quality score (default: 10)
    --genome_size    Estimated genome size for Flye (default: 10m)
    --help           Show this help message
    
    CSV format example:
    name,fastq
    sHF171,/path/to/sHF171.fastq
    S-1030_3,/path/to/S-1030_3.fastq
    """.stripIndent()
}

// Main workflow
workflow {
    // Show help message if requested and exit
    if (params.help) {
        helpMessage()
        exit 0
    }
    
    // Validate input parameters
    if (!params.samples && (!params.input_fastq || !params.name)) {
        error "Please specify either:\n" +
              "  - Both --input_fastq and --name for single sample, OR\n" +
              "  - --samples for CSV input"
    }
    
    if (params.samples && (params.input_fastq || params.name)) {
        error "Please specify either --samples OR --input_fastq/--name, not both"
    }
    
    // Run the main workflow
    ONT_FLYE()
}