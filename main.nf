#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

// Include the main workflow
include { ONT_FLYE } from './workflows/ont_flye'

// Parameters with defaults
params.help = false
params.input_fastq = null
params.outdir = "results"

// Help message function
def helpMessage() {
    log.info"""
    Usage:
    nextflow run main.nf --input_fastq <path_to_fastq>
    
    Required arguments:
      --input_fastq    Path to input FASTQ file
    
    Optional arguments:
      --outdir         Output directory (default: results)
      --min_quality    Minimum quality score (default: 10)
      --help           Show this help message
    """.stripIndent()
}

// Main workflow
workflow {
    // Show help message if requested and exit
    if (params.help) {
        helpMessage()
        exit 0
    }
    
    // Run the main workflow
    ONT_FLYE()
}