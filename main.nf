#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

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

include { ONT_FLYE } from './workflows/ont_flye'

workflow {
    // Show help and exit if requested
    if (params.help) {
        helpMessage()
        exit 0
    }
    
    ONT_FLYE()
}