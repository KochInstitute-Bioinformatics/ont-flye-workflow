#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { ONT_FLYE } from './workflows/ont_flye'

workflow {
    ONT_FLYE()
}