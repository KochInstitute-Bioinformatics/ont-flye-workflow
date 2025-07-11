include { CHOPPER } from '../modules/chopper'
include { NANOPLOT } from '../modules/nanoplot'

workflow ONT_FLYE {
    
    // Validate input parameter
    if (!params.input_fastq) {
        error "Please specify an input FASTQ file with --input_fastq"
    }
    
    // Create input channel with path validation
    input_ch = Channel.fromPath(params.input_fastq, checkIfExists: true)
    
    // Create size ranges channel
    size_ranges_ch = Channel.fromList(params.size_ranges)
    
    // Run CHOPPER for each size range
    CHOPPER(input_ch, size_ranges_ch)
    
    // Run NANOPLOT on each filtered file
    NANOPLOT(CHOPPER.out.filtered_reads)
}