include { CHOPPER } from '../modules/chopper'
include { NANOPLOT } from '../modules/nanoplot'
include { NANOPLOT as NANOPLOT_ORIGINAL } from '../modules/nanoplot'
include { PARSE_NANOSTATS } from '../modules/parse_nanostats'

workflow ONT_FLYE {
    
    // Validate input parameter
    if (!params.input_fastq) {
        error "Please specify an input FASTQ file with --input_fastq"
    }
    
    // Create input channel with path validation
    input_ch = Channel.fromPath(params.input_fastq, checkIfExists: true)
    
    // Create size ranges channel
    size_ranges_ch = Channel.fromList(params.size_ranges)
    
    // Combine input with each size range to create all combinations
    input_combinations = input_ch.combine(size_ranges_ch)
    
    // Run CHOPPER for each size range
    CHOPPER(input_combinations)
    
    // Run NANOPLOT on the original input file
    original_input_for_nanoplot = input_ch.map { file -> 
        tuple("original_input", file) 
    }
    NANOPLOT_ORIGINAL(original_input_for_nanoplot)
    
    // Run NANOPLOT on each filtered file
    NANOPLOT(CHOPPER.out.filtered_reads)
    
    // Collect all NanoPlot results and extract NanoStats.txt files
    all_nanoplot_results = NANOPLOT_ORIGINAL.out.nanoplot_results
        .mix(NANOPLOT.out.nanoplot_results)
        .collect()
    
    // Parse all NanoStats files and create summary table
    PARSE_NANOSTATS(all_nanoplot_results)
}