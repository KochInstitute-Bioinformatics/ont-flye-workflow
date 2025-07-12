include { CHOPPER } from '../modules/chopper'
include { NANOPLOT } from '../modules/nanoplot'
include { NANOPLOT as NANOPLOT_ORIGINAL } from '../modules/nanoplot'
include { PARSE_NANOSTATS } from '../modules/parse_nanostats'

workflow ONT_FLYE {
    
    // Create input channel based on input method
    if (params.samples) {
        // CSV input method
        input_ch = Channel
            .fromPath(params.samples, checkIfExists: true)
            .splitCsv(header: true)
            .map { row -> 
                [row.name, file(row.fastq, checkIfExists: true)]
            }
    } else {
        // Single sample input method
        if (!params.input_fastq || !params.name) {
            error "For single sample mode, please specify both --input_fastq and --name"
        }
        input_ch = Channel.of([params.name, file(params.input_fastq, checkIfExists: true)])
    }
    
    // Create size ranges channel
    size_ranges_ch = Channel.fromList(params.size_ranges)
    
    // Combine input with each size range to create all combinations
    // This creates tuples of [sample_name, fastq_file, size_range]
    input_combinations = input_ch.combine(size_ranges_ch)
    
    // Run CHOPPER for each size range
    CHOPPER(input_combinations)
    
    // Run NANOPLOT on the original input files
    original_input_for_nanoplot = input_ch.map { sample_name, file ->
        tuple("${sample_name}_original", file)
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