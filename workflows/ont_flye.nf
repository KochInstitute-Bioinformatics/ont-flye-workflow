include { CHOPPER } from '../modules/chopper'
include { DOWNSAMPLE_FASTQ } from '../modules/downsample_fastq'
include { NANOPLOT } from '../modules/nanoplot'
include { NANOPLOT as NANOPLOT_ORIGINAL } from '../modules/nanoplot'

process PARSE_NANOSTATS {
    publishDir "${params.outdir}/summary", mode: 'copy'
    
    input:
    path nanoplot_results
    path parse_script
    
    output:
    path "nanostats_summary.csv", emit: summary_csv
    path "nanostats_summary.json", emit: summary_json
    path "versions.yml", emit: versions
    
    script:
    """
    python ${parse_script}
    """
    
    stub:
    """
    touch nanostats_summary.csv
    touch nanostats_summary.json
    touch versions.yml
    """
}

workflow ONT_FLYE {
    
    main:
    // Create input channel from samples.csv with validation
    if (params.samples) {
        // CSV input method with file validation
        input_ch = Channel
            .fromPath(params.samples, checkIfExists: true)
            .splitCsv(header: true)
            .map { row ->
                def fastq_file = file(row.fastq, checkIfExists: true)
                
                // Validate that it's a file, not a directory
                if (!fastq_file.isFile()) {
                    error "ERROR: ${row.fastq} is not a file! Please check your samples.csv - each 'fastq' entry must point to a FASTQ file, not a directory."
                }
                
                // Validate file extension
                def valid_extensions = ['.fastq', '.fq', '.fastq.gz', '.fq.gz']
                def has_valid_ext = valid_extensions.any { ext -> fastq_file.name.toLowerCase().endsWith(ext) }
                if (!has_valid_ext) {
                    error "ERROR: ${row.fastq} does not appear to be a FASTQ file. Valid extensions: ${valid_extensions.join(', ')}"
                }
                
                def transgene = row.containsKey('transgene') ? row.transgene : params.default_transgene
                [row.name, fastq_file, transgene]
            }
    } else {
        // Single sample input method (fallback)
        if (!params.input_fastq || !params.name) {
            error "For single sample mode, please specify both --input_fastq and --name"
        }
        def transgene = params.transgene ?: params.default_transgene
        input_ch = Channel.of([params.name, file(params.input_fastq, checkIfExists: true), transgene])
    }

    // ========================================
    // PHASE 1: FILTERING - Size and quality filtering using CHOPPER
    // ========================================
    
    // Create size ranges channel from custom_config
    size_ranges_ch = Channel.fromList(params.size_ranges)
    
    // Combine input with each size range to create all combinations
    filter_combinations = input_ch.combine(size_ranges_ch)
        .map { sample_name, fastq_file, transgene_name, size_range ->
            [sample_name, fastq_file, size_range]
        }
    
    // Run CHOPPER for each size range (filtering based on length and quality)
    CHOPPER(filter_combinations)
    
    // ========================================
    // PHASE 2: DOWNSAMPLING - Multiple replicates of downsampling
    // ========================================
    
    // Get the downsample script
    downsample_script = file("${projectDir}/bin/downsample_fastq.py", checkIfExists: true)
    
    // Create channels for downsample rates and replicates
    downsample_rates_ch = Channel.fromList(params.downsample_rates)
    replicate_numbers = Channel.from(1..params.replicates)
    
    // Create all combinations: filtered_reads × downsample_rates × replicates
    downsample_combinations = CHOPPER.out.filtered_reads
        .combine(downsample_rates_ch)
        .combine(replicate_numbers)
        .map { sample_name_size, filtered_fastq, fraction, replicate ->
            // Create unique sample name with downsample rate and replicate number
            def new_sample_name = "${sample_name_size}_ds${fraction}_rep${replicate}"
            [new_sample_name, filtered_fastq, fraction]
        }
    
    // Run DOWNSAMPLE_FASTQ on the size-filtered reads with multiple replicates
    DOWNSAMPLE_FASTQ(downsample_combinations, downsample_script)
    
    // ========================================
    // COLLECT ALL PROCESSED FASTQ FILES
    // ========================================
    
    // Combine all FASTQ files (filtered + downsampled) into one channel
    all_processed_fastq = CHOPPER.out.filtered_reads
        .mix(DOWNSAMPLE_FASTQ.out.downsampled_reads)
    
    // ========================================
    // QUALITY CONTROL AND ANALYSIS - Run on original input AND all processed files
    // ========================================
    
    // Run NANOPLOT on original input files (mark them as "original")
    original_input_for_nanoplot = input_ch.map { sample_name, fastq_file, transgene_name ->
        tuple("${sample_name}_original", fastq_file)
    }
    NANOPLOT_ORIGINAL(original_input_for_nanoplot)
    
    // Run NANOPLOT on all processed files
    NANOPLOT(all_processed_fastq)
    
    // Collect ALL NanoPlot results (original input + all processed files)
    all_nanoplot_results = NANOPLOT_ORIGINAL.out.nanoplot_results
        .mix(NANOPLOT.out.nanoplot_results)
        .collect()
    
    // Get the parse_nanostats script from bin directory
    parse_nanostats_script = file("${projectDir}/bin/parse_nanostats.py", checkIfExists: true)
    
    // Parse all NanoStats files and create summary table
    PARSE_NANOSTATS(all_nanoplot_results, parse_nanostats_script)
    
    emit:
    // Emit the key outputs for FASTQ generation phase
    original_input = input_ch
    selected_fastq = all_processed_fastq
    nanoplot_results = all_nanoplot_results
    nanostats_summary = PARSE_NANOSTATS.out.summary_csv
}