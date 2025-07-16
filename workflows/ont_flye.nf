include { SCAN_MODE } from './scan_mode'
include { DOWNSAMPLE_MODE } from './downsample_mode'

workflow ONT_FLYE {
    // Create input channel based on input method
    if (params.samples) {
        // CSV input method - now supports transgene and mode columns
        input_ch = Channel
            .fromPath(params.samples, checkIfExists: true)
            .splitCsv(header: true)
            .map { row ->
                // Check if columns exist, if not use defaults
                def transgene = row.containsKey('transgene') ? row.transgene : params.default_transgene
                def mode = row.containsKey('mode') ? row.mode : params.default_mode
                [row.name, file(row.fastq, checkIfExists: true), transgene, mode]
            }
    } else {
        // Single sample input method
        if (!params.input_fastq || !params.name) {
            error "For single sample mode, please specify both --input_fastq and --name"
        }
        def transgene = params.transgene ?: params.default_transgene
        def mode = params.mode ?: params.default_mode
        input_ch = Channel.of([params.name, file(params.input_fastq, checkIfExists: true), transgene, mode])
    }

    // Split input by mode
    scan_samples = input_ch
        .filter { sample_name, fastq_file, transgene_name, mode -> mode == 'scan' }
        .map { sample_name, fastq_file, transgene_name, _mode -> [sample_name, fastq_file, transgene_name] }
    
    downsample_samples = input_ch
        .filter { sample_name, fastq_file, transgene_name, mode -> mode == 'downsample' }
        .map { sample_name, fastq_file, transgene_name, _mode -> [sample_name, fastq_file, transgene_name] }

    // Run workflows - Nextflow will automatically handle empty channels
    SCAN_MODE(scan_samples)
    DOWNSAMPLE_MODE(downsample_samples)
}