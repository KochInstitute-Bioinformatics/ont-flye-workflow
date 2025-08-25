include { CHOPPER } from '../modules/chopper'
include { DOWNSAMPLE_FASTQ } from '../modules/downsample_fastq'
include { FLYE_PREFLIGHT } from '../modules/flye_preflight'
include { PARSE_PREFLIGHT_RESULTS } from '../modules/parse_preflight_results'
include { FILTER_ASSEMBLY_CANDIDATES } from '../modules/filter_assembly_candidates'
include { FLYE } from '../modules/flye'
include { TRANSGENE_BLAST } from '../modules/transgene_blast'  // Add this line
include { PARSE_NANOSTATS } from '../modules/parse_nanostats'
include { NANOPLOT } from '../modules/nanoplot'
include { NANOPLOT as NANOPLOT_ORIGINAL } from '../modules/nanoplot'

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
        // PHASE 3: FLYE PREFLIGHT - Run preflight on all selected FASTQ files
        // ========================================
        // Run FLYE_PREFLIGHT on all processed FASTQ files
        FLYE_PREFLIGHT(all_processed_fastq)

        // Collect all preflight logs
        all_preflight_logs = FLYE_PREFLIGHT.out.preflight_logs
            .map { sample_name, log_file -> log_file }
            .collect()

        // Get the parse_preflight script from bin directory
        parse_preflight_script = file("${projectDir}/bin/parse_preflight_flyelog.py", checkIfExists: true)

        // Parse all preflight logs and create summary table
        PARSE_PREFLIGHT_RESULTS(all_preflight_logs, parse_preflight_script)

        // ========================================
        // PHASE 4: COVERAGE-BASED ASSEMBLY DECISION
        // ========================================
        // Filter assembly candidates based on coverage criteria
        FILTER_ASSEMBLY_CANDIDATES(
            PARSE_PREFLIGHT_RESULTS.out.preflight_csv
        )

        // Create channel of FASTQ files that should be assembled
        // Read the candidates CSV and match with FASTQ files
        assembly_candidates = FILTER_ASSEMBLY_CANDIDATES.out.candidates_csv
            .splitCsv(header: true)
            .map { row -> row.sample_name }
            .combine(all_processed_fastq)
            .filter { candidate_name, sample_name, fastq_file ->
                candidate_name == sample_name
            }
            .map { candidate_name, sample_name, fastq_file ->
                [sample_name, fastq_file]
            }

        // ========================================
        // PHASE 5: FLYE ASSEMBLY - Only on selected candidates
        // ========================================
        // Run FLYE assembly only on candidates that passed coverage criteria
        FLYE(assembly_candidates)

        // ========================================
        // PHASE 6: TRANSGENE BLAST ANALYSIS - Run on assembled genomes
        // ========================================

        // Create transgene channel from CSV file
        transgene_ch = Channel
            .fromPath(params.transgene_library, checkIfExists: true)
            .splitCsv(header: true)
            .map { row -> 
                [row.transgene_name, file("${params.transgeneDir}/${row.fasta_file}", checkIfExists: true)]
            }

        // Combine assemblies with transgene information from input
        assembly_with_transgene = FLYE.out.assembly_fasta
            .combine(input_ch.map { sample_name, fastq_file, transgene_name -> [sample_name, transgene_name] })
            .filter { assembly_sample, assembly_fasta, input_sample, transgene_name ->
                // Match assembly samples with their original transgene assignments
                assembly_sample.startsWith(input_sample.split('_')[0]) // Handle sample name variations
            }
            .map { assembly_sample, assembly_fasta, input_sample, transgene_name ->
                [assembly_sample, assembly_fasta, transgene_name]
            }

        // Join with transgene files
        blast_input = assembly_with_transgene
            .combine(transgene_ch)
            .filter { assembly_sample, assembly_fasta, transgene_name, transgene_file_name, transgene_file ->
                transgene_name == transgene_file_name
            }
            .map { assembly_sample, assembly_fasta, transgene_name, transgene_file_name, transgene_file ->
                [assembly_sample, assembly_fasta, transgene_name, transgene_file]
            }

        // Run TRANSGENE_BLAST
        TRANSGENE_BLAST(blast_input)

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
        preflight_logs = FLYE_PREFLIGHT.out.preflight_logs
        preflight_summary = PARSE_PREFLIGHT_RESULTS.out.preflight_csv
        assembly_candidates = FILTER_ASSEMBLY_CANDIDATES.out.candidates_csv
        assembly_filtered = FILTER_ASSEMBLY_CANDIDATES.out.filtered_csv
        assemblies = FLYE.out.assembly_fasta
        blast_results = TRANSGENE_BLAST.out.blast_results  // Add this line
        nanoplot_results = all_nanoplot_results
        nanostats_summary = PARSE_NANOSTATS.out.summary_csv
}