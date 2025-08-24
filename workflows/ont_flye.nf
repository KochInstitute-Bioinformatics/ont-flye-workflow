include { CHOPPER } from '../modules/chopper'
include { DOWNSAMPLE_FASTQ } from '../modules/downsample_fastq'
include { BOOTSTRAP_DOWNSAMPLE } from '../modules/bootstrap_downsample'
include { NANOPLOT } from '../modules/nanoplot'
include { NANOPLOT as NANOPLOT_ORIGINAL } from '../modules/nanoplot'
include { FLYE } from '../modules/flye'
include { GATHER_ASSEMBLY_STATS } from '../modules/gather_assembly_stats'
include { TRANSGENE_BLAST } from '../modules/transgene_blast'

process COPY_SELECTED_FASTQ {
    publishDir "${params.outdir}/selected_fastq", mode: 'copy'
    
    input:
    tuple val(sample_name), path(fastq_file)
    
    output:
    tuple val(sample_name), path("${sample_name}.fastq"), emit: selected_reads
    
    script:
    """
    cp ${fastq_file} ${sample_name}.fastq
    """
    
    stub:
    """
    touch ${sample_name}.fastq
    """
}

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
    // Create input channel based on input method
    if (params.samples) {
        // CSV input method
        input_ch = Channel
            .fromPath(params.samples, checkIfExists: true)
            .splitCsv(header: true)
            .map { row ->
                def transgene = row.containsKey('transgene') ? row.transgene : params.default_transgene
                [row.name, file(row.fastq, checkIfExists: true), transgene]
            }
    } else {
        // Single sample input method
        if (!params.input_fastq || !params.name) {
            error "For single sample mode, please specify both --input_fastq and --name"
        }
        def transgene = params.transgene ?: params.default_transgene
        input_ch = Channel.of([params.name, file(params.input_fastq, checkIfExists: true), transgene])
    }

    // ========================================
    // PHASE 1: SCAN MODE - Size-based filtering using CHOPPER
    // ========================================
    
    // Create size ranges channel from custom_config
    size_ranges_ch = Channel.fromList(params.size_ranges)
    
    // Combine input with each size range to create all combinations
    scan_combinations = input_ch.combine(size_ranges_ch)
        .map { sample_name, fastq_file, transgene_name, size_range ->
            [sample_name, fastq_file, size_range]
        }
    
    // Run CHOPPER for each size range (scan mode functionality)
    CHOPPER(scan_combinations)
    
    // Copy the size-filtered reads to selected_fastq directory
    COPY_SELECTED_FASTQ(CHOPPER.out.filtered_reads)
    
    // ========================================
    // PHASE 2: DOWNSAMPLE MODE - Use size-filtered reads for downsampling
    // ========================================
    
    // Get the downsample script
    downsample_script = file("${projectDir}/bin/downsample_fastq.py", checkIfExists: true)
    
    // Create fractions channel
    fractions_ch = Channel.fromList(params.downsample_rates)
    
    // Use the filtered reads from CHOPPER and combine with downsample rates
    // This creates the "selected_fastq" that will be used for downsampling
    downsample_combinations = CHOPPER.out.filtered_reads
        .combine(fractions_ch)
        .map { sample_name_size, filtered_fastq, fraction ->
            [sample_name_size, filtered_fastq, fraction]
        }
    
    // Run DOWNSAMPLE_FASTQ on the size-filtered reads
    DOWNSAMPLE_FASTQ(downsample_combinations, downsample_script)
    
    // ========================================
    // PHASE 3: BOOTSTRAP MODE - Use downsampled reads for bootstrap
    // ========================================
    
    // Generate replicate numbers from 1 to params.replicates
    replicate_numbers = Channel.from(1..params.replicates)
    
    // Use the downsampled reads for bootstrap analysis
    // Select one of the downsampled outputs (e.g., the highest fraction) for bootstrap
    selected_downsampled = DOWNSAMPLE_FASTQ.out.downsampled_reads
        .filter { sample_name, fastq_file ->
            // Select the highest downsample rate for bootstrap
            sample_name.contains("ds${params.downsample_rates.max()}")
        }
        .map { sample_name, fastq_file ->
            // Remove the downsample suffix for bootstrap naming
            def base_name = sample_name.replaceAll(/_ds[\d\.]+$/, '')
            [base_name, fastq_file]
        }
    
    // Create bootstrap combinations
    bootstrap_combinations = selected_downsampled
        .combine(Channel.of(params.fraction ?: 0.75))
        .combine(replicate_numbers)
    
    // Run bootstrap downsampling
    BOOTSTRAP_DOWNSAMPLE(bootstrap_combinations, downsample_script)
    
    // ========================================
    // QUALITY CONTROL AND ANALYSIS
    // ========================================
    
    // Run NANOPLOT on original input files
    original_input_for_nanoplot = input_ch.map { sample_name, fastq_file, transgene_name ->
        tuple("${sample_name}_original", fastq_file)
    }
    NANOPLOT_ORIGINAL(original_input_for_nanoplot)
    
    // Run NANOPLOT on filtered reads (from CHOPPER)
    NANOPLOT(CHOPPER.out.filtered_reads)
    
    // Run NANOPLOT on downsampled reads
    nanoplot_downsampled = NANOPLOT(DOWNSAMPLE_FASTQ.out.downsampled_reads)
    
    // Run NANOPLOT on bootstrap reads
    nanoplot_bootstrap = NANOPLOT(BOOTSTRAP_DOWNSAMPLE.out.bootstrap_reads)
    
    // Collect all NanoPlot results
    all_nanoplot_results = NANOPLOT_ORIGINAL.out.nanoplot_results
        .mix(NANOPLOT.out.nanoplot_results)
        .mix(nanoplot_downsampled.out.nanoplot_results)
        .mix(nanoplot_bootstrap.out.nanoplot_results)
        .collect()
    
    // Get the parse_nanostats script from bin directory
    parse_nanostats_script = file("${projectDir}/bin/parse_nanostats.py", checkIfExists: true)
    
    // Parse all NanoStats files and create summary table
    PARSE_NANOSTATS(all_nanoplot_results, parse_nanostats_script)
    
    // ========================================
    // ASSEMBLY AND ANALYSIS
    // ========================================
    
    // Combine all reads for assembly: filtered, downsampled, and bootstrap
    all_reads_for_assembly = CHOPPER.out.filtered_reads
        .mix(DOWNSAMPLE_FASTQ.out.downsampled_reads)
        .mix(BOOTSTRAP_DOWNSAMPLE.out.bootstrap_reads)
    
    // Run FLYE assembly on all read sets
    FLYE(all_reads_for_assembly)
    
    // ========================================
    // TRANSGENE BLAST ANALYSIS
    // ========================================
    
    // Get transgene files from transgeneDir if specified
    if (params.transgeneDir) {
        transgene_files = Channel.fromPath("${params.transgeneDir}/*.{fa,fasta,fna}")
            .ifEmpty { error "No FASTA files found in transgeneDir: ${params.transgeneDir}" }
    } else {
        // Fallback to default transgene (you may need to adjust this path)
        transgene_files = Channel.fromPath("${projectDir}/transgenes/A-vector_herceptin_pEY345.fa")
    }
    
    // Prepare BLAST input channel - combine assemblies with each transgene file
    blast_input_ch = FLYE.out.assembly_fasta
        .combine(transgene_files)
        .map { sample_name, assembly_fasta, transgene_fasta ->
            def transgene_name = transgene_fasta.baseName
            [sample_name, assembly_fasta, transgene_name, transgene_fasta]
        }
    
    // Run transgene BLAST analysis
    TRANSGENE_BLAST(blast_input_ch)
    
    // ========================================
    // GATHER ASSEMBLY STATISTICS
    // ========================================
    
    // Collect assembly statistics
    assembly_info_collected = FLYE.out.assembly_info
        .map { sample_name, assembly_file -> assembly_file }
        .collect()
    
    flye_log_collected = FLYE.out.flye_log
        .map { sample_name, log_file -> log_file }
        .collect()
    
    sample_names_collected = FLYE.out.assembly_info
        .map { sample_name, assembly_file -> sample_name }
        .collect()
    
    GATHER_ASSEMBLY_STATS(assembly_info_collected, flye_log_collected, sample_names_collected)
    
    emit:
    // Emit the key outputs
    selected_fastq = COPY_SELECTED_FASTQ.out.selected_reads
    downsampled_reads = DOWNSAMPLE_FASTQ.out.downsampled_reads
    bootstrap_reads = BOOTSTRAP_DOWNSAMPLE.out.bootstrap_reads
    nanoplot_results = all_nanoplot_results
    nanostats_summary = PARSE_NANOSTATS.out.summary_csv
    assembly_stats = GATHER_ASSEMBLY_STATS.out.assembly_stats
    blast_results = TRANSGENE_BLAST.out.blast_results
}