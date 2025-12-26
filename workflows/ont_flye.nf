include { CHOPPER } from '../modules/local/chopper'
include { DOWNSAMPLE_FASTQ } from '../modules/local/downsample_fastq'
include { FLYE_PREFLIGHT } from '../modules/local/flye_preflight'
include { PARSE_PREFLIGHT_RESULTS } from '../modules/local/parse_preflight_results'
include { FILTER_ASSEMBLY_CANDIDATES } from '../modules/local/filter_assembly_candidates'
include { FLYE } from '../modules/local/flye'
include { TRANSGENE_BLAST } from '../modules/local/transgene_blast'
include { PARSE_NANOSTATS } from '../modules/local/parse_nanostats'
include { NANOPLOT } from '../modules/local/nanoplot'
include { NANOPLOT as NANOPLOT_ORIGINAL } from '../modules/local/nanoplot'
include { PARSE_TRANSGENE_BLAST } from '../modules/local/parse_transgene_blast'
include { GATHER_ASSEMBLY_STATS } from '../modules/local/gather_assembly_stats'
include { SIMPLE_RESULTS_SUMMARY } from '../modules/local/simple_results_summary'
include { ALIGN_ASSEMBLY_TO_GENOME } from '../modules/local/assembly_evaluation'
include { REPAIR_ASSEMBLY } from '../modules/local/assembly_evaluation'
include { ALIGN_ANNOTATED_ASSEMBLY } from '../modules/local/assembly_evaluation'
include { CALCULATE_CHROMOSOME_COVERAGE } from '../modules/local/assembly_evaluation'
include { CONSOLIDATE_COVERAGE_SUMMARY } from '../modules/local/assembly_evaluation'
include { FINALIZE_ASSEMBLY } from '../modules/local/assembly_evaluation'
include { ALIGN_FINAL_ASSEMBLY } from '../modules/local/assembly_evaluation'
include { MAP_READS_TO_ASSEMBLY } from '../modules/local/assembly_evaluation'
include { BLAST_TRANSGENE_TO_ASSEMBLY } from '../modules/local/assembly_evaluation'
include { CONVERT_BLAST_TO_BED } from '../modules/local/assembly_evaluation'
include { MAP_TRANSCRIPTS_TO_ASSEMBLY } from '../modules/local/assembly_evaluation'
include { CONSOLIDATE_IGV_DATA } from '../modules/local/assembly_evaluation'


// ========================================
// MAIN WORKFLOW
// ========================================

workflow ONT_FLYE {

main:

// Create input channel from samples.csv with validation and per-sample parameters
if (params.samples) {
    // CSV input method with file validation and parameter parsing
    input_ch = channel
        .fromPath(params.samples, checkIfExists: true)
        .splitCsv(header: true, sep: ';')  // Use semicolon separator
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
            
            // Parse size_ranges with clean JSON
            def size_ranges
            if (row.containsKey('size_ranges') && row.size_ranges && row.size_ranges.trim() != '') {
                try {
                    log.info "Parsing size_ranges for ${row.name}: ${row.size_ranges}"
                    def parsed_ranges = new groovy.json.JsonSlurper().parseText(row.size_ranges)
                    
                    // Convert to the expected format with max: null
                    size_ranges = parsed_ranges.collect { range ->
                        [min: range.min as Integer, max: null, name: range.name as String]
                    }
                    log.info "Successfully parsed size_ranges for ${row.name}: ${size_ranges}"
                } catch (Exception e) {
                    log.warn "Failed to parse size_ranges for ${row.name}: ${e.message}. Using defaults."
                    log.warn "Raw size_ranges string: '${row.size_ranges}'"
                    size_ranges = [
                        [min: 40000, max: null, name: "40k_Plus"],
                        [min: 50000, max: null, name: "50k_Plus"]
                    ]
                }
            } else {
                // Default size ranges
                size_ranges = [
                    [min: 40000, max: null, name: "40k_Plus"],
                    [min: 50000, max: null, name: "50k_Plus"]
                ]
            }
            
            // Parse downsample_rates with clean JSON
            def downsample_rates
            if (row.containsKey('downsample_rates') && row.downsample_rates && row.downsample_rates.trim() != '') {
                try {
                    log.info "Parsing downsample_rates for ${row.name}: ${row.downsample_rates}"
                    def parsed_rates = new groovy.json.JsonSlurper().parseText(row.downsample_rates)
                    
                    // Convert to list of doubles
                    downsample_rates = parsed_rates.collect { rate -> rate as Double }
                    log.info "Successfully parsed downsample_rates for ${row.name}: ${downsample_rates}"
                } catch (Exception e) {
                    log.warn "Failed to parse downsample_rates for ${row.name}: ${e.message}. Using defaults."
                    log.warn "Raw downsample_rates string: '${row.downsample_rates}'"
                    downsample_rates = [0.25, 0.5]
                }
            } else {
                // Default downsample rates
                downsample_rates = [0.25, 0.5]
            }
            
            [row.name, fastq_file, transgene, size_ranges, downsample_rates]
        }
} else {
    // Single sample input method (fallback) - use global defaults
    if (!params.input_fastq || !params.name) {
        error "For single sample mode, please specify both --input_fastq and --name"
    }
    def transgene = params.transgene ?: params.default_transgene
    def default_size_ranges = [
        [min: 40000, max: null, name: "40k_Plus"],
        [min: 50000, max: null, name: "50k_Plus"]
    ]
    def default_downsample_rates = [0.25, 0.5]
    input_ch = channel.of([params.name, file(params.input_fastq, checkIfExists: true), transgene, default_size_ranges, default_downsample_rates])
}

// ========================================
// PHASE 1: FILTERING - Size and quality filtering using CHOPPER
// ========================================

// Create combinations of samples with their specific size ranges
filter_combinations = input_ch
    .flatMap { sample_name, fastq_file, transgene_name, size_ranges, downsample_rates ->
        // Create a combination for each size range for this sample
        size_ranges.collect { size_range ->
            [sample_name, fastq_file, size_range, transgene_name, downsample_rates]
        }
    }
    .map { sample_name, fastq_file, size_range, _transgene_name, _downsample_rates ->
        [sample_name, fastq_file, size_range]
    }

// Run CHOPPER for each size range (filtering based on length and quality)
CHOPPER(filter_combinations)

// ========================================
// PHASE 2: DOWNSAMPLING - Multiple replicates of downsampling
// ========================================

// Get the downsample script
downsample_script = file("${projectDir}/bin/downsample_fastq.py", checkIfExists: true)

// Create replicate numbers channel
replicate_numbers = channel.from(1..params.replicates)

// Create downsample combinations using per-sample downsample rates
// First, create a lookup map for sample-specific downsample rates
sample_downsample_map = input_ch
    .map { sample_name, _fastq_file, _transgene_name, _size_ranges, downsample_rates ->
        [sample_name, downsample_rates]
    }

// Create downsample combinations
downsample_combinations = CHOPPER.out.filtered_reads
    .combine(sample_downsample_map)
    .filter { filtered_sample_name, _filtered_fastq, original_sample_name, _downsample_rates ->
        // Match filtered samples back to their original sample for downsample rates
        filtered_sample_name.startsWith(original_sample_name)
    }
    .flatMap { filtered_sample_name, filtered_fastq, _original_sample_name, downsample_rates ->
        // Create combinations with each downsample rate for this sample
        downsample_rates.collect { rate ->
            [filtered_sample_name, filtered_fastq, rate]
        }
    }
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
    .map { _sample_name, log_file -> log_file }
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
    PARSE_PREFLIGHT_RESULTS.out.preflight_csv,
    params.min_assembly_depth,
    params.max_assembly_depth
)

// Create channel of FASTQ files that should be assembled
// Read the candidates CSV and match with FASTQ files
assembly_candidates = FILTER_ASSEMBLY_CANDIDATES.out.candidates_csv
    .splitCsv(header: true)
    .map { row -> row.sample_name }
    .combine(all_processed_fastq)
    .filter { candidate_name, sample_name, _fastq_file ->
        candidate_name == sample_name
    }
    .map { _candidate_name, sample_name, fastq_file ->
        [sample_name, fastq_file]
    }

// ========================================
// PHASE 5: FLYE ASSEMBLY - Only on selected candidates
// ========================================

// Run FLYE assembly only on candidates that passed coverage criteria
FLYE(assembly_candidates)

// Store assembly with its source FASTQ for downstream processes
assembly_with_source_fastq = FLYE.out.assembly_fasta
    .join(assembly_candidates, by: 0)  // Join by sample_name to get the FASTQ that was used

// ========================================
// PHASE 6: TRANSGENE BLAST ANALYSIS - Run on assembled genomes
// ========================================

// Create transgene channel from CSV file
transgene_ch = channel
    .fromPath(params.transgene_library, checkIfExists: true)
    .splitCsv(header: true)
    .map { row ->
        [row.transgene_name, file("${params.transgene_dir}/${row.fasta_file}", checkIfExists: true)]
    }

// Combine assemblies with transgene information from input
assembly_with_transgene = FLYE.out.assembly_fasta
    .combine(input_ch.map { sample_name, _fastq_file, transgene_name, _size_ranges, _downsample_rates -> [sample_name, transgene_name] })
    .filter { assembly_sample, _assembly_fasta, input_sample, _transgene_name ->
        // Match assembly samples with their original transgene assignments
        assembly_sample.startsWith(input_sample.split('_')[0]) // Handle sample name variations
    }
    .map { assembly_sample, assembly_fasta, _input_sample, transgene_name ->
        [assembly_sample, assembly_fasta, transgene_name]
    }

// Join with transgene files
blast_input = assembly_with_transgene
    .combine(transgene_ch)
    .filter { _assembly_sample, _assembly_fasta, transgene_name, transgene_file_name, _transgene_file ->
        transgene_name == transgene_file_name
    }
    .map { assembly_sample, assembly_fasta, transgene_name, _transgene_file_name, transgene_file ->
        [assembly_sample, assembly_fasta, transgene_name, transgene_file]
    }

// Run TRANSGENE_BLAST
TRANSGENE_BLAST(blast_input)

// ========================================
// PHASE 7: PARSE TRANSGENE BLAST RESULTS
// ========================================

// Collect all BLAST result files with unique names to avoid collisions
all_blast_results = TRANSGENE_BLAST.out.blast_results
    .map { _sample_name, blast_file -> blast_file }
    .collect()

// Get the parse script from bin directory
parse_transgene_script = file("${projectDir}/bin/parse_transgene_blast.py", checkIfExists: true)

// Parse all BLAST results
PARSE_TRANSGENE_BLAST(
    all_blast_results,
    file(params.transgene_library),
    parse_transgene_script
)

// ========================================
// PHASE 8: GATHER ASSEMBLY STATISTICS
// ========================================

// Collect assembly info files and flye logs from successful assemblies
assembly_info_files = FLYE.out.assembly_info
    .map { _sample_name, info_file -> info_file }
    .collect()

flye_log_files = FLYE.out.flye_log
    .map { _sample_name, log_file -> log_file }
    .collect()

// Extract sample names from successful assemblies
assembly_sample_names = FLYE.out.assembly_fasta
    .map { sample_name, _fasta_file -> sample_name }
    .collect()

// Run GATHER_ASSEMBLY_STATS
GATHER_ASSEMBLY_STATS(
    assembly_info_files,
    flye_log_files,
    assembly_sample_names,
    file("${projectDir}/bin/gather_assembly_stats.py")  // Pass script as input
)

// ========================================
// QUALITY CONTROL AND ANALYSIS - Run on original input AND all processed files
// ========================================

// Run NANOPLOT on original input files (mark them as "original")
original_input_for_nanoplot = input_ch.map { sample_name, fastq_file, _transgene_name, _size_ranges, _downsample_rates ->
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

// results summary (simple version)
simple_summary_script = file("${projectDir}/bin/simple_summary.py", checkIfExists: true)

// Update the SIMPLE_RESULTS_SUMMARY call
SIMPLE_RESULTS_SUMMARY(
    PARSE_NANOSTATS.out.summary_json,
    PARSE_PREFLIGHT_RESULTS.out.preflight_json,
    GATHER_ASSEMBLY_STATS.out.assembly_stats,
    PARSE_TRANSGENE_BLAST.out.json_results,
    simple_summary_script
)

// ========================================
// PHASE 9: ASSEMBLY EVALUATION (OPTIONAL)
// ========================================

if (params.run_assembly_evaluation && params.reference_genome) {
    
    log.info "Assembly evaluation enabled - running complete refinement pipeline"
    
    // Prepare reference genome
    reference_genome = file(params.reference_genome, checkIfExists: true)
    
    // Prepare transcripts FASTA (optional)
    transcripts_fasta = params.transcripts_fasta ? 
        file(params.transcripts_fasta, checkIfExists: true) : 
        null
    
    // Prepare assembly input with reference genome
    assembly_with_ref = FLYE.out.assembly_fasta
        .map { sample_name, assembly_fasta -> 
            tuple(sample_name, assembly_fasta, reference_genome)
        }
    
    // STEP 1: Align assembly to reference genome
    ALIGN_ASSEMBLY_TO_GENOME(
        assembly_with_ref
    )
    
    // STEP 2: Repair and annotate assembly based on alignment
    REPAIR_ASSEMBLY(
        ALIGN_ASSEMBLY_TO_GENOME.out.alignment
    )
    
    // STEP 3: Align annotated assembly to reference
    annotated_with_ref = REPAIR_ASSEMBLY.out.annotated_assembly
        .map { sample_name, annotated_fasta ->
            tuple(sample_name, annotated_fasta, reference_genome)
        }
    
    ALIGN_ANNOTATED_ASSEMBLY(
        annotated_with_ref
    )
    
    // STEP 3a: Calculate chromosome coverage (if reference_chromosomes is provided)
    if (params.reference_chromosomes) {
        reference_chromosomes = file(params.reference_chromosomes, checkIfExists: true)
        
        annotated_with_chr_ref = REPAIR_ASSEMBLY.out.annotated_assembly
            .map { sample_name, annotated_fasta ->
                tuple(sample_name, annotated_fasta, reference_chromosomes)
            }
        
        CALCULATE_CHROMOSOME_COVERAGE(
            annotated_with_chr_ref
        )
        
        // Collect all coverage files and create summary
        all_coverage_files = CALCULATE_CHROMOSOME_COVERAGE.out.coverage
            .map { _sample_name, coverage_file -> coverage_file }
            .collect()
        
        CONSOLIDATE_COVERAGE_SUMMARY(
            all_coverage_files
        )
    }
    
    // STEP 4: Finalize assembly (concatenate chromosomal contigs)
    FINALIZE_ASSEMBLY(
        ALIGN_ANNOTATED_ASSEMBLY.out.alignment
    )
    
    // STEP 5: Align final assembly to reference
    final_with_ref = FINALIZE_ASSEMBLY.out.final_assembly
        .map { sample_name, final_fasta ->
            tuple(sample_name, final_fasta, reference_genome)
        }
    
    ALIGN_FINAL_ASSEMBLY(
        final_with_ref
    )
    
    // STEP 6: Map ONT reads to final assembly
    // Use the tracked FASTQ file that was actually used for assembly
    final_with_reads = FINALIZE_ASSEMBLY.out.final_assembly
        .join(assembly_with_source_fastq, by: 0)
        .map { sample_name, final_fasta, _assembly_fasta, source_fastq ->
            tuple(sample_name, final_fasta, source_fastq)
        }
    
    MAP_READS_TO_ASSEMBLY(
        final_with_reads
    )
    
    // STEP 7: BLAST transgene against final assembly
    // Create a map of base sample names to transgene info from input CSV
    transgene_map = input_ch
        .map { sample_name, _fastq_file, transgene_name, _size_ranges, _downsample_rates ->
            tuple(sample_name, transgene_name)
        }
        .unique()
    
    // Join final assemblies with transgene info
    // Extract base sample name to match with original input
    final_with_transgene = FINALIZE_ASSEMBLY.out.final_assembly
        .map { sample_name, final_fasta ->
            // Extract base sample name (e.g., S-1077-1 from S-1077-1_70k_Plus_ds0.5_rep1)
            def base_sample = sample_name.replaceAll(/_\d+k_Plus.*/, '')
            tuple(base_sample, sample_name, final_fasta)
        }
        .combine(transgene_map, by: 0)  // Join by base_sample
        .map { _base_sample, full_sample_name, final_fasta, transgene_name ->
            def transgene_file = file("${params.transgene_dir}/${transgene_name}.fa", checkIfExists: true)
            tuple(full_sample_name, final_fasta, transgene_name, transgene_file)
        }
    
    BLAST_TRANSGENE_TO_ASSEMBLY(
        final_with_transgene
    )
    
    // STEP 8: Convert BLAST results to BED format
    CONVERT_BLAST_TO_BED(
        BLAST_TRANSGENE_TO_ASSEMBLY.out.blast_results
            .map { sample_name, blast_file, _final_assembly ->
                tuple(sample_name, blast_file)
            }
    )
    
    // STEP 9: Map transcripts to final assembly (if provided)
    if (params.transcripts_fasta) {
        MAP_TRANSCRIPTS_TO_ASSEMBLY(
            FINALIZE_ASSEMBLY.out.final_assembly,
            transcripts_fasta
        )
    }
    
    // STEP 10: Consolidate all IGV data into strain-specific directories
    // Combine all the data needed for IGV visualization
    
    // Start with final assembly and BAM files
    igv_data_base = FINALIZE_ASSEMBLY.out.final_assembly
        .join(MAP_READS_TO_ASSEMBLY.out.mapped_reads, by: 0)
        .map { sample_name, final_asm, _assembly_copy, bam_file ->
            tuple(sample_name, final_asm, bam_file)
        }
        .join(MAP_READS_TO_ASSEMBLY.out.bam_index.map { sample_name, bai -> tuple(sample_name, bai) }, by: 0)
    
    // Add transgene BED files
    igv_with_transgene = igv_data_base
        .join(CONVERT_BLAST_TO_BED.out.bed_file, by: 0)
    
    // Add transcript BED files (if available)
    if (params.transcripts_fasta) {
        igv_complete = igv_with_transgene
            .join(MAP_TRANSCRIPTS_TO_ASSEMBLY.out.transcript_bed, by: 0)
            .map { sample_name, final_asm, bam, bam_idx, trans_bed, transcript_bed ->
                def base_strain = sample_name.replaceAll(/_\d+k_Plus.*/, '')
                tuple(base_strain, sample_name, final_asm, bam, bam_idx, trans_bed, transcript_bed)
            }
    } else {
        // Create a dummy transcript file for optional input
        def dummy_transcript = file('transcript_optional.bed')
        igv_complete = igv_with_transgene
            .map { sample_name, final_asm, bam, bam_idx, trans_bed ->
                def base_strain = sample_name.replaceAll(/_\d+k_Plus.*/, '')
                tuple(base_strain, sample_name, final_asm, bam, bam_idx, trans_bed, dummy_transcript)
            }
    }
    
    CONSOLIDATE_IGV_DATA(igv_complete)
}

emit:
    // Emit the key outputs for FASTQ generation phase
    original_input = input_ch
    selected_fastq = all_processed_fastq
    preflight_logs = FLYE_PREFLIGHT.out.preflight_logs
    preflight_summary = PARSE_PREFLIGHT_RESULTS.out.preflight_csv
    assembly_candidates = FILTER_ASSEMBLY_CANDIDATES.out.candidates_csv
    assembly_filtered = FILTER_ASSEMBLY_CANDIDATES.out.filtered_csv
    assemblies = FLYE.out.assembly_fasta
    blast_results = TRANSGENE_BLAST.out.blast_results
    transgene_summary_json = PARSE_TRANSGENE_BLAST.out.json_results
    transgene_summary_csv = PARSE_TRANSGENE_BLAST.out.csv_results
    nanoplot_results = all_nanoplot_results
    nanostats_summary = PARSE_NANOSTATS.out.summary_csv
    assembly_summary = GATHER_ASSEMBLY_STATS.out.assembly_stats
    
    // Assembly evaluation outputs (only if enabled)
    // Complete refinement pipeline outputs
    initial_alignments = params.run_assembly_evaluation && params.reference_genome ? 
        ALIGN_ASSEMBLY_TO_GENOME.out.alignment : channel.empty()
    annotated_assemblies = params.run_assembly_evaluation && params.reference_genome ? 
        REPAIR_ASSEMBLY.out.annotated_assembly : channel.empty()
    annotated_alignments = params.run_assembly_evaluation && params.reference_genome ? 
        ALIGN_ANNOTATED_ASSEMBLY.out.alignment : channel.empty()
    final_assemblies = params.run_assembly_evaluation && params.reference_genome ? 
        FINALIZE_ASSEMBLY.out.final_assembly : channel.empty()
    final_alignments = params.run_assembly_evaluation && params.reference_genome ? 
        ALIGN_FINAL_ASSEMBLY.out.alignment : channel.empty()
    read_mappings = params.run_assembly_evaluation && params.reference_genome ? 
        MAP_READS_TO_ASSEMBLY.out.mapped_reads : channel.empty()
    transgene_blast_beds = params.run_assembly_evaluation && params.reference_genome ? 
        CONVERT_BLAST_TO_BED.out.bed_file : channel.empty()
    transcript_alignments = params.run_assembly_evaluation && params.reference_genome && params.transcripts_fasta ? 
        MAP_TRANSCRIPTS_TO_ASSEMBLY.out.transcript_bed : channel.empty()
    coverage_summary = params.run_assembly_evaluation && params.reference_genome && params.reference_chromosomes ? 
        CONSOLIDATE_COVERAGE_SUMMARY.out.coverage_summary : channel.empty()
}