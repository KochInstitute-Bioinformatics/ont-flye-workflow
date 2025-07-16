include { CHOPPER } from '../modules/chopper'
include { NANOPLOT } from '../modules/nanoplot'
include { NANOPLOT as NANOPLOT_ORIGINAL } from '../modules/nanoplot'
include { PARSE_NANOSTATS } from '../modules/parse_nanostats'
include { FLYE } from '../modules/flye'
include { GATHER_ASSEMBLY_STATS } from '../modules/gather_assembly_stats'
include { TRANSGENE_BLAST } from '../modules/transgene_blast'

workflow ONT_FLYE {
    // Create input channel based on input method
    if (params.samples) {
        // CSV input method - now supports transgene column
        input_ch = Channel
            .fromPath(params.samples, checkIfExists: true)
            .splitCsv(header: true)
            .map { row ->
                // Check if transgene column exists, if not use default
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

    // Create size ranges channel
    size_ranges_ch = Channel.fromList(params.size_ranges)

    // Combine input with each size range to create all combinations
    input_combinations = input_ch.combine(size_ranges_ch)
        .map { sample_name, fastq_file, _transgene_name, size_range ->
            // Prefix with _ to suppress warning since we're not using transgene_name here
            [sample_name, fastq_file, size_range]
        }

    // Run CHOPPER for each size range
    CHOPPER(input_combinations)

    // Run NANOPLOT on the original input files
    original_input_for_nanoplot = input_ch.map { sample_name, file, _transgene_name ->
        // Prefix with _ to suppress warning since we're not using transgene_name here
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

    // Run FLYE assembly on all filtered reads
    FLYE(CHOPPER.out.filtered_reads)

    // SIMPLE BLAST APPROACH - Just use the default transgene for all samples for now
    transgene_fasta = file("${params.transgene_dir}/A-vector_herceptin_pEY345.fa", checkIfExists: true)
    
    // Prepare BLAST input
    blast_input_ch = FLYE.out.assembly_fasta
        .map { sample_name_size, assembly_fasta ->
            [sample_name_size, assembly_fasta, "A-vector_herceptin_pEY345", transgene_fasta]
        }

    // Run TRANSGENE_BLAST
    TRANSGENE_BLAST(blast_input_ch)

    // Gather assembly statistics from FLYE outputs
    assembly_info_with_names = FLYE.out.assembly_info
        .map { sample_name, file ->
            [sample_name, file]
        }
    flye_log_with_names = FLYE.out.flye_log
        .map { sample_name, file ->
            [sample_name, file]
        }

    // Collect all files but preserve sample names in filenames
    all_assembly_info = assembly_info_with_names
        .map { sample_name, file ->
            file.copyTo("${sample_name}.assembly_info.txt")
        }
        .collect()
    all_flye_logs = flye_log_with_names
        .map { sample_name, file ->
            file.copyTo("${sample_name}.flye.log")
        }
        .collect()

    GATHER_ASSEMBLY_STATS(all_assembly_info, all_flye_logs)
}