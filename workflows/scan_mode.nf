include { CHOPPER } from '../modules/chopper'
include { NANOPLOT } from '../modules/nanoplot'
include { NANOPLOT as NANOPLOT_ORIGINAL } from '../modules/nanoplot'
include { PARSE_NANOSTATS } from '../modules/parse_nanostats'
include { FLYE } from '../modules/flye'
include { GATHER_ASSEMBLY_STATS } from '../modules/gather_assembly_stats'
include { TRANSGENE_BLAST } from '../modules/transgene_blast'

workflow SCAN_MODE {
    take:
    input_ch  // [sample_name, fastq_file, transgene_name]
    
    main:
    // Create size ranges channel
    size_ranges_ch = Channel.fromList(params.size_ranges)

    // Combine input with each size range to create all combinations
    input_combinations = input_ch.combine(size_ranges_ch)
        .map { sample_name, fastq_file, _transgene_name, size_range ->
            [sample_name, fastq_file, size_range]
        }

    // Run CHOPPER for each size range
    CHOPPER(input_combinations)

    // Run NANOPLOT on the original input files
    original_input_for_nanoplot = input_ch.map { sample_name, fastq_file, _transgene_name ->
        tuple("${sample_name}_original", fastq_file)
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

    // Filter out failed assemblies before BLAST
    successful_assemblies = FLYE.out.assembly_fasta
        .filter { sample_name, assembly_fasta ->
            // Check if assembly file exists and is not a failure marker
            assembly_fasta.exists() && 
            assembly_fasta.size() > 100 &&  // Minimum size check
            !assembly_fasta.text.contains("ASSEMBLY_FAILED")
    }


    // BLAST analysis
    transgene_fasta = file("${params.transgene_dir}/A-vector_herceptin_pEY345.fa", checkIfExists: true)
    
    blast_input_ch = FLYE.out.assembly_fasta
        .map { sample_name_size, assembly_fasta ->
            [sample_name_size, assembly_fasta, "A-vector_herceptin_pEY345", transgene_fasta]
        }

    TRANSGENE_BLAST(blast_input_ch)

    // Gather assembly statistics
    assembly_info_with_names = FLYE.out.assembly_info
        .map { sample_name, assembly_file ->
            [sample_name, assembly_file]
        }
    flye_log_with_names = FLYE.out.flye_log
        .map { sample_name, log_file ->
            [sample_name, log_file]
        }

    all_assembly_info = assembly_info_with_names
        .map { sample_name, assembly_file ->
            assembly_file.copyTo("${sample_name}.assembly_info.txt")
        }
        .collect()
    all_flye_logs = flye_log_with_names
        .map { sample_name, log_file ->
            log_file.copyTo("${sample_name}.flye.log")
        }
        .collect()

    GATHER_ASSEMBLY_STATS(all_assembly_info, all_flye_logs)
    
    emit:
    nanoplot_results = all_nanoplot_results
    assembly_stats = GATHER_ASSEMBLY_STATS.out.assembly_stats
    blast_results = TRANSGENE_BLAST.out.blast_results
}