include { DOWNSAMPLE_FASTQ } from '../modules/downsample_fastq'
include { NANOPLOT } from '../modules/nanoplot'
include { NANOPLOT as NANOPLOT_ORIGINAL } from '../modules/nanoplot'
include { PARSE_NANOSTATS } from '../modules/parse_nanostats'
include { FLYE } from '../modules/flye'
include { GATHER_ASSEMBLY_STATS } from '../modules/gather_assembly_stats'
include { TRANSGENE_BLAST } from '../modules/transgene_blast'

workflow DOWNSAMPLE_MODE {
    take:
    input_ch  // [sample_name, fastq_file, transgene_name]
    
    main:
    // Get the downsample script
    downsample_script = file("${projectDir}/bin/downsample_fastq.py", checkIfExists: true)
    
    // Create fractions channel
    fractions_ch = Channel.fromList(params.fractions)
    
    // Combine input with each fraction to create all combinations
    input_combinations = input_ch.combine(fractions_ch)
        .map { sample_name, fastq_file, _transgene_name, fraction ->
            [sample_name, fastq_file, fraction]
        }
    
    // Run DOWNSAMPLE_FASTQ for each fraction
    DOWNSAMPLE_FASTQ(input_combinations, downsample_script)
    
    // Run NANOPLOT on the original input files
    original_input_for_nanoplot = input_ch.map { sample_name, fastq_file, _transgene_name ->
        tuple("${sample_name}_original", fastq_file)
    }
    NANOPLOT_ORIGINAL(original_input_for_nanoplot)
    
    // Run NANOPLOT on each downsampled file
    NANOPLOT(DOWNSAMPLE_FASTQ.out.downsampled_reads)
    
    // Collect all NanoPlot results and extract NanoStats.txt files
    all_nanoplot_results = NANOPLOT_ORIGINAL.out.nanoplot_results
        .mix(NANOPLOT.out.nanoplot_results)
        .collect()
    
    // Parse all NanoStats files and create summary table
    PARSE_NANOSTATS(all_nanoplot_results)
    
    // Run FLYE assembly on all downsampled reads
    FLYE(DOWNSAMPLE_FASTQ.out.downsampled_reads)
    
    // BLAST analysis - only on successful assemblies
    transgene_fasta = file("${params.transgene_dir}/A-vector_herceptin_pEY345.fa", checkIfExists: true)
    blast_input_ch = FLYE.out.assembly_fasta
        .map { sample_name_fraction, assembly_fasta ->
            [sample_name_fraction, assembly_fasta, "A-vector_herceptin_pEY345", transgene_fasta]
        }
    
    TRANSGENE_BLAST(blast_input_ch)
    
    // Gather assembly statistics - separate files and sample names
    assembly_info_collected = FLYE.out.assembly_info
        .map { _sample_name, assembly_file -> assembly_file }
        .collect()
    
    flye_log_collected = FLYE.out.flye_log
        .map { _sample_name, log_file -> log_file }
        .collect()
    
    sample_names_collected = FLYE.out.assembly_info
        .map { sample_name, _assembly_file -> sample_name }
        .collect()
    
    GATHER_ASSEMBLY_STATS(assembly_info_collected, flye_log_collected, sample_names_collected)
    
    emit:
    nanoplot_results = all_nanoplot_results
    assembly_stats = GATHER_ASSEMBLY_STATS.out.assembly_stats
    blast_results = TRANSGENE_BLAST.out.blast_results
}