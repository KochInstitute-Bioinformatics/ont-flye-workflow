include { DOWNSAMPLE_FASTQ } from '../modules/downsample_fastq'
include { NANOPLOT } from '../modules/nanoplot'
include { PARSE_NANOSTATS } from '../modules/parse_nanostats'
include { FLYE } from '../modules/flye'
include { GATHER_ASSEMBLY_STATS } from '../modules/gather_assembly_stats'
include { TRANSGENE_BLAST } from '../modules/transgene_blast'

workflow DOWNSAMPLE_MODE {
    take:
    input_ch  // [sample_name, fastq_file, transgene_name]
    
    main:
    // Define downsampling rates
    downsample_rates = Channel.fromList(params.downsample_rates)
    
    // Get the downsample script
    downsample_script = file("${projectDir}/bin/downsample_fastq.py", checkIfExists: true)
    
    // Create combinations of samples and downsample rates
    downsample_input = input_ch
        .map { sample_name, fastq_file, _transgene_name ->
            [sample_name, fastq_file]
        }
        .combine(downsample_rates)
    
    // Downsample the FASTQ files
    DOWNSAMPLE_FASTQ(downsample_input, downsample_script)
    
    // Run NANOPLOT on each downsampled file
    NANOPLOT(DOWNSAMPLE_FASTQ.out.downsampled_reads)
    
    // Collect all NanoPlot results
    all_nanoplot_results = NANOPLOT.out.nanoplot_results.collect()
    
    // Parse all NanoStats files and create summary table
    PARSE_NANOSTATS(all_nanoplot_results)
    
    // Run FLYE assembly on downsampled reads
    FLYE(DOWNSAMPLE_FASTQ.out.downsampled_reads)
    
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
        .map { sample_name_ds, assembly_fasta ->
            [sample_name_ds, assembly_fasta, "A-vector_herceptin_pEY345", transgene_fasta]
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
    downsampled_reads = DOWNSAMPLE_FASTQ.out.downsampled_reads
    nanoplot_results = all_nanoplot_results
    assembly_stats = GATHER_ASSEMBLY_STATS.out.assembly_stats
    blast_results = TRANSGENE_BLAST.out.blast_results
}