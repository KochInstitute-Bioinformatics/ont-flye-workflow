include { CHOPPER } from '../modules/chopper'
include { NANOPLOT } from '../modules/nanoplot'

workflow ONT_FLYE {  // ← Changed from NANOPORE_QC to ONT_FLYE
    
    // Create input channel
    input_ch = Channel.fromPath(params.input_fastq, checkIfExists: true)
    
    // Create size ranges channel
    size_ranges_ch = Channel.fromList(params.size_ranges)
    
    // Run CHOPPER for each size range
    CHOPPER(input_ch, size_ranges_ch)
    
    // Run NANOPLOT on each filtered file
    NANOPLOT(CHOPPER.out.filtered_reads)
    
    // Create combined output for 40k+ reads (combining 40k_50k and 50k_Plus)
    combined_40k_plus = CHOPPER.out.filtered_reads
        .filter { sample_name, fastq -> 
            sample_name == "40k_50k" || sample_name == "50k_Plus" 
        }
        .map { sample_name, fastq -> fastq }
        .collectFile(name: "sHF171_40k_Plus.fastq", newLine: false)
        .map { file -> ["40k_Plus", file] }
    
    // Run NANOPLOT on combined 40k+ file
    NANOPLOT(combined_40k_plus)
}