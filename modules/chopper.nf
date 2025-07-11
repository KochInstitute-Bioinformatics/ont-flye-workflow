process CHOPPER {
    tag "${size_range.name}"
    publishDir "${params.outdir}/filtered_reads", mode: 'copy'
    
    input:
    path input_fastq
    val size_range
    
    output:
    tuple val(size_range.name), path("sHF171_${size_range.name}.fastq"), emit: filtered_reads
    
    script:
    def max_length_arg = size_range.max ? "--maxlength ${size_range.max}" : ""
    """
    chopper -q ${params.min_quality} \\
            -l ${size_range.min} \\
            ${max_length_arg} \\
            -i ${input_fastq} > sHF171_${size_range.name}.fastq
    """
}