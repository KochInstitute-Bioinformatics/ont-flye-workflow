process CHOPPER {
    tag "${sample_name}_${size_range.name}"
    publishDir "${params.outdir}/filtered_reads", mode: 'copy'
    
    input:
    tuple val(sample_name), path(input_fastq), val(size_range)
    
    output:
    tuple val("${sample_name}_${size_range.name}"), path("${sample_name}_${size_range.name}.fastq"), emit: filtered_reads
    
    script:
    def max_length_arg = size_range.max ? "--maxlength ${size_range.max}" : ""
    """
    chopper -q ${params.min_quality} \\
        -l ${size_range.min} \\
        ${max_length_arg} \\
        -i ${input_fastq} > ${sample_name}_${size_range.name}.fastq
    """
}