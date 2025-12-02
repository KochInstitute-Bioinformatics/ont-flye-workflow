process CHOPPER {
    publishDir "${params.outdir}/selected_fastq", mode: 'copy'
    
    input:
    tuple val(sample_name), path(fastq_file), val(size_range)
    
    output:
    tuple val("${sample_name}_${size_range.name}"), path("${sample_name}_${size_range.name}.fastq"), emit: filtered_reads
    
    script:
    def min_length = size_range.min ?: 0
    def max_length = size_range.max ? "-u ${size_range.max}" : ""
    def quality = params.min_quality ?: 10
    
    """
    chopper -q ${quality} \\
        -l ${min_length} \\
        ${max_length} \\
        -i ${fastq_file} > ${sample_name}_${size_range.name}.fastq
    """
    
    stub:
    """
    touch ${sample_name}_${size_range.name}.fastq
    """
}