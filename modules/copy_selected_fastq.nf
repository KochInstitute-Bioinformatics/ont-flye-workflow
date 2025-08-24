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