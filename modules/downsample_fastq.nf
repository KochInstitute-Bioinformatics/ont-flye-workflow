process DOWNSAMPLE_FASTQ {
    publishDir "${params.outdir}/selected_fastq", mode: 'symlink'
    
    input:
    tuple val(sample_name), path(fastq_file), val(fraction)
    path downsample_script
    
    output:
    tuple val(sample_name), path("${sample_name}.fastq"), emit: downsampled_reads
    
    script:
    """
    python ${downsample_script} ${fastq_file} ${sample_name}.fastq ${fraction}
    """
    
    stub:
    """
    touch ${sample_name}.fastq
    """
}