process DOWNSAMPLE_FASTQ {
    publishDir "${params.outdir}/selected_fastq", mode: 'copy'
    
    input:
    tuple val(sample_name), path(fastq_file), val(fraction)
    path downsample_script
    
    output:
    tuple val(sample_name), path("${sample_name}.fastq"), emit: downsampled_reads
    
    script:
    """
    python ${downsample_script} ${fastq_file} ${fraction} ${sample_name}.fastq
    """
    
    stub:
    """
    touch ${sample_name}.fastq
    """
}