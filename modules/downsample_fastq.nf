process DOWNSAMPLE_FASTQ {
    conda "bioconda::biopython=1.81"
    
    publishDir "${params.outdir}/downsampled_reads", mode: 'copy'
    
    input:
    tuple val(sample_name), path(fastq_file), val(downsample_rate)
    
    output:
    tuple val("${sample_name}_ds${downsample_rate}"), path("${sample_name}_ds${downsample_rate}.fastq"), emit: downsampled_reads
    path "versions.yml", emit: versions
    
    script:
    """
    python3 ${projectDir}/bin/downsample_fastq.py \\
        ${fastq_file} \\
        ${downsample_rate} \\
        ${sample_name}_ds${downsample_rate}.fastq
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //')
        biopython: \$(python -c "import Bio; print(Bio.__version__)")
    END_VERSIONS
    """
    
    stub:
    """
    touch ${sample_name}_ds${downsample_rate}.fastq
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: "3.9.0"
        biopython: "1.81"
    END_VERSIONS
    """
}