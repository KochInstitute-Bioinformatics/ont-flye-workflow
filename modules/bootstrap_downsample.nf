process BOOTSTRAP_DOWNSAMPLE {
    conda '/net/ostrom/data/bcc/charliew/.conda/envs/nf-core_v24'
    
    publishDir "${params.outdir}/bootstrap_reads", mode: 'copy'
    
    input:
    tuple val(sample_name), path(fastq_file), val(fraction), val(replicate)
    path downsample_script
    
    output:
    tuple val("${sample_name}_bs${fraction}_rep${replicate}"), path("${sample_name}_bs${fraction}_rep${replicate}.fastq"), emit: bootstrap_reads
    path "versions.yml", emit: versions
    
    script:
    """
    python3 ${downsample_script} \\
        ${fastq_file} \\
        ${fraction} \\
        ${sample_name}_bs${fraction}_rep${replicate}.fastq
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //')
        biopython: \$(python -c "import Bio; print(Bio.__version__)")
    END_VERSIONS
    """
    
    stub:
    """
    touch ${sample_name}_bs${fraction}_rep${replicate}.fastq
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: "3.9.0"
        biopython: "1.81"
    END_VERSIONS
    """
}