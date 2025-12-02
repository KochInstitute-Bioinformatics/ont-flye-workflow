process FLYE_PREFLIGHT {
    tag "${sample_name}"
    publishDir "${params.outdir}/preflight_assemblies", mode: 'copy'
    
    input:
    tuple val(sample_name), path(fastq_file)
    
    output:
    tuple val(sample_name), path("${sample_name}_flye.log"), emit: preflight_logs
    path "versions.yml", emit: versions
    
    script:
    """
    flye --nano-raw ${fastq_file} \\
         --genome-size ${params.genome_size} \\
         --threads ${task.cpus} \\
         --stop-after configure \\
         --out-dir flye_temp_${sample_name} > ${sample_name}_flye.log 2>&1
    
    # Clean up temporary directory
    rm -rf flye_temp_${sample_name}
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        flye: \$(flye --version | head -n1 | sed 's/^.*flye //')
    END_VERSIONS
    """
    
    stub:
    """
    echo "Mock flye preflight log for ${sample_name}" > ${sample_name}_flye.log
    echo "Total read length: 1,000,000" >> ${sample_name}_flye.log
    echo "Estimated coverage: 10.5" >> ${sample_name}_flye.log
    echo "Reads N50/N90: 15,000 / 8,000" >> ${sample_name}_flye.log
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        flye: "2.9.1"
    END_VERSIONS
    """
}