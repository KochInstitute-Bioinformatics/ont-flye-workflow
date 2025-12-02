process PARSE_TRANSGENE_BLAST {
    
    publishDir "${params.outdir}/summary", mode: 'copy'
    
    input:
    path blast_files
    path transgene_library
    path parse_script
    
    output:
    path "transgene_count.json", emit: json_results
    path "transgene_count.csv", emit: csv_results
    path "versions.yml", emit: versions
    
    when:
    task.ext.when == null || task.ext.when
    
    script:
    def args = task.ext.args ?: ''
    """
    python ${parse_script} \\
        --blast-files ${blast_files} \\
        --transgene-library ${transgene_library} \\
        --output-dir . \\
        --tolerance 100 \\
        ${args}
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
    
    stub:
    """
    touch transgene_count.json
    touch transgene_count.csv
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}