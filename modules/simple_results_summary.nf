process SIMPLE_RESULTS_SUMMARY {
    publishDir "${params.outdir}/summary", mode: 'copy'
    
    input:
    path nanostats_summary
    path assembly_filtered
    path assembly_candidates  
    path assembly_summary
    path transgene_count
    
    output:
    path "simple_results_summary.csv", emit: summary_csv
    path "versions.yml", emit: versions
    
    script:
    """
    simple_summary.py \\
        --nanostats ${nanostats_summary} \\
        --assembly-filtered ${assembly_filtered} \\
        --assembly-candidates ${assembly_candidates} \\
        --assembly-summary ${assembly_summary} \\
        --transgene-count ${transgene_count} \\
        --output simple_results_summary.csv
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}