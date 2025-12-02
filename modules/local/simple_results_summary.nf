process SIMPLE_RESULTS_SUMMARY {
    publishDir "${params.outdir}/summary", mode: 'copy'
    
    input:
    path nanostats_summary
    path preflight_summary
    path assembly_summary
    path transgene_count
    path simple_summary_script  // <-- Add this input
    
    output:
    path "simple_results_summary.csv", emit: summary_csv
    path "versions.yml", emit: versions
    
    script:
    """
    python ${simple_summary_script} \\
        --nanostats ${nanostats_summary} \\
        --preflight-summary ${preflight_summary} \\
        --assembly-summary ${assembly_summary} \\
        --transgene-count ${transgene_count} \\
        --output simple_results_summary.csv
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
        pandas: \$(python -c "import pandas; print(pandas.__version__)")
    END_VERSIONS
    """
}