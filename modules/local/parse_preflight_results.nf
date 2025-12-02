process PARSE_PREFLIGHT_RESULTS {
    publishDir "${params.outdir}/summary", mode: 'copy'
    
    input:
    path preflight_logs
    path parse_script
    
    output:
    path "flye_preflight_summary.csv", emit: preflight_csv
    path "flye_preflight_summary.json", emit: preflight_json
    path "versions.yml", emit: versions
    
    script:
    """
    python3 ${parse_script}
    """
    
    stub:
    """
    echo "BaseSample,FullSample,Category,TotalReadLength,EstimatedCoverage,ReadsN50,ReadsN90" > flye_preflight_summary.csv
    echo "sample1,sample1_original,original,1000000,10.5,15000,8000" >> flye_preflight_summary.csv
    
    echo '[]' > flye_preflight_summary.json
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: "3.9.0"
    END_VERSIONS
    """
}