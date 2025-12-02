process GATHER_ASSEMBLY_STATS {
    publishDir "${params.outdir}/assembly_summary", mode: 'copy'
    
    input:
    path assembly_info_files, stageAs: "assembly_info_*.txt"
    path flye_log_files, stageAs: "flye_log_*.log"
    val sample_names
    path gather_script
    
    output:
    path "assembly_summary.json", emit: assembly_stats
    path "versions.yml", emit: versions
    
    script:
    """
    python ${gather_script} \\
        --sample-names '${sample_names}' \\
        --output-json assembly_summary.json \\
        --output-versions versions.yml
    """
    
    stub:
    """
    echo '{}' > assembly_summary.json
    echo 'ONT_FLYE:GATHER_ASSEMBLY_STATS:' > versions.yml
    """
}