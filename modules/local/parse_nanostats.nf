process PARSE_NANOSTATS {
    publishDir "${params.outdir}/summary", mode: 'copy'
    
    input:
    path nanoplot_results
    path parse_script
    
    output:
    path "nanostats_summary.csv", emit: summary_csv
    path "nanostats_summary.json", emit: summary_json
    path "versions.yml", emit: versions
    
    script:
    """
    python ${parse_script}
    """
    
    stub:
    """
    touch nanostats_summary.csv
    touch nanostats_summary.json
    touch versions.yml
    """
}