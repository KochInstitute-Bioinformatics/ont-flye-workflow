process FILTER_ASSEMBLY_CANDIDATES {
    publishDir "${params.outdir}/summary", mode: 'copy'
    
    input:
    path preflight_summary_csv
    val min_depth
    val max_depth
    
    output:
    path "assembly_candidates.csv", emit: candidates_csv
    path "assembly_filtered.csv", emit: filtered_csv
    
    script:
    """
    filter_assembly_candidates.py ${preflight_summary_csv} ${min_depth} ${max_depth}
    """
    
    stub:
    """
    echo "sample_name,estimated_coverage,status" > assembly_candidates.csv
    echo "test_sample,50.5,selected_for_assembly" >> assembly_candidates.csv
    echo "sample_name,estimated_coverage,status" > assembly_filtered.csv
    echo "bad_sample,5.2,filtered_out_coverage_too_low" >> assembly_filtered.csv
    """
}