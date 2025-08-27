process FILTER_ASSEMBLY_CANDIDATES {
    publishDir "${params.outdir}/summary", mode: 'copy'
    
    input:
    path preflight_summary_csv
    
    output:
    path "assembly_candidates.csv", emit: candidates_csv
    path "assembly_filtered.csv", emit: filtered_csv
    
    script:
    """
    #!/usr/bin/env python3
    import pandas as pd
    import csv
    
    # Read preflight results
    preflight_df = pd.read_csv('${preflight_summary_csv}')
    
    # Filter based on coverage criteria (10 < coverage < 300)
    candidates = []
    filtered_out = []
    
    for _, row in preflight_df.iterrows():
        sample_name = row['FullSample']
        coverage = float(row['EstimatedCoverage']) if row['EstimatedCoverage'] else 0
        
        if 10 < coverage < 300:
            candidates.append({
                'sample_name': sample_name,
                'estimated_coverage': coverage,
                'status': 'selected_for_assembly'
            })
        else:
            reason = 'coverage_too_low' if coverage <= 10 else 'coverage_too_high'
            filtered_out.append({
                'sample_name': sample_name,
                'estimated_coverage': coverage,
                'status': f'filtered_out_{reason}'
            })
    
    # Write candidates file
    with open('assembly_candidates.csv', 'w', newline='') as f:
        if candidates:
            writer = csv.DictWriter(f, fieldnames=['sample_name', 'estimated_coverage', 'status'])
            writer.writeheader()
            writer.writerows(candidates)
        else:
            # Create empty file with headers
            writer = csv.DictWriter(f, fieldnames=['sample_name', 'estimated_coverage', 'status'])
            writer.writeheader()
    
    # Write filtered out file
    with open('assembly_filtered.csv', 'w', newline='') as f:
        if filtered_out:
            writer = csv.DictWriter(f, fieldnames=['sample_name', 'estimated_coverage', 'status'])
            writer.writeheader()
            writer.writerows(filtered_out)
        else:
            # Create empty file with headers
            writer = csv.DictWriter(f, fieldnames=['sample_name', 'estimated_coverage', 'status'])
            writer.writeheader()
    
    print(f"Selected {len(candidates)} samples for assembly")
    print(f"Filtered out {len(filtered_out)} samples")
    """
    
    stub:
    """
    echo "sample_name,estimated_coverage,status" > assembly_candidates.csv
    echo "test_sample,50.5,selected_for_assembly" >> assembly_candidates.csv
    echo "sample_name,estimated_coverage,status" > assembly_filtered.csv
    echo "bad_sample,5.2,filtered_out_coverage_too_low" >> assembly_filtered.csv
    """
}