#!/usr/bin/env python3
"""
Filter assembly candidates based on coverage criteria
Usage: filter_assembly_candidates.py <preflight_summary.csv> <min_depth> <max_depth>
"""

import sys
import pandas as pd
import csv

def filter_candidates(preflight_csv, min_depth, max_depth):
    """Filter samples based on coverage thresholds"""
    
    # Read preflight results
    preflight_df = pd.read_csv(preflight_csv)
    
    print(f"Filtering samples with coverage thresholds: {min_depth} < coverage < {max_depth}")
    
    candidates = []
    filtered_out = []
    
    for _, row in preflight_df.iterrows():
        sample_name = row['FullSample']
        coverage = float(row['EstimatedCoverage']) if row['EstimatedCoverage'] else 0
        
        if min_depth < coverage < max_depth:
            candidates.append({
                'sample_name': sample_name,
                'estimated_coverage': coverage,
                'status': 'selected_for_assembly'
            })
        else:
            reason = 'coverage_too_low' if coverage <= min_depth else 'coverage_too_high'
            filtered_out.append({
                'sample_name': sample_name,
                'estimated_coverage': coverage,
                'status': f'filtered_out_{reason}'
            })
    
    # Write candidates file
    with open('assembly_candidates.csv', 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=['sample_name', 'estimated_coverage', 'status'])
        writer.writeheader()
        if candidates:
            writer.writerows(candidates)
    
    # Write filtered out file
    with open('assembly_filtered.csv', 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=['sample_name', 'estimated_coverage', 'status'])
        writer.writeheader()
        if filtered_out:
            writer.writerows(filtered_out)
    
    print(f"Selected {len(candidates)} samples for assembly")
    print(f"Filtered out {len(filtered_out)} samples")

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: filter_assembly_candidates.py <preflight_summary.csv> <min_depth> <max_depth>", file=sys.stderr)
        sys.exit(1)
    
    preflight_csv = sys.argv[1]
    min_depth = float(sys.argv[2])
    max_depth = float(sys.argv[3])
    
    filter_candidates(preflight_csv, min_depth, max_depth)