#!/usr/bin/env python3

import json
import pandas as pd
import argparse
import sys
from pathlib import Path

def parse_nanostats_summary(file_path):
    """Parse nanostats_summary.json for BaseSample, FullSample, Category, ReadCount, MeanLength"""
    try:
        with open(file_path, 'r') as f:
            data = json.load(f)
        
        results = []
        for entry in data:
            results.append({
                'BaseSample': entry.get('BaseSample', ''),
                'FullSample': entry.get('FullSample', ''),
                'Category': entry.get('Category', ''),
                'ReadCount': entry.get('ReadCount', 0),
                'MeanLength': entry.get('MeanLength', 0.0)
            })
        return pd.DataFrame(results)
    except Exception as e:
        print(f"Error parsing nanostats_summary.json: {e}", file=sys.stderr)
        return pd.DataFrame()

def parse_preflight_summary(file_path):
    """Parse flye_preflight_summary.json for coverage data"""
    try:
        with open(file_path, 'r') as f:
            data = json.load(f)
        
        results = []
        for entry in data:
            results.append({
                'FullSample': entry.get('sample_name', ''),
                'estimated_coverage': entry.get('estimated_coverage', 0)
            })
        return pd.DataFrame(results)
    except Exception as e:
        print(f"Error parsing flye_preflight_summary.json: {e}", file=sys.stderr)
        return pd.DataFrame()

def parse_assembly_summary(file_path):
    """Parse assembly_summary.json for Fragments and Mean coverage from flye.log section"""
    try:
        with open(file_path, 'r') as f:
            data = json.load(f)
        
        results = []
        for entry in data:
            flye_log = entry.get('flye.log', {})
            results.append({
                'FullSample': entry.get('sample_name', ''),
                'Fragments': flye_log.get('Fragments', 0),
                'Mean_coverage': flye_log.get('Mean coverage', 0.0)
            })
        return pd.DataFrame(results)
    except Exception as e:
        print(f"Error parsing assembly_summary.json: {e}", file=sys.stderr)
        return pd.DataFrame()

def parse_transgene_count(file_path):
    """Parse transgene_count.json for full_length_count and contig_names"""
    try:
        with open(file_path, 'r') as f:
            data = json.load(f)
        
        results = []
        for entry in data:
            # Handle contig_names - join list into string if it's a list
            contig_names = entry.get('contig_names', [])
            if isinstance(contig_names, list):
                contig_names_str = ';'.join(contig_names)
            else:
                contig_names_str = str(contig_names)
            
            results.append({
                'FullSample': entry.get('sample_name', ''),
                'full_length_count': entry.get('full_length_count', 0),
                'contig_names': contig_names_str
            })
        return pd.DataFrame(results)
    except Exception as e:
        print(f"Error parsing transgene_count.json: {e}", file=sys.stderr)
        return pd.DataFrame()

def main():
    parser = argparse.ArgumentParser(description='Create simple results summary CSV')
    parser.add_argument('--nanostats', required=True, help='Path to nanostats_summary.json')
    parser.add_argument('--preflight-summary', required=True, help='Path to flye_preflight_summary.json')
    parser.add_argument('--assembly-summary', required=True, help='Path to assembly_summary.json')
    parser.add_argument('--transgene-count', required=True, help='Path to transgene_count.json')
    parser.add_argument('--output', required=True, help='Output CSV file path')
    
    args = parser.parse_args()
    
    # Parse all input files
    print("Parsing nanostats summary...")
    nanostats_df = parse_nanostats_summary(args.nanostats)
    
    print("Parsing preflight summary...")
    preflight_df = parse_preflight_summary(args.preflight_summary)
    
    print("Parsing assembly summary...")
    assembly_summary_df = parse_assembly_summary(args.assembly_summary)
    
    print("Parsing transgene count...")
    transgene_df = parse_transgene_count(args.transgene_count)
    
    # Start with nanostats data as the base
    if nanostats_df.empty:
        print("Warning: No nanostats data found", file=sys.stderr)
        final_df = pd.DataFrame()
    else:
        final_df = nanostats_df.copy()
    
    # Join with preflight data
    if not preflight_df.empty and not final_df.empty:
        final_df = final_df.merge(preflight_df, on='FullSample', how='left')
    
    # Join with assembly summary data
    if not assembly_summary_df.empty and not final_df.empty:
        final_df = final_df.merge(assembly_summary_df, on='FullSample', how='left')
    
    # Join with transgene data
    if not transgene_df.empty and not final_df.empty:
        final_df = final_df.merge(transgene_df, on='FullSample', how='left')
    
    # Fill NaN values with appropriate defaults
    numeric_columns = ['ReadCount', 'MeanLength', 'estimated_coverage', 'Fragments', 'Mean_coverage', 'full_length_count']
    for col in numeric_columns:
        if col in final_df.columns:
            final_df[col] = final_df[col].fillna(0)
    
    string_columns = ['BaseSample', 'FullSample', 'Category', 'contig_names']
    for col in string_columns:
        if col in final_df.columns:
            final_df[col] = final_df[col].fillna('')
    
    # Save to CSV
    if not final_df.empty:
        final_df.to_csv(args.output, index=False)
        print(f"Summary saved to {args.output}")
        print(f"Total rows: {len(final_df)}")
    else:
        print("Warning: No data to write to summary file", file=sys.stderr)
        # Create empty CSV with headers
        empty_df = pd.DataFrame(columns=[
            'BaseSample', 'FullSample', 'Category', 'ReadCount', 'MeanLength',
            'estimated_coverage', 'Fragments', 'Mean_coverage',
            'full_length_count', 'contig_names'
        ])
        empty_df.to_csv(args.output, index=False)

if __name__ == '__main__':
    main()