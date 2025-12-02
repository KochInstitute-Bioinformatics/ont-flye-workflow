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
                'ReadCount': float(entry.get('ReadCount', 0)),  # Convert to float
                'MeanLength': float(entry.get('MeanLength', 0.0))  # Convert to float
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
                'FullSample': entry.get('FullSample', ''),  # Note: using FullSample, not sample_name
                'estimated_coverage': float(entry.get('EstimatedCoverage', 0))  # Note: EstimatedCoverage not estimated_coverage
            })
        
        return pd.DataFrame(results)
    except Exception as e:
        print(f"Error parsing flye_preflight_summary.json: {e}", file=sys.stderr)
        return pd.DataFrame()

def parse_assembly_summary(file_path):
    """Parse assembly_summary.json for Fragments, Mean coverage, and N50 from flye_log section"""
    try:
        with open(file_path, 'r') as f:
            data = json.load(f)
        
        results = []
        # Handle the nested structure - data is a dict with sample names as keys
        for sample_key, sample_data in data.items():
            # Get the sample_name from within the sample_data
            sample_name = sample_data.get('sample_name', sample_key)
            
            # Look for flye_log data (note: underscore, not dot!)
            flye_log = sample_data.get('flye_log', {})
            
            # Extract Fragments, Mean coverage, and N50
            fragments = 0
            mean_coverage = 0.0
            n50 = 0
            
            if isinstance(flye_log, dict):
                fragments = int(flye_log.get('Fragments', 0)) if flye_log.get('Fragments') else 0
                mean_coverage = float(flye_log.get('Mean coverage', 0.0)) if flye_log.get('Mean coverage') else 0.0
                n50 = int(flye_log.get('N50', 0)) if flye_log.get('N50') else 0
            
            results.append({
                'FullSample': sample_name,
                'Fragments': fragments,
                'Mean_coverage': mean_coverage,
                'N50': n50
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
        # Group by assembly_name since there can be multiple transgenes per assembly
        assembly_data = {}
        
        for key, entry in data.items():
            assembly_name = entry.get('assembly_name', '')
            if assembly_name not in assembly_data:
                assembly_data[assembly_name] = {
                    'FullSample': assembly_name,
                    'full_length_count': 0,
                    'contig_names': set()
                }
            
            # Sum up full_length_count across all transgenes for this assembly
            assembly_data[assembly_name]['full_length_count'] += int(entry.get('full_length_count', 0))
            
            # Collect all contig names
            contig_names = entry.get('contig_names', '')
            if contig_names:
                if isinstance(contig_names, list):
                    assembly_data[assembly_name]['contig_names'].update(contig_names)
                else:
                    assembly_data[assembly_name]['contig_names'].add(str(contig_names))
        
        # Convert to list format
        for assembly_name, data_dict in assembly_data.items():
            results.append({
                'FullSample': data_dict['FullSample'],
                'full_length_count': data_dict['full_length_count'],
                'contig_names': ';'.join(sorted(data_dict['contig_names']))
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
    numeric_columns = ['ReadCount', 'MeanLength', 'estimated_coverage', 'Fragments', 'Mean_coverage', 'N50', 'full_length_count']
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
            'estimated_coverage', 'Fragments', 'Mean_coverage', 'N50',
            'full_length_count', 'contig_names'
        ])
        empty_df.to_csv(args.output, index=False)

if __name__ == '__main__':
    main()