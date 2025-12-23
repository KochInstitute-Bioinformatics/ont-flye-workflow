#!/usr/bin/env python3

import os
import json
import re
import subprocess
import yaml
import glob
import argparse
import sys

def extract_assembly_stats(assembly_info_file, flye_log_file, sample_name):
    """Extract assembly statistics from Flye output files."""
    stats = {
        "sample_name": sample_name,
        "assembly_info": {},
        "flye_log": {}
    }

    # Extract from assembly_info.txt
    if os.path.exists(assembly_info_file):
        with open(assembly_info_file, 'r') as f:
            lines = f.readlines()
            if len(lines) > 1:  # Skip header line
                # Parse the assembly info table
                for line in lines[1:]:  # Skip header
                    if line.strip():
                        parts = line.strip().split('\t')
                        if len(parts) >= 4:
                            contig_name = parts[0]
                            length = int(parts[1])
                            coverage = float(parts[2])
                            circular = parts[3] == 'Y'

                            if contig_name not in stats["assembly_info"]:
                                stats["assembly_info"][contig_name] = {
                                    "length": length,
                                    "coverage": coverage,
                                    "circular": circular
                                }

    # Extract from flye.log
    if os.path.exists(flye_log_file):
        with open(flye_log_file, 'r') as f:
            content = f.read()

            # Extract key statistics from flye.log
            patterns = {
                'Total length': r'Total length:\s+([0-9,]+)',
                'Fragments': r'Fragments:\s+(\d+)',
                'Mean coverage': r'Mean coverage:\s+([0-9.]+)',
                'N50': r'N50:\s+([0-9,]+)',
                'N90': r'N90:\s+([0-9,]+)',
                'Largest contig': r'Largest contig:\s+([0-9,]+)'
            }

            for key, pattern in patterns.items():
                match = re.search(pattern, content)
                if match:
                    value = match.group(1).replace(',', '')
                    try:
                        # Try to convert to int first, then float
                        if '.' in value:
                            stats["flye_log"][key] = float(value)
                        else:
                            stats["flye_log"][key] = int(value)
                    except ValueError:
                        stats["flye_log"][key] = value

    return stats

def parse_sample_names(sample_names_str):
    """Parse sample names from Nextflow with multiple fallback methods."""
    print(f"Raw sample names string: {sample_names_str}")
    
    sample_names_list = []

    # Method 1: Try ast.literal_eval
    try:
        import ast
        sample_names_list = ast.literal_eval(sample_names_str)
        print(f"Method 1 - ast.literal_eval successful: {sample_names_list}")
        return sample_names_list
    except Exception as e:
        print(f"Method 1 - ast.literal_eval failed: {e}")
    
    # Method 2: Try JSON parsing
    try:
        sample_names_list = json.loads(sample_names_str.replace("'", '"'))
        print(f"Method 2 - JSON parsing successful: {sample_names_list}")
        return sample_names_list
    except Exception as e:
        print(f"Method 2 - JSON parsing failed: {e}")
    
    # Method 3: Manual string parsing
    try:
        # Remove brackets and split by comma
        clean_str = sample_names_str.strip('[]')
        if clean_str:
            sample_names_list = [name.strip().strip("'").strip('"') for name in clean_str.split(',')]
            print(f"Method 3 - Manual parsing successful: {sample_names_list}")
            return sample_names_list
        else:
            return []
    except Exception as e:
        print(f"Method 3 - Manual parsing failed: {e}")
        return []

def main():
    parser = argparse.ArgumentParser(description='Gather assembly statistics from Flye outputs')
    parser.add_argument('--sample-names', required=True, help='Sample names as string')
    parser.add_argument('--output-json', default='assembly_summary.json', help='Output JSON file')
    parser.add_argument('--output-versions', default='versions.yml', help='Output versions file')
    
    args = parser.parse_args()
    
    # Parse sample names
    sample_names_list = parse_sample_names(args.sample_names)
    
    # Find all staged files
    assembly_info_files = sorted(glob.glob("assembly_info_*.txt"))
    flye_log_files = sorted(glob.glob("flye_log_*.log"))

    print(f"Found {len(assembly_info_files)} assembly_info files: {assembly_info_files}")
    print(f"Found {len(flye_log_files)} flye_log files: {flye_log_files}")
    print(f"Sample names count: {len(sample_names_list)}")
    print(f"Final sample names list: {sample_names_list}")

    all_assembly_stats = {}

    # Process files in order with corresponding sample names
    for i in range(min(len(assembly_info_files), len(flye_log_files))):
        assembly_file = assembly_info_files[i]
        log_file = flye_log_files[i]
        
        # Use provided sample name or generate one
        if i < len(sample_names_list) and sample_names_list[i]:
            sample_name = sample_names_list[i]
        else:
            sample_name = f"sample_{i+1}"
            print(f"Warning: Using fallback name {sample_name} for index {i}")
        
        print(f"Processing {sample_name} (files: {assembly_file}, {log_file})...")
        
        if os.path.exists(assembly_file) and os.path.exists(log_file):
            stats = extract_assembly_stats(assembly_file, log_file, sample_name)
            all_assembly_stats[sample_name] = stats
        else:
            print(f"Warning: Files not found for {sample_name}")

    # Write consolidated JSON file
    with open(args.output_json, 'w') as json_file:
        json.dump(all_assembly_stats, json_file, indent=4)

    print(f"Assembly statistics extracted for {len(all_assembly_stats)} samples")
    print(f"Output written to: {args.output_json}")

    # Create versions file
    try:
        python_version = subprocess.check_output(['python', '--version'],
                                               stderr=subprocess.STDOUT,
                                               text=True).strip().replace('Python ', '')
    except Exception:
        python_version = "unknown"

    versions_data = {
        "ONT_FLYE:GATHER_ASSEMBLY_STATS": {
            "python": python_version
        }
    }

    with open(args.output_versions, "w") as f:
        yaml.dump(versions_data, f, default_flow_style=False)

if __name__ == "__main__":
    main()