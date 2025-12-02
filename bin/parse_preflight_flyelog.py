#!/usr/bin/env python3
import os
import glob
import json
import csv
import re
import yaml
import sys

def extract_base_sample_name(sample_name):
    """Extract the base sample name without suffixes"""
    base_name = re.sub(r'_original.*', '', sample_name)
    base_name = re.sub(r'_\d+k_Plus.*', '', base_name)
    base_name = re.sub(r'_ds[\d.]+.*', '', base_name)
    base_name = re.sub(r'_bs[\d.]+.*', '', base_name)
    return base_name

def categorize_sample(sample_name):
    """Categorize sample based on naming convention"""
    if '_original' in sample_name:
        if '_ds' in sample_name:
            match = re.search(r'_ds([\d.]+)', sample_name)
            rate = match.group(1) if match else 'unknown'
            return f"downsample_{rate}"
        elif '_bs' in sample_name:
            match = re.search(r'_bs([\d.]+)_rep(\d+)', sample_name)
            if match:
                fraction, rep = match.groups()
                return f"bootstrap_{fraction}_rep{rep}"
            return "bootstrap"
        else:
            return "original"
    elif re.search(r'_\d+k_Plus', sample_name):
        match = re.search(r'_(\d+k)_Plus', sample_name)
        size = match.group(1) if match else 'unknown'
        if '_ds' in sample_name:
            ds_match = re.search(r'_ds([\d.]+)', sample_name)
            rate = ds_match.group(1) if ds_match else 'unknown'
            return f"filtered_{size}_downsample_{rate}"
        elif '_bs' in sample_name:
            bs_match = re.search(r'_bs([\d.]+)_rep(\d+)', sample_name)
            if bs_match:
                fraction, rep = bs_match.groups()
                return f"filtered_{size}_bootstrap_{fraction}_rep{rep}"
            return f"filtered_{size}_bootstrap"
        else:
            return f"filtered_{size}"
    else:
        return "unknown"

def parse_flye_log(file_path):
    """Parse a single flye.log file to extract preflight information"""
    sample_name = os.path.basename(file_path).replace('_flye.log', '')
    base_sample = extract_base_sample_name(sample_name)
    category = categorize_sample(sample_name)
    
    total_read_length = estimated_coverage = reads_n50 = reads_n90 = None
    
    with open(file_path, 'r') as f:
        content = f.read()
        
        # Extract Total read length
        total_length_match = re.search(r'Total read length:\s+([0-9,]+)', content)
        if total_length_match:
            total_read_length = total_length_match.group(1).replace(',', '')
        
        # Extract Estimated coverage
        coverage_match = re.search(r'Estimated coverage:\s+([0-9.]+)', content)
        if coverage_match:
            estimated_coverage = coverage_match.group(1)
        
        # Extract Reads N50/N90
        n50_n90_match = re.search(r'Reads N50/N90:\s+([0-9,]+)\s*/\s*([0-9,]+)', content)
        if n50_n90_match:
            reads_n50 = n50_n90_match.group(1).replace(',', '')
            reads_n90 = n50_n90_match.group(2).replace(',', '')
    
    return {
        "BaseSample": base_sample,
        "FullSample": sample_name,
        "Category": category,
        "TotalReadLength": total_read_length,
        "EstimatedCoverage": estimated_coverage,
        "ReadsN50": reads_n50,
        "ReadsN90": reads_n90
    }

def main():
    # Find all flye.log files
    flye_log_files = glob.glob("*_flye.log")
    print(f"Found {len(flye_log_files)} flye.log files")
    
    if not flye_log_files:
        print("No flye.log files found in current directory")
        sys.exit(1)
    
    # Parse all files and collect results
    results = [parse_flye_log(file_path) for file_path in flye_log_files]
    
    # Sort results by base sample name, then by category
    results.sort(key=lambda x: (x["BaseSample"], x["Category"]))
    
    # Write comprehensive summary to JSON file
    with open("flye_preflight_summary.json", "w") as f:
        json.dump(results, f, indent=4)
    
    # Write CSV summary for easier analysis
    if results:
        with open("flye_preflight_summary.csv", "w", newline='') as f:
            fieldnames = ["BaseSample", "FullSample", "Category", "TotalReadLength", 
                         "EstimatedCoverage", "ReadsN50", "ReadsN90"]
            writer = csv.DictWriter(f, fieldnames=fieldnames)
            writer.writeheader()
            writer.writerows(results)
    
    # Create versions file
    import subprocess
    try:
        python_version = subprocess.check_output(['python', '--version'], 
                                               stderr=subprocess.STDOUT, 
                                               text=True).strip().replace('Python ', '')
    except:
        python_version = "unknown"
    
    versions_data = {
        "ONT_FLYE:PARSE_PREFLIGHT_FLYELOG": {
            "python": python_version
        }
    }
    
    with open("versions.yml", "w") as f:
        yaml.dump(versions_data, f, default_flow_style=False)
    
    print(f"Processed {len(results)} flye preflight logs")
    print("Results saved to flye_preflight_summary.json and flye_preflight_summary.csv")

if __name__ == "__main__":
    main()