#!/usr/bin/env python3
import os
import glob
import json
import csv
import re
import yaml
import sys

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

def extract_base_sample_name(sample_name):
    """Extract the base sample name without suffixes"""
    base_name = re.sub(r'_original.*', '', sample_name)
    base_name = re.sub(r'_\d+k_Plus.*', '', base_name)
    return base_name

def parse_nanostats(file_path):
    """Parse a single NanoStats.txt file"""
    sample = os.path.basename(os.path.dirname(file_path))
    base_sample = extract_base_sample_name(sample)
    category = categorize_sample(sample)
    
    read_count = mean_length = n50 = total_bases = None
    
    with open(file_path, 'r') as f:
        for line in f:
            if line.startswith("Mean read length:"):
                mean_length = line.split(":")[1].strip().replace(",", "")
            elif line.startswith("Number of reads:"):
                read_count = line.split(":")[1].strip().replace(",", "")
            elif line.startswith("Read length N50:"):
                n50 = line.split(":")[1].strip().replace(",", "")
            elif line.startswith("Total bases:"):
                total_bases = line.split(":")[1].strip().replace(",", "")
    
    return {
        "BaseSample": base_sample,
        "FullSample": sample,
        "Category": category,
        "ReadCount": read_count,
        "MeanLength": mean_length,
        "N50": n50,
        "TotalBases": total_bases
    }

def main():
    # Find all NanoStats.txt files
    nanostats_files = glob.glob("*/NanoStats.txt")
    print(f"Found {len(nanostats_files)} NanoStats.txt files")
    
    # Parse all files and collect results
    results = [parse_nanostats(file_path) for file_path in nanostats_files]
    
    # Sort results by base sample name, then by category
    results.sort(key=lambda x: (x["BaseSample"], x["Category"]))
    
    # Write comprehensive summary to JSON file
    with open("nanostats_summary.json", "w") as f:
        json.dump(results, f, indent=4)
    
    # Write CSV summary for easier analysis
    if results:
        with open("nanostats_summary.csv", "w", newline='') as f:
            fieldnames = ["BaseSample", "FullSample", "Category", "ReadCount", "MeanLength", "N50", "TotalBases"]
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
        "ONT_FLYE:PARSE_NANOSTATS": {
            "python": python_version
        }
    }
    
    with open("versions.yml", "w") as f:
        yaml.dump(versions_data, f, default_flow_style=False)
    
    print(f"Processed {len(results)} samples across all workflow modes")

if __name__ == "__main__":
    main()