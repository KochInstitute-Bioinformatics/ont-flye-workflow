process GATHER_ASSEMBLY_STATS {
    publishDir "${params.outdir}/assembly_summary", mode: 'copy'
    
    input:
    path assembly_info_files, stageAs: "assembly_info_*.txt"
    path flye_log_files, stageAs: "flye_log_*.log"
    val sample_names
    
    output:
    path "assembly_summary.json", emit: assembly_stats
    path "versions.yml", emit: versions
    
    script:
    """
    # Create and run Python script
    cat > gather_stats.py << 'EOF'
import os
import json
import re
import subprocess
import yaml
import glob

def extract_assembly_stats(assembly_info_file, flye_log_file, sample_name):
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
                        parts = line.strip().split('\\t')
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
                'Total length': r'Total length:\\s+([\\d,]+)',
                'Fragments': r'Fragments:\\s+(\\d+)',
                'Mean coverage': r'Mean coverage:\\s+([\\d.]+)',
                'N50': r'N50:\\s+([\\d,]+)',
                'N90': r'N90:\\s+([\\d,]+)',
                'Largest contig': r'Largest contig:\\s+([\\d,]+)'
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

# Get sample names from Nextflow - improved parsing
sample_names_str = '''${sample_names}'''
print(f"Raw sample names string: {sample_names_str}")

# Parse sample names with multiple fallback methods
sample_names_list = []

# Method 1: Try ast.literal_eval
try:
    import ast
    sample_names_list = ast.literal_eval(sample_names_str)
    print(f"Method 1 - ast.literal_eval successful: {sample_names_list}")
except Exception as e:
    print(f"Method 1 - ast.literal_eval failed: {e}")
    
    # Method 2: Try JSON parsing
    try:
        sample_names_list = json.loads(sample_names_str.replace("'", '"'))
        print(f"Method 2 - JSON parsing successful: {sample_names_list}")
    except Exception as e:
        print(f"Method 2 - JSON parsing failed: {e}")
        
        # Method 3: Manual string parsing
        try:
            # Remove brackets and split by comma
            clean_str = sample_names_str.strip('[]')
            if clean_str:
                sample_names_list = [name.strip().strip("'").strip('"') for name in clean_str.split(',')]
                print(f"Method 3 - Manual parsing successful: {sample_names_list}")
            else:
                sample_names_list = []
        except Exception as e:
            print(f"Method 3 - Manual parsing failed: {e}")
            sample_names_list = []

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
with open("assembly_summary.json", 'w') as json_file:
    json.dump(all_assembly_stats, json_file, indent=4)

print(f"Assembly statistics extracted for {len(all_assembly_stats)} samples")
print("Output written to: assembly_summary.json")

# Create versions file
python_version = subprocess.check_output(['python', '--version'],
                                       stderr=subprocess.STDOUT,
                                       text=True).strip().replace('Python ', '')

versions_data = {
    "ONT_FLYE:BOOTSTRAP_MODE:GATHER_ASSEMBLY_STATS": {
        "python": python_version
    }
}

with open("versions.yml", "w") as f:
    yaml.dump(versions_data, f, default_flow_style=False)
EOF

    # Run the Python script
    python gather_stats.py
    """
}