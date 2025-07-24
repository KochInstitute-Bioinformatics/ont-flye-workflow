process GATHER_ASSEMBLY_STATS {
    publishDir "${params.outdir}/assembly_summary", mode: 'copy'
    
    input:
    val sample_data_list
    
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

# Parse the sample data from Nextflow
sample_data_str = '''${sample_data_list}'''
print(f"Raw sample data: {sample_data_str}")

all_assembly_stats = {}

# Parse the sample data list
# The format should be: [[sample_name, assembly_file, log_file], ...]
import ast
try:
    sample_data = ast.literal_eval(sample_data_str)
    print(f"Parsed sample data: {sample_data}")
    
    for item in sample_data:
        if len(item) >= 3:
            sample_name = item[0]
            assembly_file = item[1]
            log_file = item[2]
            
            print(f"Processing {sample_name} (files: {assembly_file}, {log_file})...")
            
            if os.path.exists(assembly_file) and os.path.exists(log_file):
                stats = extract_assembly_stats(assembly_file, log_file, sample_name)
                all_assembly_stats[sample_name] = stats
            else:
                print(f"Warning: Files not found for {sample_name}")
                print(f"  Assembly file exists: {os.path.exists(assembly_file)}")
                print(f"  Log file exists: {os.path.exists(log_file)}")
                
except Exception as e:
    print(f"Error parsing sample data: {e}")
    print("Falling back to file discovery method...")
    
    # Fallback: try to discover files and extract names from paths
    import glob
    assembly_files = glob.glob("**/assembly_info.txt", recursive=True)
    log_files = glob.glob("**/flye.log", recursive=True)
    
    print(f"Found assembly files: {assembly_files}")
    print(f"Found log files: {log_files}")
    
    for assembly_file in assembly_files:
        # Try to extract sample name from path
        path_parts = assembly_file.split('/')
        sample_name = "unknown_sample"
        for part in path_parts:
            if 'assembly' in part and part != 'assembly_info.txt':
                sample_name = part.replace('.assembly', '')
                break
        
        if sample_name == "unknown_sample":
            sample_name = f"unknown_sample_{len(all_assembly_stats) + 1}"
        
        # Find corresponding log file
        log_file = assembly_file.replace('assembly_info.txt', 'flye.log')
        if os.path.exists(log_file):
            print(f"Processing {sample_name} (files: {assembly_file}, {log_file})...")
            stats = extract_assembly_stats(assembly_file, log_file, sample_name)
            all_assembly_stats[sample_name] = stats
        else:
            print(f"Warning: No corresponding log file found for {assembly_file}")

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