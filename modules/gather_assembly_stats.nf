process GATHER_ASSEMBLY_STATS {
    publishDir "${params.outdir}/assembly_summary", mode: 'copy'
    
    input:
    path assembly_info_files, stageAs: "assembly_info_*.txt"
    path flye_log_files, stageAs: "flye_log_*.log"
    
    output:
    path "assembly_summary.json", emit: summary
    path "versions.yml", emit: versions
    
    script:
    """
    #!/usr/bin/env python3
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
    
    # Find all staged files
    assembly_info_files = sorted(glob.glob("assembly_info_*.txt"))
    flye_log_files = sorted(glob.glob("flye_log_*.log"))
    
    print(f"Found {len(assembly_info_files)} assembly_info files")
    print(f"Found {len(flye_log_files)} flye_log files")
    
    all_assembly_stats = {}
    
    # Process files in pairs (they should be in the same order due to sorting)
    for i, assembly_info_file in enumerate(assembly_info_files):
        if i < len(flye_log_files):
            flye_log_file = flye_log_files[i]
            
            # Extract sample name from the staged filename
            # assembly_info_1.txt -> sample name from index 1
            sample_index = assembly_info_file.replace('assembly_info_', '').replace('.txt', '')
            
            # We need to get the actual sample name - let's extract it from the file path in a different way
            # For now, use the index as identifier, but we'll improve this
            sample_name = f"sample_{sample_index}"
            
            print(f"Processing {sample_name} (files: {assembly_info_file}, {flye_log_file})...")
            stats = extract_assembly_stats(assembly_info_file, flye_log_file, sample_name)
            all_assembly_stats[sample_name] = stats
    
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
        "ONT_FLYE:GATHER_ASSEMBLY_STATS": {
            "python": python_version
        }
    }
    
    with open("versions.yml", "w") as f:
        yaml.dump(versions_data, f, default_flow_style=False)
    """
}