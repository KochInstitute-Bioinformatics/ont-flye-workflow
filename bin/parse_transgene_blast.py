#!/usr/bin/env python3
"""
Parse transgene BLAST results to count full-length alignments.

This script processes BLAST output files and transgene library information
to identify full-length transgene alignments in assemblies.
"""

import argparse
import csv
import json
import os
import sys
from pathlib import Path
from collections import defaultdict


def load_transgene_library(library_file):
    """Load transgene library with lengths."""
    transgenes = {}
    with open(library_file, 'r') as f:
        reader = csv.DictReader(f)
        for row in reader:
            transgenes[row['transgene_name']] = {
                'length': int(row['length']),
                'description': row['description']
            }
    return transgenes


def parse_blast_file(blast_file, transgene_length, tolerance=100):
    """
    Parse a single BLAST output file and identify full-length alignments.
    
    BLAST tabular format columns:
    0: qseqid, 1: sseqid, 2: pident, 3: length, 4: mismatch, 5: gapopen,
    6: qstart, 7: qend, 8: sstart, 9: send, 10: evalue, 11: bitscore
    """
    full_length_alignments = []
    
    if not os.path.exists(blast_file) or os.path.getsize(blast_file) == 0:
        return full_length_alignments
    
    with open(blast_file, 'r') as f:
        for line in f:
            if line.strip() and not line.startswith('#'):
                fields = line.strip().split('\t')
                if len(fields) >= 12:
                    alignment_length = int(fields[3])  # Column 3: alignment length
                    contig_name = fields[1]  # Column 1: subject sequence ID (contig)
                    
                    # Check if alignment is full-length (within tolerance)
                    length_diff = abs(transgene_length - alignment_length)
                    if length_diff <= tolerance:
                        full_length_alignments.append({
                            'contig': contig_name,
                            'alignment_length': alignment_length,
                            'length_diff': length_diff,
                            'identity': float(fields[2]),
                            'evalue': float(fields[10])
                        })
    
    return full_length_alignments


def process_blast_results(blast_files, transgene_library, tolerance=100):
    """Process all BLAST result files."""
    results = []
    
    for blast_file in blast_files:
        blast_path = Path(blast_file)
        
        # Extract sample and transgene names from filename
        # Expected format: {sample_name}_{transgene_name}_blast.txt
        filename_parts = blast_path.stem.replace('_blast', '').split('_')
        
        if len(filename_parts) < 2:
            print(f"Warning: Could not parse filename {blast_path.name}", file=sys.stderr)
            continue
        
        # Find transgene name by checking against library
        transgene_name = None
        sample_name = None
        
        for tg_name in transgene_library.keys():
            if tg_name in blast_path.stem:
                transgene_name = tg_name
                sample_name = blast_path.stem.replace(f'_{tg_name}_blast', '')
                break
        
        if not transgene_name:
            print(f"Warning: Could not identify transgene for {blast_path.name}", file=sys.stderr)
            continue
        
        transgene_length = transgene_library[transgene_name]['length']
        
        # Parse BLAST results
        full_length_alignments = parse_blast_file(blast_file, transgene_length, tolerance)
        
        # Count full-length alignments and identify contigs
        full_length_count = len(full_length_alignments)
        contigs = list(set([aln['contig'] for aln in full_length_alignments]))
        contig_names = ','.join(sorted(contigs)) if contigs else 'None'
        
        results.append({
            'assembly_name': sample_name,
            'transgene_name': transgene_name,
            'transgene_length': transgene_length,
            'full_length_count': full_length_count,
            'contig_names': contig_names,
            'alignments': full_length_alignments
        })
    
    return results


def write_outputs(results, output_dir):
    """Write JSON and CSV output files."""
    os.makedirs(output_dir, exist_ok=True)
    
    # Prepare data for output
    summary_data = []
    detailed_data = {}
    
    for result in results:
        summary_data.append({
            'assembly_name': result['assembly_name'],
            'transgene_name': result['transgene_name'],
            'transgene_length': result['transgene_length'],
            'full_length_count': result['full_length_count'],
            'contig_names': result['contig_names']
        })
        
        detailed_data[f"{result['assembly_name']}_{result['transgene_name']}"] = {
            'assembly_name': result['assembly_name'],
            'transgene_name': result['transgene_name'],
            'transgene_length': result['transgene_length'],
            'full_length_count': result['full_length_count'],
            'contig_names': result['contig_names'],
            'alignments': result['alignments']
        }
    
    # Write JSON file
    json_file = os.path.join(output_dir, 'transgene_count.json')
    with open(json_file, 'w') as f:
        json.dump(detailed_data, f, indent=2)
    
    # Write CSV file
    csv_file = os.path.join(output_dir, 'transgene_count.csv')
    with open(csv_file, 'w', newline='') as f:
        if summary_data:
            fieldnames = ['assembly_name', 'transgene_name', 'transgene_length', 
                         'full_length_count', 'contig_names']
            writer = csv.DictWriter(f, fieldnames=fieldnames)
            writer.writeheader()
            writer.writerows(summary_data)
    
    print(f"Results written to:")
    print(f"  JSON: {json_file}")
    print(f"  CSV: {csv_file}")


def main():
    parser = argparse.ArgumentParser(description='Parse transgene BLAST results')
    parser.add_argument('--blast-files', nargs='+', required=True,
                       help='BLAST result files to process')
    parser.add_argument('--transgene-library', required=True,
                       help='Transgene library CSV file')
    parser.add_argument('--output-dir', required=True,
                       help='Output directory for results')
    parser.add_argument('--tolerance', type=int, default=100,
                       help='Length tolerance for full-length alignment (default: 100)')
    
    args = parser.parse_args()
    
    # Load transgene library
    try:
        transgene_library = load_transgene_library(args.transgene_library)
        print(f"Loaded {len(transgene_library)} transgenes from library")
    except Exception as e:
        print(f"Error loading transgene library: {e}", file=sys.stderr)
        sys.exit(1)
    
    # Process BLAST results
    try:
        results = process_blast_results(args.blast_files, transgene_library, args.tolerance)
        print(f"Processed {len(results)} BLAST result files")
    except Exception as e:
        print(f"Error processing BLAST results: {e}", file=sys.stderr)
        sys.exit(1)
    
    # Write outputs
    try:
        write_outputs(results, args.output_dir)
    except Exception as e:
        print(f"Error writing outputs: {e}", file=sys.stderr)
        sys.exit(1)


if __name__ == '__main__':
    main()