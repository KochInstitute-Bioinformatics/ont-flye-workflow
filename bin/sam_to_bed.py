#!/usr/bin/env python3
"""
Convert SAM alignment file to BED format
Usage: sam_to_bed.py <input.sam> <output.bed>
"""

import sys
import re

def sam_to_bed(sam_file, bed_file):
    """Convert SAM to BED format"""
    with open(sam_file, 'r') as sam_in:
        with open(bed_file, 'w') as bed_out:
            for line in sam_in:
                # Skip header lines
                if line.startswith('@'):
                    continue
                
                fields = line.strip().split('\t')
                if len(fields) < 11:
                    continue
                
                qname = fields[0]   # Query name
                flag = int(fields[1])
                rname = fields[2]   # Reference name
                pos = int(fields[3]) - 1  # Convert to 0-based
                mapq = fields[4]
                cigar = fields[5]
                
                # Skip unmapped reads
                if flag & 4:
                    continue
                
                # Calculate end position from CIGAR string
                cigar_ops = re.findall(r'(\d+)([MIDNSHP=X])', cigar)
                length = sum(int(n) for n, op in cigar_ops if op in 'MDN=X')
                end = pos + length
                
                # Determine strand
                strand = '-' if flag & 16 else '+'
                
                # Write BED6 format (chrom, start, end, name, score, strand)
                bed_out.write(f"{rname}\t{pos}\t{end}\t{qname}\t{mapq}\t{strand}\n")

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: sam_to_bed.py <input.sam> <output.bed>", file=sys.stderr)
        sys.exit(1)
    
    sam_file = sys.argv[1]
    bed_file = sys.argv[2]
    
    sam_to_bed(sam_file, bed_file)