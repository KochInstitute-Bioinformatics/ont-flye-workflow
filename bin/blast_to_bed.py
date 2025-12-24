#!/usr/bin/env python3

import sys

def blast_to_bed(blast_output_file):
    with open(blast_output_file, 'r') as blast_file:
        for line in blast_file:
            fields = line.strip().split('\t')
            if len(fields) < 12:
                continue
            
            query_name = fields[0]
            subject_name = fields[1]
            percent_identity = float(fields[2])
            alignment_length = int(fields[3])
            query_start = int(fields[6])
            query_end = int(fields[7])
            subject_start = int(fields[8])
            subject_end = int(fields[9])
            e_value = float(fields[10])
            bit_score = float(fields[11])
            
            if subject_start > subject_end:
                subject_start, subject_end = subject_end, subject_start
                strand = '-'
            else:
                strand = '+'
            
            bed_line = f"{subject_name}\t{subject_start - 1}\t{subject_end}\t{query_name}\t0\t{strand}\n"
            print(bed_line, end='')

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python script.py blast_output_file")
        sys.exit(1)
    
    blast_output_file = sys.argv[1]
    blast_to_bed(blast_output_file)
