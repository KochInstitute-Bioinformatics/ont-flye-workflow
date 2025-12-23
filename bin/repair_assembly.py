#!/usr/bin/env python3
"""
Repair and annotate assembly based on minimap2 alignment data.
Uses BioPython for FASTA I/O and sequence operations.

Reads assembly.fasta and assembly_to_genome.txt to:
1. Rename contigs to their target chromosomes
2. Reverse complement sequences where needed (FLAG=16)
3. Skip sequences not mapping to chromosomes (col3 doesn't start with 'chr')
"""

from Bio import SeqIO
from Bio.Seq import Seq
import sys

def main():
    assembly_file = 'assembly.fasta'
    mapping_file = 'assembly_to_genome.txt'
    output_file = 'annotated_assembly.fasta'
    
    # Read assembly into memory
    print(f"Reading assembly from {assembly_file}...")
    sequences = SeqIO.to_dict(SeqIO.parse(assembly_file, "fasta"))
    print(f"✓ Loaded {len(sequences)} sequences\n")
    
    # Process mapping file and write output
    print(f"Processing mapping from {mapping_file}...\n")
    processed_count = 0
    skipped_count = 0
    
    with open(output_file, 'w') as out_handle:
        with open(mapping_file, 'r') as f:
            for line_num, line in enumerate(f, 1):
                line = line.rstrip('\n')
                if not line.strip():
                    continue
                
                try:
                    fields = line.split('\t')
                    contig_name = fields[0]
                    flag = int(fields[1])
                    chromosome = fields[2]
                    
                except (IndexError, ValueError) as e:
                    print(f"✗ Error parsing line {line_num}: {e}")
                    continue
                
                # Skip if chromosome doesn't start with 'chr'
                if not chromosome.startswith('chr'):
                    skipped_count += 1
                    continue
                
                # Get sequence record
                if contig_name not in sequences:
                    print(f"✗ Warning: {contig_name} not found in assembly")
                    continue
                
                seq_record = sequences[contig_name]
                seq = seq_record.seq
                
                # Process based on FLAG
                if flag == 0:
                    # Forward strand - just rename
                    new_id = f"{chromosome}_{contig_name}"
                    new_seq = seq
                    strand = "+"
                    
                elif flag == 16:
                    # Reverse complement
                    new_id = f"{chromosome}_{contig_name}rc"
                    new_seq = seq.reverse_complement()
                    strand = "-"
                    
                else:
                    print(f"✗ Unexpected FLAG value {flag} on line {line_num}")
                    continue
                
                # Create new SeqRecord and write
                new_record = seq_record
                new_record.id = new_id
                new_record.description = ""  # Clear description
                new_record.seq = new_seq
                
                SeqIO.write(new_record, out_handle, "fasta")
                
                processed_count += 1
                print(f"✓ {processed_count:3d}. {new_id:40s} {len(new_seq):7d} bp [{strand}]")
    
    print(f"\n{'='*70}")
    print(f"Complete!")
    print(f"  Processed: {processed_count} sequences")
    print(f"  Skipped:   {skipped_count} sequences (non-chr)")
    print(f"  Output:    {output_file}")
    print(f"{'='*70}")

if __name__ == '__main__':
    main()