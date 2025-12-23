#!/usr/bin/env python3
"""
Assemble chromosomal contigs into final sequences based on alignment data.
Uses BioPython for FASTA I/O and sequence operations.

Reads annotated_assembly.fasta and annotated_assembly_to_genome.txt to:
1. Group contigs by chromosome
2. Concatenate multiple contigs per chromosome in order
3. Write final sequences with appropriate naming
"""

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from collections import OrderedDict

def main():
    assembly_file = 'annotated_assembly.fasta'
    mapping_file = 'annotated_assembly_to_genome.txt'
    output_file = 'final_assembly.fasta'
    
    # Read assembly into memory
    print(f"Reading assembly from {assembly_file}...")
    sequences = SeqIO.to_dict(SeqIO.parse(assembly_file, "fasta"))
    print(f"✓ Loaded {len(sequences)} sequences\n")
    
    # Parse mapping file and group by chromosome
    print(f"Reading alignments from {mapping_file}...")
    chr_groups = OrderedDict()
    
    with open(mapping_file, 'r') as f:
        for line_num, line in enumerate(f, 1):
            line = line.rstrip('\n')
            if not line.strip():
                continue
            
            try:
                fields = line.split('\t')
                contig_name = fields[0]
                chromosome = fields[2]
                
            except IndexError as e:
                print(f"✗ Error parsing line {line_num}: {e}")
                continue
            
            # Group contigs by chromosome
            if chromosome not in chr_groups:
                chr_groups[chromosome] = []
            chr_groups[chromosome].append(contig_name)
    
    print(f"✓ Found {len(chr_groups)} chromosomes\n")
    
    # Process and write output
    print(f"Processing and writing to {output_file}...\n")
    written_count = 0
    
    with open(output_file, 'w') as out_handle:
        for chromosome, contig_list in chr_groups.items():
            
            # Case 1: Single contig for this chromosome
            if len(contig_list) == 1:
                contig_name = contig_list[0]
                
                if contig_name not in sequences:
                    print(f"✗ Warning: {contig_name} not found in assembly")
                    continue
                
                seq_record = sequences[contig_name]
                SeqIO.write(seq_record, out_handle, "fasta")
                
                written_count += 1
                print(f"✓ {written_count:3d}. {contig_name:40s} {len(seq_record.seq):7d} bp [single]")
            
            # Case 2: Multiple contigs for this chromosome
            else:
                concatenated_seq = Seq("")
                contig_parts = []
                
                for contig_name in contig_list:
                    if contig_name not in sequences:
                        print(f"✗ Warning: {contig_name} not found in assembly")
                        continue
                    
                    seq_record = sequences[contig_name]
                    concatenated_seq += seq_record.seq
                    contig_parts.append(contig_name)
                
                # Build new header: first_contig__contig2_contig3...
                if len(contig_parts) > 1:
                    first_contig = contig_parts[0]
                    # Extract contig names (remove chromosome prefix)
                    chr_prefix = f"{chromosome}_"
                    remaining_contigs = []
                    for contig in contig_parts[1:]:
                        if contig.startswith(chr_prefix):
                            remaining_contigs.append(contig[len(chr_prefix):])
                        else:
                            remaining_contigs.append(contig)
                    
                    new_id = f"{first_contig}_{('_').join(remaining_contigs)}"
                else:
                    new_id = contig_parts[0]
                
                # Create new SeqRecord
                new_record = SeqRecord(
                    concatenated_seq,
                    id=new_id,
                    description=""
                )
                
                SeqIO.write(new_record, out_handle, "fasta")
                
                written_count += 1
                contig_str = " + ".join(contig_parts)
                print(f"✓ {written_count:3d}. {new_id:40s} {len(concatenated_seq):7d} bp [concat: {contig_str}]")
    
    print(f"\n{'='*70}")
    print(f"Complete!")
    print(f"  Written: {written_count} sequences")
    print(f"  Output:  {output_file}")
    print(f"{'='*70}")

if __name__ == '__main__':
    main()