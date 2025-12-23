// Assembly Evaluation Module
// This module contains processes for evaluating and refining genome assemblies
// through iterative alignment, repair, and finalization steps

// Process 1: Align assembly to reference genome
process ALIGN_ASSEMBLY_TO_GENOME {
    tag "${sample_name}"
    publishDir "${params.outdir}/assembly_evaluation/${sample_name}/alignments", mode: 'copy'
    
    input:
    tuple val(sample_name), path(assembly_fasta), path(reference_genome)
    
    output:
    tuple val(sample_name), path(assembly_fasta), path("${sample_name}_assembly_to_genome.txt"), emit: alignment
    path "versions.yml", emit: versions
    
    script:
    """
    # Align assembly to WT reference and filter
    minimap2 -a -x asm5 ${reference_genome} ${assembly_fasta} \\
        | samtools view -F 2048 -F 256 \\
        | grep -v '^@' \\
        | cut -f 1-5,12,14 \\
        | sort -k3 > ${sample_name}_assembly_to_genome.txt
    
    # Create versions file
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        minimap2: \$(minimap2 --version 2>&1)
        samtools: \$(samtools --version 2>&1 | head -n1 | sed 's/samtools //')
    END_VERSIONS
    """
}

// Process 2: Repair and annotate assembly based on alignment
process REPAIR_ASSEMBLY {
    tag "${sample_name}"
    publishDir "${params.outdir}/assembly_evaluation/${sample_name}/repaired", mode: 'copy'
    
    input:
    tuple val(sample_name), path(assembly_fasta), path(alignment_txt)
    
    output:
    tuple val(sample_name), path("${sample_name}_annotated_assembly.fasta"), emit: annotated_assembly
    path "${sample_name}_repair_assembly.log", emit: log
    path "versions.yml", emit: versions
    
    script:
    """
    # Create working copies with standard names for the Python script
    cp ${assembly_fasta} assembly.fasta
    cp ${alignment_txt} assembly_to_genome.txt
    
    # Run repair script (outputs to annotated_assembly.fasta)
    repair_assembly.py > ${sample_name}_repair_assembly.log
    
    # Rename output to include sample name
    mv annotated_assembly.fasta ${sample_name}_annotated_assembly.fasta
    
    # Create versions file
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version 2>&1 | sed 's/Python //')
        biopython: \$(python -c "import Bio; print(Bio.__version__)")
    END_VERSIONS
    """
}

// Process 3: Align annotated assembly to reference
process ALIGN_ANNOTATED_ASSEMBLY {
    tag "${sample_name}"
    publishDir "${params.outdir}/assembly_evaluation/${sample_name}/alignments", mode: 'copy'
    
    input:
    tuple val(sample_name), path(annotated_assembly), path(reference_genome)
    
    output:
    tuple val(sample_name), path(annotated_assembly), path("${sample_name}_annotated_assembly_to_genome.txt"), emit: alignment
    path "versions.yml", emit: versions
    
    script:
    """
    # Align annotated assembly to reference
    minimap2 -a -x asm5 ${reference_genome} ${annotated_assembly} \\
        | samtools view -F 2048 -F 256 \\
        | grep -v '^@' \\
        | cut -f 1-5,12,14 \\
        | sort -k3,3 -k4,4n > ${sample_name}_annotated_assembly_to_genome.txt
    
    # Create versions file
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        minimap2: \$(minimap2 --version 2>&1)
        samtools: \$(samtools --version 2>&1 | head -n1 | sed 's/samtools //')
    END_VERSIONS
    """
}

// Process 4: Finalize assembly by concatenating chromosomal contigs
process FINALIZE_ASSEMBLY {
    tag "${sample_name}"
    publishDir "${params.outdir}/assembly_evaluation/${sample_name}/final", mode: 'copy'
    
    input:
    tuple val(sample_name), path(annotated_assembly), path(alignment_txt)
    
    output:
    tuple val(sample_name), path("${sample_name}_final_assembly.fasta"), emit: final_assembly
    path "${sample_name}_finalize_assembly.log", emit: log
    path "versions.yml", emit: versions
    
    script:
    """
    # Create working copies with standard names for the Python script
    cp ${annotated_assembly} annotated_assembly.fasta
    cp ${alignment_txt} annotated_assembly_to_genome.txt
    
    # Run finalize script (outputs to final_assembly.fasta)
    final_assembly.py > ${sample_name}_finalize_assembly.log
    
    # Rename output to include sample name
    mv final_assembly.fasta ${sample_name}_final_assembly.fasta
    
    # Create versions file
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version 2>&1 | sed 's/Python //')
        biopython: \$(python -c "import Bio; print(Bio.__version__)")
    END_VERSIONS
    """
}

// Process 5: Align final assembly to reference
process ALIGN_FINAL_ASSEMBLY {
    tag "${sample_name}"
    publishDir "${params.outdir}/assembly_evaluation/${sample_name}/alignments", mode: 'copy'
    
    input:
    tuple val(sample_name), path(final_assembly), path(reference_genome)
    
    output:
    tuple val(sample_name), path(final_assembly), path("${sample_name}_final_assembly_to_genome.txt"), emit: alignment
    path "versions.yml", emit: versions
    
    script:
    """
    # Align final assembly to reference
    minimap2 -a -x asm5 ${reference_genome} ${final_assembly} \\
        | samtools view -F 2048 -F 256 \\
        | grep -v '^@' \\
        | cut -f 1-5,12,14 \\
        | sort -k3,3 -k4,4n > ${sample_name}_final_assembly_to_genome.txt
    
    # Create versions file
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        minimap2: \$(minimap2 --version 2>&1)
        samtools: \$(samtools --version 2>&1 | head -n1 | sed 's/samtools //')
    END_VERSIONS
    """
}

// Process 6: Map ONT reads to final assembly
process MAP_READS_TO_ASSEMBLY {
    tag "${sample_name}"
    publishDir "${params.outdir}/assembly_evaluation/${sample_name}/read_mapping", mode: 'copy'
    
    input:
    tuple val(sample_name), path(final_assembly), path(query_fastq)
    
    output:
    tuple val(sample_name), path(final_assembly), path("${sample_name}_ont_to_assembled.sorted.bam"), emit: mapped_reads
    path "${sample_name}_ont_to_assembled.sorted.bam.bai", emit: bam_index
    path "versions.yml", emit: versions
    
    script:
    """
    # Map ONT reads to final assembly
    minimap2 -ax map-ont ${final_assembly} ${query_fastq} > ont_to_assembled.sam
    
    # Convert to BAM, sort, and index
    samtools view -b ont_to_assembled.sam -o ont_to_assembled.bam
    samtools sort -@ ${task.cpus} ont_to_assembled.bam -o ${sample_name}_ont_to_assembled.sorted.bam
    samtools index ${sample_name}_ont_to_assembled.sorted.bam
    
    # Clean up intermediate files
    rm ont_to_assembled.sam ont_to_assembled.bam
    
    # Create versions file
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        minimap2: \$(minimap2 --version 2>&1)
        samtools: \$(samtools --version 2>&1 | head -n1 | sed 's/samtools //')
    END_VERSIONS
    """
}

// Process 7: BLAST transgene against final assembly
process BLAST_TRANSGENE_TO_ASSEMBLY {
    tag "${sample_name}_${transgene_name}"
    publishDir "${params.outdir}/assembly_evaluation/${sample_name}/transgene_blast", mode: 'copy'
    
    input:
    tuple val(sample_name), path(final_assembly), val(transgene_name), path(transgene_fasta)
    
    output:
    tuple val(sample_name), path("${sample_name}_${transgene_name}_transgene_blast.txt"), emit: blast_results
    path "versions.yml", emit: versions
    
    script:
    """
    # Create BLAST database from final assembly
    makeblastdb -in ${final_assembly} -dbtype nucl -out final_assembly_db
    
    # Run BLAST search
    blastn -query ${transgene_fasta} \\
        -db final_assembly_db \\
        -outfmt 6 \\
        > ${sample_name}_${transgene_name}_transgene_blast.txt
    
    # Create versions file
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        blast: \$(blastn -version 2>&1 | head -n1 | sed 's/blastn: //')
    END_VERSIONS
    """
}

// Process 8: Convert BLAST results to BED format
process CONVERT_BLAST_TO_BED {
    tag "${sample_name}_${transgene_name}"
    publishDir "${params.outdir}/assembly_evaluation/${sample_name}/transgene_blast", mode: 'copy'
    
    input:
    tuple val(sample_name), path(blast_txt), val(transgene_name)
    
    output:
    tuple val(sample_name), path("${sample_name}_${transgene_name}_transgene_blast.bed"), emit: bed_file
    path "versions.yml", emit: versions
    
    script:
    """
    # Convert BLAST output to BED format
    blast_to_bed.py ${blast_txt} > ${sample_name}_${transgene_name}_transgene_blast.bed
    
    # Create versions file
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version 2>&1 | sed 's/Python //')
    END_VERSIONS
    """
}

// Process 9: Map transcripts to final assembly (optional)
process MAP_TRANSCRIPTS_TO_ASSEMBLY {
    tag "${sample_name}"
    publishDir "${params.outdir}/assembly_evaluation/${sample_name}/transcript_mapping", mode: 'copy'
    
    input:
    tuple val(sample_name), path(final_assembly), path(transcripts_fasta)
    
    output:
    tuple val(sample_name), path("${sample_name}_transcripts_to_assembly.bed"), emit: transcript_bed
    path "${sample_name}_transcripts_to_assembly.sam", emit: transcript_sam
    path "versions.yml", emit: versions
    
    script:
    """
    # Map transcripts to assembly
    minimap2 -a ${final_assembly} ${transcripts_fasta} > ${sample_name}_transcripts_to_assembly.sam
    
    # Convert SAM to BED
    sam2bed < ${sample_name}_transcripts_to_assembly.sam > ${sample_name}_transcripts_to_assembly.bed
    
    # Create versions file
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        minimap2: \$(minimap2 --version 2>&1)
        bedops: \$(sam2bed --version 2>&1 || echo "2.4.41")
    END_VERSIONS
    """
}