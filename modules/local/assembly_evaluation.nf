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
    tuple val(sample_name), path('input_assembly.fasta'), path('input_alignment.txt')
    
    output:
    tuple val(sample_name), path("${sample_name}_annotated_assembly.fasta"), emit: annotated_assembly
    path "${sample_name}_repair_assembly.log", emit: log
    path "versions.yml", emit: versions
    
    script:
    """
    # Create working copies with standard names for the Python script
    cp input_assembly.fasta assembly.fasta
    cp input_alignment.txt assembly_to_genome.txt
    
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

// Process 3a: Calculate chromosome coverage for annotated assembly
process CALCULATE_CHROMOSOME_COVERAGE {
    tag "${sample_name}"
    publishDir "${params.outdir}/assembly_evaluation/${sample_name}/coverage", mode: 'copy'
    
    input:
    tuple val(sample_name), path(annotated_assembly), path(reference_chromosomes)
    
    output:
    tuple val(sample_name), path("${sample_name}_chromosome_coverage.csv"), emit: coverage
    path "versions.yml", emit: versions
    
    script:
    """
    # Align annotated assembly to chromosome reference
    minimap2 -a -x asm5 ${reference_chromosomes} ${annotated_assembly} \\
        | samtools view -F 2048 -F 256 -b \\
        | samtools sort -o temp.bam -

    # Index the BAM file
    samtools index temp.bam

    # Create contig-to-chromosome mapping with alignment counts
    # Step 1: Extract chromosome and contig name from alignments, count occurrences
    samtools view temp.bam | awk '{print \$3 "\\t" \$1}' | \\
        sort | uniq -c | \\
        awk '{print \$2 "\\t" \$3 "\\t" \$1}' > contig_chr_count.txt

    # Get coverage data: chromosome, length, meandepth
    samtools coverage temp.bam | tail -n +2 | cut -f1,3,7 | sort -k1,1 > coverage_data.txt

    # ========================================
    # COVERAGE REPORT (all contigs with alignment counts)
    # ========================================
    echo "sample,chromosome,contig,length,coverage,alignment_count" > ${sample_name}_chromosome_coverage.csv

    # For each line in contig_chr_count.txt, look up the chromosome coverage
    # contig_chr_count.txt: chr contig count
    # coverage_data.txt: chr length meandepth
    while IFS=\$'\\t' read -r chr contig count; do
        # Look up coverage for this chromosome
        cov_line=\$(grep "^\${chr}"\$'\\t' coverage_data.txt)
        if [ -n "\$cov_line" ]; then
            length=\$(echo "\$cov_line" | cut -f2)
            coverage=\$(echo "\$cov_line" | cut -f3)
            echo "${sample_name},\${chr},\${contig},\${length},\${coverage},\${count}"
        fi
    done < contig_chr_count.txt >> ${sample_name}_chromosome_coverage.csv

    # Clean up
    rm -f temp.bam temp.bam.bai contig_chr_count.txt coverage_data.txt

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
    tuple val(sample_name), path('input_annotated_assembly.fasta'), path('input_alignment.txt')
    
    output:
    tuple val(sample_name), path("${sample_name}_final_assembly.fasta"), emit: final_assembly
    path "${sample_name}_finalize_assembly.log", emit: log
    path "versions.yml", emit: versions
    
    script:
    """
    # Create working copies with standard names for the Python script
    cp input_annotated_assembly.fasta annotated_assembly.fasta
    cp input_alignment.txt annotated_assembly_to_genome.txt
    
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
    
    input:
    tuple val(sample_name), path(final_assembly), path(query_fastq)
    
    output:
    tuple val(sample_name), path(final_assembly), path("${sample_name}_ont_to_assembled.sorted.bam"), emit: mapped_reads
    tuple val(sample_name), path("${sample_name}_ont_to_assembled.sorted.bam.bai"), emit: bam_index
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
    
    input:
    tuple val(sample_name), path(final_assembly), val(transgene_name), path(transgene_file)
    
    output:
    tuple val(sample_name), path("${sample_name}_transgene_blast.txt"), path(final_assembly), emit: blast_results
    path "versions.yml", emit: versions
    
    script:
    """
    # Create BLAST database from final assembly
    makeblastdb -in ${final_assembly} -dbtype nucl -out assembly_db
    
    # BLAST transgene against final assembly
    blastn -query ${transgene_file} -db assembly_db \\
        -outfmt 6 -out ${sample_name}_transgene_blast.txt
    
    # Create versions file
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        blast: \$(blastn -version 2>&1 | grep "blastn")
    END_VERSIONS
    """
}

// Workflow to handle transgene BLAST processing
workflow ASSEMBLY_EVALUATION_WORKFLOW {
    take:
    final_assembly_input
    
    main:
    // Read sample info from CSV if available
    if (params.samples) {
        sample_info = channel
            .fromPath(params.samples)
            .splitCsv(header: true)
            .map { row ->
                tuple(row.sample, row.transgene ?: params.default_transgene)
            }
    } else {
        // Use default for all
        sample_info = channel.value([null, params.default_transgene])
    }
    
    final_with_transgene = final_assembly_input
        .map { sample_name, final_fasta ->
            // Extract base sample (e.g., S-1077-1 from S-1077-1_70k_Plus_ds0.5_rep1)
            def base_sample = sample_name.replaceAll(/_\d+k_Plus.*/, '')
            tuple(base_sample, sample_name, final_fasta)
        }
        .combine(sample_info)
        .filter { base_sample, _sample_name, _final_fasta, csv_sample, _transgene ->
            csv_sample == null || csv_sample == base_sample
        }
        .map { _base_sample, sample_name, final_fasta, _csv_sample, transgene ->
            def transgene_file = file("${params.transgene_dir}/${transgene}.fa", checkIfExists: true)
            tuple(sample_name, final_fasta, transgene, transgene_file)
        }
    
    BLAST_TRANSGENE_TO_ASSEMBLY(final_with_transgene)
    
    emit:
    blast_results = BLAST_TRANSGENE_TO_ASSEMBLY.out.blast_results
}

// Process 8: Convert BLAST results to BED format
process CONVERT_BLAST_TO_BED {
    tag "${sample_name}"
    
    input:
    tuple val(sample_name), path(blast_txt)
    
    output:
    tuple val(sample_name), path("${sample_name}_transgene_blast.bed"), emit: bed_file
    path "versions.yml", emit: versions
    
    script:
    """
    # Convert BLAST output to BED format
    blast_to_bed.py ${blast_txt} > ${sample_name}_transgene_blast.bed
    
    # Create versions file
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version 2>&1 | sed 's/Python //')
    END_VERSIONS
    """
}

// Process 9: Map transcripts to final assembly
process MAP_TRANSCRIPTS_TO_ASSEMBLY {
    tag "${sample_name}"
    
    input:
    tuple val(sample_name), path(final_assembly)
    path transcripts_fasta
    
    output:
    tuple val(sample_name), path("${sample_name}_transcripts_to_assembly.bed"), emit: transcript_bed
    path "versions.yml", emit: versions
    
    script:
    """
    # Map transcripts to assembly
    minimap2 -a ${final_assembly} ${transcripts_fasta} > transcripts_to_assembly.sam
    
    # Convert SAM to BED using a local script
    
    sam_to_bed.py transcripts_to_assembly.sam ${sample_name}_transcripts_to_assembly.bed

    # Clean up SAM file
    rm transcripts_to_assembly.sam
    
    # Create versions file
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        minimap2: \$(minimap2 --version 2>&1)
        bedops: \$(sam2bed --version 2>&1 || echo "2.4.35")
    END_VERSIONS
    """
}

// Process 10a: Consolidate PRIMARY coverage summary from all samples
process CONSOLIDATE_COVERAGE_SUMMARY {
    publishDir "${params.outdir}/summary", mode: 'copy'
    
    input:
    path coverage_files
    
    output:
    path "coverage_summary_primary.csv", emit: coverage_summary
    
    script:
    """
    # Create header with contig column
    echo "sample,chromosome,contig,length,coverage,is_primary" > coverage_summary_primary.csv
    
    # Debug: List all input files
    echo "Processing PRIMARY coverage files:" >&2
    ls -lh ${coverage_files} >&2
    echo "" >&2
    
    # Concatenate all coverage files (skip headers)
    for file in ${coverage_files}; do
        echo "Processing: \${file}" >&2
        if [ -f "\${file}" ] && [ -s "\${file}" ]; then
            echo "  File exists and has content" >&2
            tail -n +2 "\${file}" >> coverage_summary_primary.csv
        else
            echo "  WARNING: File is empty or doesn't exist!" >&2
        fi
    done
    
    # Report final line count
    LINES=\$(wc -l < coverage_summary_primary.csv)
    echo "" >&2
    echo "Final coverage_summary_primary.csv has \${LINES} lines (including header)" >&2
    """
}

// Process 10b: Consolidate DETAILED coverage summary from all samples
process CONSOLIDATE_DETAILED_COVERAGE_SUMMARY {
    publishDir "${params.outdir}/summary", mode: 'copy'
    
    input:
    path coverage_files
    
    output:
    path "coverage_summary_detailed.csv", emit: coverage_summary_detailed
    
    script:
    """
    # Create header with all columns
    echo "sample,chromosome,contig,length,coverage,alignment_count,is_primary" > coverage_summary_detailed.csv
    
    # Debug: List all input files
    echo "Processing DETAILED coverage files:" >&2
    ls -lh ${coverage_files} >&2
    echo "" >&2
    
    # Concatenate all coverage files (skip headers)
    for file in ${coverage_files}; do
        echo "Processing: \${file}" >&2
        if [ -f "\${file}" ] && [ -s "\${file}" ]; then
            echo "  File exists and has content" >&2
            tail -n +2 "\${file}" >> coverage_summary_detailed.csv
        else
            echo "  WARNING: File is empty or doesn't exist!" >&2
        fi
    done
    
    # Report final line count
    LINES=\$(wc -l < coverage_summary_detailed.csv)
    echo "" >&2
    echo "Final coverage_summary_detailed.csv has \${LINES} lines (including header)" >&2
    """
}

// Process 10: Consolidate IGV data into strain-specific directories
process CONSOLIDATE_IGV_DATA {
    tag "${sample_name}"
    publishDir "${params.outdir}/igv_data/${base_strain}", mode: 'copy'
    
    input:
    tuple val(base_strain),
          val(sample_name), 
          path(final_assembly), 
          path(bam_file), 
          path(bam_index),
          path(transgene_bed),
          path(transcript_bed, stageAs: 'transcript_optional.bed')
    
    output:
    tuple val(base_strain), 
          path("${sample_name}_final_assembly.fasta"),
          path("${sample_name}_ont_to_assembled.sorted.bam"),
          path("${sample_name}_ont_to_assembled.sorted.bam.bai"),
          path("${sample_name}_transgene_blast.bed"),
          path("${sample_name}_transcripts_to_assembly.bed", optional: true),
          emit: igv_files
    
    script:
    """
    # Nextflow stages files as symlinks, but we need actual files for outputs
    # Copy files to ensure they exist as real files in the work directory
    
    # Handle final assembly - always create new file for output matching
    TARGET="${sample_name}_final_assembly.fasta"
    TEMP=\$(mktemp)
    cp -L ${final_assembly} \${TEMP}
    mv \${TEMP} \${TARGET}
    
    # Handle BAM file - always create new file for output matching
    TARGET="${sample_name}_ont_to_assembled.sorted.bam"
    TEMP=\$(mktemp)
    cp -L ${bam_file} \${TEMP}
    mv \${TEMP} \${TARGET}
    
    # Handle BAM index - always create new file for output matching
    TARGET="${sample_name}_ont_to_assembled.sorted.bam.bai"
    TEMP=\$(mktemp)
    cp -L ${bam_index} \${TEMP}
    mv \${TEMP} \${TARGET}
    
    # Handle transgene BED - copy with standardized filename
    TARGET="${sample_name}_transgene_blast.bed"
    TEMP=\$(mktemp)
    cp -L ${transgene_bed} \${TEMP}
    mv \${TEMP} \${TARGET}
    
    # Handle optional transcript file - always create new file for output matching
    if [ -f "${transcript_bed}" ] && [ -s "${transcript_bed}" ]; then
        TARGET="${sample_name}_transcripts_to_assembly.bed"
        TEMP=\$(mktemp)
        cp -L ${transcript_bed} \${TEMP}
        mv \${TEMP} \${TARGET}
        echo "✓ Transcript BED included: ${transcript_bed}"
    else
        echo "○ No transcript BED file provided (optional)"
    fi
    
    # Set proper permissions: rw-rw-r-- (664) for owner+group write, world read
    chmod 664 ${sample_name}_final_assembly.fasta || true
    chmod 664 ${sample_name}_ont_to_assembled.sorted.bam || true
    chmod 664 ${sample_name}_ont_to_assembled.sorted.bam.bai || true
    chmod 664 ${sample_name}_transgene_blast.bed || true
    if [ -f "${sample_name}_transcripts_to_assembly.bed" ]; then
        chmod 664 ${sample_name}_transcripts_to_assembly.bed || true
    fi
    
    # Log what we're consolidating
    echo ""
    echo "═══════════════════════════════════════"
    echo "Consolidating IGV data for ${sample_name}"
    echo "  Base strain: ${base_strain}"
    echo "  Output directory: igv_data/${base_strain}/"
    echo "═══════════════════════════════════════"
    echo ""
    echo "Files prepared:"
    ls -lh ${sample_name}* 2>/dev/null || echo "  (listing files...)"
    echo ""
    """
}