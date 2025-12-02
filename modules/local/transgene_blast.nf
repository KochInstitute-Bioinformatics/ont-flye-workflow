process TRANSGENE_BLAST {
    container "ncbi/blast:latest"
    
    publishDir "${params.outdir}/blast", mode: 'copy'
    
    input:
    tuple val(sample_name), path(assembly_fasta), val(transgene_name), path(transgene_fasta)
    
    output:
    tuple val(sample_name), path("${sample_name}_${transgene_name}_blast.txt"), emit: blast_results
    path "versions.yml", emit: versions
    
    script:
    def args = task.ext.args ?: ''
    """
    # Create BLAST database from assembly
    makeblastdb -in ${assembly_fasta} -dbtype nucl -out assembly_db
    
    # Run BLAST search
    blastn -query ${transgene_fasta} -db assembly_db -outfmt 6 ${args} > ${sample_name}_${transgene_name}_blast.txt
    
    # Create versions file
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        blast: \$(blastn -version 2>&1 | head -n1 | sed 's/blastn: //')
    END_VERSIONS
    """
    
    stub:
    """
    touch ${sample_name}_${transgene_name}_blast.txt
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        blast: "2.14.0+"
    END_VERSIONS
    """
}