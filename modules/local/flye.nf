process FLYE {
    tag "${sample_name}"
    publishDir "${params.outdir}/assemblies", mode: 'copy'
    
    input:
    tuple val(sample_name), path(fastq_file)
    
    output:
    tuple val(sample_name), path("${sample_name}.assembly"), emit: assembly_dir
    tuple val(sample_name), path("${sample_name}.assembly/assembly.fasta"), emit: assembly_fasta
    tuple val(sample_name), path("${sample_name}.assembly/assembly_info.txt"), emit: assembly_info
    tuple val(sample_name), path("${sample_name}.assembly/flye.log"), emit: flye_log
    
    script:
    """
    flye --nano-hq ${fastq_file} \\
        --threads ${task.cpus} \\
        --genome-size ${params.genome_size} \\
        --out-dir ${sample_name}.assembly
    """
}