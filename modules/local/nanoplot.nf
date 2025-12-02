process NANOPLOT {
    tag "${sample_name}"
    publishDir "${params.outdir}/nanoplot", mode: 'copy'
    
    input:
    tuple val(sample_name), path(fastq_file)
    
    output:
    path "${sample_name}", emit: nanoplot_results
    
    script:
    """
    NanoPlot -o ${sample_name} --fastq ${fastq_file}
    """
}