process PARSE_NANOSTATS {
    publishDir "${params.outdir}/summary", mode: 'copy'
    
    input:
    path nanostats_files
    
    output:
    path "nanostats_summary.tsv", emit: summary_table
    
    script:
    """
    #!/usr/bin/env python3
    import os
    import glob
    
    def parse_nanostats(file_path):
        sample = os.path.basename(os.path.dirname(file_path))
        read_count = mean_length = n50 = None
        
        with open(file_path, 'r') as f:
            for line in f:
                if line.startswith("Mean read length:"):
                    mean_length = line.split(":")[1].strip().replace(",", "")
                elif line.startswith("Number of reads:"):
                    read_count = line.split(":")[1].strip().replace(",", "")
                elif line.startswith("Read length N50:"):
                    n50 = line.split(":")[1].strip().replace(",", "")
        
        return sample, read_count, mean_length, n50
    
    # Find all NanoStats.txt files
    nanostats_files = glob.glob("*/NanoStats.txt")
    
    # Parse all files and collect results
    results = []
    for file_path in nanostats_files:
        sample, read_count, mean_length, n50 = parse_nanostats(file_path)
        results.append((sample, read_count, mean_length, n50))
    
    # Sort results by sample name for consistent output
    results.sort(key=lambda x: x[0])
    
    # Write summary table
    with open("nanostats_summary.tsv", "w") as f:
        f.write("Sample\\tReadCount\\tMeanLength\\tN50\\n")
        for sample, read_count, mean_length, n50 in results:
            f.write(f"{sample}\\t{read_count}\\t{mean_length}\\t{n50}\\n")
    """
}