process PARSE_NANOSTATS {
    publishDir "${params.outdir}/nanostats_summary", mode: 'copy'
    
    input:
    path nanostats_files
    
    output:
    path "nanostats_summary.json", emit: summary
    path "versions.yml", emit: versions
    
    script:
    """
    #!/usr/bin/env python3
    import os
    import glob
    import json

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

        return {
            "Sample": sample,
            "ReadCount": read_count,
            "MeanLength": mean_length,
            "N50": n50
        }

    # Find all NanoStats.txt files
    nanostats_files = glob.glob("*/NanoStats.txt")

    # Parse all files and collect results
    results = [parse_nanostats(file_path) for file_path in nanostats_files]

    # Sort results by sample name
    results.sort(key=lambda x: x["Sample"])

    # Write summary to JSON file
    with open("nanostats_summary.json", "w") as f:
        json.dump(results, f, indent=4)
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}