include { CHOPPER } from '../modules/local/chopper'
include { DOWNSAMPLE_FASTQ } from '../modules/local/downsample_fastq'
include { FLYE_PREFLIGHT } from '../modules/local/flye_preflight'
include { PARSE_PREFLIGHT_RESULTS } from '../modules/local/parse_preflight_results'
include { FILTER_ASSEMBLY_CANDIDATES } from '../modules/local/filter_assembly_candidates'
include { FLYE } from '../modules/local/flye'
include { TRANSGENE_BLAST } from '../modules/local/transgene_blast'
include { PARSE_NANOSTATS } from '../modules/local/parse_nanostats'
include { NANOPLOT } from '../modules/local/nanoplot'
include { NANOPLOT as NANOPLOT_ORIGINAL } from '../modules/local/nanoplot'
include { PARSE_TRANSGENE_BLAST } from '../modules/local/parse_transgene_blast'
include { GATHER_ASSEMBLY_STATS } from '../modules/local/gather_assembly_stats'
include { SIMPLE_RESULTS_SUMMARY } from '../modules/local/simple_results_summary'

// ========================================
// ASSEMBLY EVALUATION PROCESSES
// These must be defined BEFORE the workflow block
// ========================================

process ALIGN_ASSEMBLY_TO_GENOME {
    tag "${sample_name}"
    publishDir "${params.outdir}/assembly_evaluation/alignments", mode: 'copy'
    
    container 'bumproo/general_genomics:latest'
    
    input:
    tuple val(sample_name), path(assembly_fasta)
    path reference_genome
    
    output:
    tuple val(sample_name), path("${sample_name}.bam"), emit: alignment_bam
    tuple val(sample_name), path("${sample_name}.bam.bai"), emit: alignment_bai
    
    script:
    """
    # Align assembly to reference using minimap2
    minimap2 -ax asm5 -t ${task.cpus} \
        ${reference_genome} \
        ${assembly_fasta} \
        | samtools sort -@ ${task.cpus} -o ${sample_name}.bam -
    
    # Index the BAM file
    samtools index ${sample_name}.bam
    """
}

process CALCULATE_ALIGNMENT_STATS {
    tag "${sample_name}"
    publishDir "${params.outdir}/assembly_evaluation/alignment_stats", mode: 'copy'
    
    container 'bumproo/general_genomics:latest'
    
    input:
    tuple val(sample_name), path(alignment_bam)
    
    output:
    tuple val(sample_name), path("${sample_name}_alignment_stats.json"), emit: stats_json
    tuple val(sample_name), path("${sample_name}_coverage.txt"), emit: coverage_txt
    
    script:
    """
    # Calculate basic alignment statistics
    samtools flagstat ${alignment_bam} > ${sample_name}_flagstat.txt
    samtools stats ${alignment_bam} > ${sample_name}_samtools_stats.txt
    samtools coverage ${alignment_bam} > ${sample_name}_coverage.txt
    
    # Parse and convert to JSON
    cat > parse_stats.py <<'PYSCRIPT'
import json
import re
import sys

sample_name = sys.argv[1]

stats = {}

# Parse flagstat
with open(f'{sample_name}_flagstat.txt', 'r') as f:
    for line in f:
        if 'mapped' in line and '%' in line:
            match = re.search(r'(\\d+)\\s+\\+\\s+\\d+\\s+mapped\\s+\\(([\\d.]+)%', line)
            if match:
                stats['mapped_reads'] = int(match.group(1))
                stats['mapping_rate'] = float(match.group(2))

# Parse coverage
with open(f'{sample_name}_coverage.txt', 'r') as f:
    next(f)  # Skip header
    coverage_sum = 0
    coverage_count = 0
    for line in f:
        parts = line.strip().split('\\t')
        if len(parts) >= 7:
            try:
                coverage_sum += float(parts[6])  # meandepth column
                coverage_count += 1
            except ValueError:
                pass
    if coverage_count > 0:
        stats['mean_coverage'] = coverage_sum / coverage_count

# Write JSON
with open(f'{sample_name}_alignment_stats.json', 'w') as f:
    json.dump({
        'sample': sample_name,
        'stats': stats
    }, f, indent=2)
PYSCRIPT

    python3 parse_stats.py "${sample_name}"
    """
}

process CALCULATE_ASSEMBLY_CONTIGUITY {
    tag "${sample_name}"
    publishDir "${params.outdir}/assembly_evaluation/contiguity", mode: 'copy'
    
    container 'bumproo/general_genomics:latest'
    
    input:
    tuple val(sample_name), path(assembly_fasta)
    
    output:
    tuple val(sample_name), path("${sample_name}_contiguity.json"), emit: contiguity_json
    
    script:
    """
    cat > calc_contiguity.py <<'PYSCRIPT'
import json
import sys
from collections import defaultdict

def parse_fasta(filename):
    # Parse FASTA file and return list of sequence lengths
    lengths = []
    current_seq = []
    
    with open(filename, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                if current_seq:
                    lengths.append(len(''.join(current_seq)))
                    current_seq = []
            else:
                current_seq.append(line)
        
        if current_seq:
            lengths.append(len(''.join(current_seq)))
    
    return lengths

def calculate_nx(lengths, x=50):
    # Calculate NX value (e.g., N50, N90)
    sorted_lengths = sorted(lengths, reverse=True)
    total_length = sum(sorted_lengths)
    target_length = total_length * (x / 100.0)
    
    cumsum = 0
    for length in sorted_lengths:
        cumsum += length
        if cumsum >= target_length:
            return length
    return 0

def calculate_lx(lengths, x=50):
    # Calculate LX value - number of contigs needed to reach X% of total length
    sorted_lengths = sorted(lengths, reverse=True)
    total_length = sum(sorted_lengths)
    target_length = total_length * (x / 100.0)
    
    cumsum = 0
    for idx, length in enumerate(sorted_lengths, 1):
        cumsum += length
        if cumsum >= target_length:
            return idx
    return len(sorted_lengths)

sample_name = sys.argv[1]
fasta_file = sys.argv[2]

# Parse assembly
lengths = parse_fasta(fasta_file)

# Calculate statistics
stats = {
    'sample': sample_name,
    'num_contigs': len(lengths),
    'total_length': sum(lengths),
    'mean_length': sum(lengths) / len(lengths) if lengths else 0,
    'min_length': min(lengths) if lengths else 0,
    'max_length': max(lengths) if lengths else 0,
    'n50': calculate_nx(lengths, 50),
    'n90': calculate_nx(lengths, 90),
    'l50': calculate_lx(lengths, 50),
    'l90': calculate_lx(lengths, 90)
}

# Write JSON
with open(f'{sample_name}_contiguity.json', 'w') as f:
    json.dump(stats, f, indent=2)
PYSCRIPT

    python3 calc_contiguity.py "${sample_name}" "${assembly_fasta}"
    """
}

process IDENTIFY_STRUCTURAL_VARIANTS {
    tag "${sample_name}"
    publishDir "${params.outdir}/assembly_evaluation/structural_variants", mode: 'copy'
    
    container 'quay.io/biocontainers/sniffles:2.4--pyhdfd78af_0'
    
    input:
    tuple val(sample_name), path(alignment_bam), path(alignment_bai)
    
    output:
    tuple val(sample_name), path("${sample_name}_sv.vcf"), emit: sv_vcf
    
    script:
    """
    # BAM index file is already provided as input
    sniffles --input ${alignment_bam} \
        --vcf ${sample_name}_sv.vcf \
        --threads ${task.cpus}
    """
}

process ALIGN_TRANSCRIPTS_TO_ASSEMBLY {
    tag "${sample_name}"
    publishDir "${params.outdir}/assembly_evaluation/transcript_alignments", mode: 'copy'
    
    container 'bumproo/general_genomics:latest'
    
    input:
    tuple val(sample_name), path(assembly_fasta)
    path transcripts_fasta
    
    output:
    tuple val(sample_name), path("${sample_name}_transcripts.bam"), emit: transcript_bam
    tuple val(sample_name), path("${sample_name}_transcript_stats.json"), emit: transcript_stats
    
    script:
    """
    # Align transcripts to assembly
    minimap2 -ax splice -t ${task.cpus} \
        ${assembly_fasta} \
        ${transcripts_fasta} \
        | samtools sort -@ ${task.cpus} -o ${sample_name}_transcripts.bam -
    
    # Calculate alignment statistics
    samtools index ${sample_name}_transcripts.bam
    samtools flagstat ${sample_name}_transcripts.bam > ${sample_name}_transcript_flagstat.txt
    
    # Parse stats to JSON
    cat > parse_transcript_stats.py <<'PYSCRIPT'
import json
import re
import sys

sample_name = sys.argv[1]

with open(f'{sample_name}_transcript_flagstat.txt', 'r') as f:
    content = f.read()
    match = re.search(r'(\\d+)\\s+\\+\\s+\\d+\\s+mapped\\s+\\(([\\d.]+)%', content)
    if match:
        stats = {
            'sample': sample_name,
            'mapped_transcripts': int(match.group(1)),
            'transcript_mapping_rate': float(match.group(2))
        }
    else:
        stats = {'sample': sample_name, 'error': 'Could not parse stats'}

with open(f'{sample_name}_transcript_stats.json', 'w') as f:
    json.dump(stats, f, indent=2)
PYSCRIPT

    python3 parse_transcript_stats.py "${sample_name}"
    """
}

process GENERATE_EVALUATION_REPORT {
    publishDir "${params.outdir}/assembly_evaluation", mode: 'copy'
    
    container 'bumproo/general_genomics:latest'
    
    input:
    path alignment_stats_files
    path contiguity_stats_files
    path sv_vcf_files
    path transcript_stats_files
    
    output:
    path "assembly_evaluation_report.html", emit: report_html
    path "assembly_evaluation_summary.json", emit: report_json
    
    script:
    """
    # Create a simple HTML report combining all metrics
    cat > generate_report.py <<'PYSCRIPT'
import json
import glob
from datetime import datetime

# Collect all stats
alignment_stats = []
for f in glob.glob('*_alignment_stats.json'):
    with open(f) as fh:
        alignment_stats.append(json.load(fh))

contiguity_stats = []
for f in glob.glob('*_contiguity.json'):
    with open(f) as fh:
        contiguity_stats.append(json.load(fh))

sv_counts = {}
for f in glob.glob('*_sv.vcf'):
    sample = f.replace('_sv.vcf', '')
    with open(f) as fh:
        count = sum(1 for line in fh if not line.startswith('#'))
    sv_counts[sample] = count

transcript_stats = []
for f in glob.glob('*_transcript_stats.json'):
    try:
        with open(f) as fh:
            transcript_stats.append(json.load(fh))
    except:
        pass

# Generate HTML report
html = '''
<!DOCTYPE html>
<html>
<head>
    <title>Assembly Evaluation Report</title>
    <style>
        body { font-family: Arial, sans-serif; margin: 20px; }
        h1 { color: #333; }
        table { border-collapse: collapse; width: 100%; margin-top: 20px; }
        th, td { border: 1px solid #ddd; padding: 8px; text-align: left; }
        th { background-color: #4CAF50; color: white; }
        tr:nth-child(even) { background-color: #f2f2f2; }
        .metric-section { margin-top: 30px; }
    </style>
</head>
<body>
    <h1>Assembly Evaluation Report</h1>
    <p>Generated: ''' + datetime.now().strftime('%Y-%m-%d %H:%M:%S') + '''</p>
    
    <div class="metric-section">
        <h2>Contiguity Metrics</h2>
        <table>
            <tr>
                <th>Sample</th>
                <th>Contigs</th>
                <th>Total Length</th>
                <th>N50</th>
                <th>L50</th>
                <th>Max Length</th>
            </tr>
'''

for stats in contiguity_stats:
    html += f'''
            <tr>
                <td>{stats['sample']}</td>
                <td>{stats['num_contigs']}</td>
                <td>{stats['total_length']:,}</td>
                <td>{stats['n50']:,}</td>
                <td>{stats['l50']}</td>
                <td>{stats['max_length']:,}</td>
            </tr>
'''

html += '''
        </table>
    </div>
    
    <div class="metric-section">
        <h2>Alignment Statistics</h2>
        <table>
            <tr>
                <th>Sample</th>
                <th>Mapped Reads</th>
                <th>Mapping Rate (%)</th>
                <th>Mean Coverage</th>
            </tr>
'''

for stats in alignment_stats:
    s = stats.get('stats', {})
    html += f'''
            <tr>
                <td>{stats['sample']}</td>
                <td>{s.get('mapped_reads', 'N/A')}</td>
                <td>{s.get('mapping_rate', 'N/A')}</td>
                <td>{s.get('mean_coverage', 'N/A'):.2f}</td>
            </tr>
'''

html += '''
        </table>
    </div>
    
    <div class="metric-section">
        <h2>Structural Variants</h2>
        <table>
            <tr>
                <th>Sample</th>
                <th>SV Count</th>
            </tr>
'''

for sample, count in sv_counts.items():
    html += f'''
            <tr>
                <td>{sample}</td>
                <td>{count}</td>
            </tr>
'''

html += '''
        </table>
    </div>
'''

if transcript_stats:
    html += '''
    <div class="metric-section">
        <h2>Transcript Alignment</h2>
        <table>
            <tr>
                <th>Sample</th>
                <th>Mapped Transcripts</th>
                <th>Mapping Rate (%)</th>
            </tr>
'''
    for stats in transcript_stats:
        html += f'''
            <tr>
                <td>{stats['sample']}</td>
                <td>{stats.get('mapped_transcripts', 'N/A')}</td>
                <td>{stats.get('transcript_mapping_rate', 'N/A')}</td>
            </tr>
'''
    html += '''
        </table>
    </div>
'''

html += '''
</body>
</html>
'''

# Write HTML report
with open('assembly_evaluation_report.html', 'w') as f:
    f.write(html)

# Write JSON summary
summary = {
    'generated': datetime.now().isoformat(),
    'alignment_stats': alignment_stats,
    'contiguity_stats': contiguity_stats,
    'sv_counts': sv_counts,
    'transcript_stats': transcript_stats
}

with open('assembly_evaluation_summary.json', 'w') as f:
    json.dump(summary, f, indent=2)
PYSCRIPT

    python3 generate_report.py
    """
}

// ========================================
// MAIN WORKFLOW
// ========================================

workflow ONT_FLYE {

main:

// Create input channel from samples.csv with validation and per-sample parameters
if (params.samples) {
    // CSV input method with file validation and parameter parsing
    input_ch = channel
        .fromPath(params.samples, checkIfExists: true)
        .splitCsv(header: true, sep: ';')  // Use semicolon separator
        .map { row ->
            def fastq_file = file(row.fastq, checkIfExists: true)
            
            // Validate that it's a file, not a directory
            if (!fastq_file.isFile()) {
                error "ERROR: ${row.fastq} is not a file! Please check your samples.csv - each 'fastq' entry must point to a FASTQ file, not a directory."
            }
            
            // Validate file extension
            def valid_extensions = ['.fastq', '.fq', '.fastq.gz', '.fq.gz']
            def has_valid_ext = valid_extensions.any { ext -> fastq_file.name.toLowerCase().endsWith(ext) }
            if (!has_valid_ext) {
                error "ERROR: ${row.fastq} does not appear to be a FASTQ file. Valid extensions: ${valid_extensions.join(', ')}"
            }
            
            def transgene = row.containsKey('transgene') ? row.transgene : params.default_transgene
            
            // Parse size_ranges with clean JSON
            def size_ranges
            if (row.containsKey('size_ranges') && row.size_ranges && row.size_ranges.trim() != '') {
                try {
                    log.info "Parsing size_ranges for ${row.name}: ${row.size_ranges}"
                    def parsed_ranges = new groovy.json.JsonSlurper().parseText(row.size_ranges)
                    
                    // Convert to the expected format with max: null
                    size_ranges = parsed_ranges.collect { range ->
                        [min: range.min as Integer, max: null, name: range.name as String]
                    }
                    log.info "Successfully parsed size_ranges for ${row.name}: ${size_ranges}"
                } catch (Exception e) {
                    log.warn "Failed to parse size_ranges for ${row.name}: ${e.message}. Using defaults."
                    log.warn "Raw size_ranges string: '${row.size_ranges}'"
                    size_ranges = [
                        [min: 40000, max: null, name: "40k_Plus"],
                        [min: 50000, max: null, name: "50k_Plus"]
                    ]
                }
            } else {
                // Default size ranges
                size_ranges = [
                    [min: 40000, max: null, name: "40k_Plus"],
                    [min: 50000, max: null, name: "50k_Plus"]
                ]
            }
            
            // Parse downsample_rates with clean JSON
            def downsample_rates
            if (row.containsKey('downsample_rates') && row.downsample_rates && row.downsample_rates.trim() != '') {
                try {
                    log.info "Parsing downsample_rates for ${row.name}: ${row.downsample_rates}"
                    def parsed_rates = new groovy.json.JsonSlurper().parseText(row.downsample_rates)
                    
                    // Convert to list of doubles
                    downsample_rates = parsed_rates.collect { rate -> rate as Double }
                    log.info "Successfully parsed downsample_rates for ${row.name}: ${downsample_rates}"
                } catch (Exception e) {
                    log.warn "Failed to parse downsample_rates for ${row.name}: ${e.message}. Using defaults."
                    log.warn "Raw downsample_rates string: '${row.downsample_rates}'"
                    downsample_rates = [0.25, 0.5]
                }
            } else {
                // Default downsample rates
                downsample_rates = [0.25, 0.5]
            }
            
            [row.name, fastq_file, transgene, size_ranges, downsample_rates]
        }
} else {
    // Single sample input method (fallback) - use global defaults
    if (!params.input_fastq || !params.name) {
        error "For single sample mode, please specify both --input_fastq and --name"
    }
    def transgene = params.transgene ?: params.default_transgene
    def default_size_ranges = [
        [min: 40000, max: null, name: "40k_Plus"],
        [min: 50000, max: null, name: "50k_Plus"]
    ]
    def default_downsample_rates = [0.25, 0.5]
    input_ch = channel.of([params.name, file(params.input_fastq, checkIfExists: true), transgene, default_size_ranges, default_downsample_rates])
}

// ========================================
// PHASE 1: FILTERING - Size and quality filtering using CHOPPER
// ========================================

// Create combinations of samples with their specific size ranges
filter_combinations = input_ch
    .flatMap { sample_name, fastq_file, transgene_name, size_ranges, downsample_rates ->
        // Create a combination for each size range for this sample
        size_ranges.collect { size_range ->
            [sample_name, fastq_file, size_range, transgene_name, downsample_rates]
        }
    }
    .map { sample_name, fastq_file, size_range, _transgene_name, _downsample_rates ->
        [sample_name, fastq_file, size_range]
    }

// Run CHOPPER for each size range (filtering based on length and quality)
CHOPPER(filter_combinations)

// ========================================
// PHASE 2: DOWNSAMPLING - Multiple replicates of downsampling
// ========================================

// Get the downsample script
downsample_script = file("${projectDir}/bin/downsample_fastq.py", checkIfExists: true)

// Create replicate numbers channel
replicate_numbers = channel.from(1..params.replicates)

// Create downsample combinations using per-sample downsample rates
// First, create a lookup map for sample-specific downsample rates
sample_downsample_map = input_ch
    .map { sample_name, _fastq_file, _transgene_name, _size_ranges, downsample_rates ->
        [sample_name, downsample_rates]
    }

// Create downsample combinations
downsample_combinations = CHOPPER.out.filtered_reads
    .combine(sample_downsample_map)
    .filter { filtered_sample_name, _filtered_fastq, original_sample_name, _downsample_rates ->
        // Match filtered samples back to their original sample for downsample rates
        filtered_sample_name.startsWith(original_sample_name)
    }
    .flatMap { filtered_sample_name, filtered_fastq, _original_sample_name, downsample_rates ->
        // Create combinations with each downsample rate for this sample
        downsample_rates.collect { rate ->
            [filtered_sample_name, filtered_fastq, rate]
        }
    }
    .combine(replicate_numbers)
    .map { sample_name_size, filtered_fastq, fraction, replicate ->
        // Create unique sample name with downsample rate and replicate number
        def new_sample_name = "${sample_name_size}_ds${fraction}_rep${replicate}"
        [new_sample_name, filtered_fastq, fraction]
    }

// Run DOWNSAMPLE_FASTQ on the size-filtered reads with multiple replicates
DOWNSAMPLE_FASTQ(downsample_combinations, downsample_script)

// ========================================
// COLLECT ALL PROCESSED FASTQ FILES
// ========================================

// Combine all FASTQ files (filtered + downsampled) into one channel
all_processed_fastq = CHOPPER.out.filtered_reads
    .mix(DOWNSAMPLE_FASTQ.out.downsampled_reads)

// ========================================
// PHASE 3: FLYE PREFLIGHT - Run preflight on all selected FASTQ files
// ========================================

// Run FLYE_PREFLIGHT on all processed FASTQ files
FLYE_PREFLIGHT(all_processed_fastq)

// Collect all preflight logs
all_preflight_logs = FLYE_PREFLIGHT.out.preflight_logs
    .map { _sample_name, log_file -> log_file }
    .collect()

// Get the parse_preflight script from bin directory
parse_preflight_script = file("${projectDir}/bin/parse_preflight_flyelog.py", checkIfExists: true)

// Parse all preflight logs and create summary table
PARSE_PREFLIGHT_RESULTS(all_preflight_logs, parse_preflight_script)

// ========================================
// PHASE 4: COVERAGE-BASED ASSEMBLY DECISION
// ========================================

// Filter assembly candidates based on coverage criteria
FILTER_ASSEMBLY_CANDIDATES(
    PARSE_PREFLIGHT_RESULTS.out.preflight_csv,
    params.min_assembly_depth,
    params.max_assembly_depth
)

// Create channel of FASTQ files that should be assembled
// Read the candidates CSV and match with FASTQ files
assembly_candidates = FILTER_ASSEMBLY_CANDIDATES.out.candidates_csv
    .splitCsv(header: true)
    .map { row -> row.sample_name }
    .combine(all_processed_fastq)
    .filter { candidate_name, sample_name, _fastq_file ->
        candidate_name == sample_name
    }
    .map { _candidate_name, sample_name, fastq_file ->
        [sample_name, fastq_file]
    }

// ========================================
// PHASE 5: FLYE ASSEMBLY - Only on selected candidates
// ========================================

// Run FLYE assembly only on candidates that passed coverage criteria
FLYE(assembly_candidates)

// ========================================
// PHASE 6: TRANSGENE BLAST ANALYSIS - Run on assembled genomes
// ========================================

// Create transgene channel from CSV file
transgene_ch = channel
    .fromPath(params.transgene_library, checkIfExists: true)
    .splitCsv(header: true)
    .map { row ->
        [row.transgene_name, file("${params.transgene_dir}/${row.fasta_file}", checkIfExists: true)]
    }

// Combine assemblies with transgene information from input
assembly_with_transgene = FLYE.out.assembly_fasta
    .combine(input_ch.map { sample_name, _fastq_file, transgene_name, _size_ranges, _downsample_rates -> [sample_name, transgene_name] })
    .filter { assembly_sample, _assembly_fasta, input_sample, _transgene_name ->
        // Match assembly samples with their original transgene assignments
        assembly_sample.startsWith(input_sample.split('_')[0]) // Handle sample name variations
    }
    .map { assembly_sample, assembly_fasta, _input_sample, transgene_name ->
        [assembly_sample, assembly_fasta, transgene_name]
    }

// Join with transgene files
blast_input = assembly_with_transgene
    .combine(transgene_ch)
    .filter { _assembly_sample, _assembly_fasta, transgene_name, transgene_file_name, _transgene_file ->
        transgene_name == transgene_file_name
    }
    .map { assembly_sample, assembly_fasta, transgene_name, _transgene_file_name, transgene_file ->
        [assembly_sample, assembly_fasta, transgene_name, transgene_file]
    }

// Run TRANSGENE_BLAST
TRANSGENE_BLAST(blast_input)

// ========================================
// PHASE 7: PARSE TRANSGENE BLAST RESULTS
// ========================================

// Collect all BLAST result files with unique names to avoid collisions
all_blast_results = TRANSGENE_BLAST.out.blast_results
    .collectFile() { sample_name, blast_file ->
        // Use sample name to create unique filenames
        ["${sample_name}.blast.txt", blast_file.text]
    }
    .collect()

// Get the parse script from bin directory
parse_transgene_script = file("${projectDir}/bin/parse_transgene_blast.py", checkIfExists: true)

// Parse all BLAST results
PARSE_TRANSGENE_BLAST(
    all_blast_results,
    file(params.transgene_library),
    parse_transgene_script
)

// ========================================
// PHASE 8: GATHER ASSEMBLY STATISTICS
// ========================================

// Collect assembly info files and flye logs from successful assemblies
assembly_info_files = FLYE.out.assembly_info
    .map { _sample_name, info_file -> info_file }
    .collect()

flye_log_files = FLYE.out.flye_log
    .map { _sample_name, log_file -> log_file }
    .collect()

// Extract sample names from successful assemblies
assembly_sample_names = FLYE.out.assembly_fasta
    .map { sample_name, _fasta_file -> sample_name }
    .collect()

// Run GATHER_ASSEMBLY_STATS
GATHER_ASSEMBLY_STATS(
    assembly_info_files,
    flye_log_files,
    assembly_sample_names,
    file("${projectDir}/bin/gather_assembly_stats.py")  // Pass script as input
)

// ========================================
// QUALITY CONTROL AND ANALYSIS - Run on original input AND all processed files
// ========================================

// Run NANOPLOT on original input files (mark them as "original")
original_input_for_nanoplot = input_ch.map { sample_name, fastq_file, _transgene_name, _size_ranges, _downsample_rates ->
    tuple("${sample_name}_original", fastq_file)
}
NANOPLOT_ORIGINAL(original_input_for_nanoplot)

// Run NANOPLOT on all processed files
NANOPLOT(all_processed_fastq)

// Collect ALL NanoPlot results (original input + all processed files)
all_nanoplot_results = NANOPLOT_ORIGINAL.out.nanoplot_results
    .mix(NANOPLOT.out.nanoplot_results)
    .collect()

// Get the parse_nanostats script from bin directory
parse_nanostats_script = file("${projectDir}/bin/parse_nanostats.py", checkIfExists: true)

// Parse all NanoStats files and create summary table
PARSE_NANOSTATS(all_nanoplot_results, parse_nanostats_script)

// results summary (simple version)
simple_summary_script = file("${projectDir}/bin/simple_summary.py", checkIfExists: true)

// Update the SIMPLE_RESULTS_SUMMARY call
SIMPLE_RESULTS_SUMMARY(
    PARSE_NANOSTATS.out.summary_json,
    PARSE_PREFLIGHT_RESULTS.out.preflight_json,
    GATHER_ASSEMBLY_STATS.out.assembly_stats,
    PARSE_TRANSGENE_BLAST.out.json_results,
    simple_summary_script
)

// ========================================
// PHASE 9: ASSEMBLY EVALUATION (OPTIONAL)
// ========================================

// Only run assembly evaluation if enabled and reference genome is provided
if (params.run_assembly_evaluation && params.reference_genome) {
    
    log.info "Assembly evaluation enabled - will align assemblies to reference genome"
    
    // Prepare reference genome
    reference_genome = file(params.reference_genome, checkIfExists: true)
    
    // Prepare transcripts FASTA (optional)
    transcripts_fasta = params.transcripts_fasta ? 
        file(params.transcripts_fasta, checkIfExists: true) : 
        file('NO_FILE')
    
    // STEP 1: Align each assembly to the reference genome using minimap2
    ALIGN_ASSEMBLY_TO_GENOME(
        FLYE.out.assembly_fasta,
        reference_genome
    )
    
    // STEP 2: Calculate alignment statistics from BAM files
    CALCULATE_ALIGNMENT_STATS(
        ALIGN_ASSEMBLY_TO_GENOME.out.alignment_bam
    )
    
    // STEP 3: Calculate assembly contiguity metrics (N50, L50, etc.)
    CALCULATE_ASSEMBLY_CONTIGUITY(
        FLYE.out.assembly_fasta
    )
    
    IDENTIFY_STRUCTURAL_VARIANTS(
    ALIGN_ASSEMBLY_TO_GENOME.out.alignment_bam
        .join(ALIGN_ASSEMBLY_TO_GENOME.out.alignment_bai)
    )
    
    // STEP 5: Align transcripts to assembly (if transcripts provided)
    if (params.transcripts_fasta) {
        ALIGN_TRANSCRIPTS_TO_ASSEMBLY(
            FLYE.out.assembly_fasta,
            transcripts_fasta
        )
    }
    
    // STEP 6: Generate combined evaluation report
    // Collect all evaluation metrics
    // Extract only the files from the tuples (discard sample names)
    alignment_stats = CALCULATE_ALIGNMENT_STATS.out.stats_json
        .map { sample_name, file -> file }
        .collect()
    contiguity_stats = CALCULATE_ASSEMBLY_CONTIGUITY.out.contiguity_json
        .map { sample_name, file -> file }
        .collect()
    sv_vcf_files = IDENTIFY_STRUCTURAL_VARIANTS.out.sv_vcf
        .map { sample_name, file -> file }
        .collect()

    transcript_alignment_stats = params.transcripts_fasta ? 
        ALIGN_TRANSCRIPTS_TO_ASSEMBLY.out.transcript_stats
            .map { sample_name, file -> file }
            .collect() :
        channel.value(file('NO_FILE'))
    
    GENERATE_EVALUATION_REPORT(
        alignment_stats,
        contiguity_stats,
        sv_vcf_files,
        transcript_alignment_stats
    )
}

emit:
    // Emit the key outputs for FASTQ generation phase
    original_input = input_ch
    selected_fastq = all_processed_fastq
    preflight_logs = FLYE_PREFLIGHT.out.preflight_logs
    preflight_summary = PARSE_PREFLIGHT_RESULTS.out.preflight_csv
    assembly_candidates = FILTER_ASSEMBLY_CANDIDATES.out.candidates_csv
    assembly_filtered = FILTER_ASSEMBLY_CANDIDATES.out.filtered_csv
    assemblies = FLYE.out.assembly_fasta
    blast_results = TRANSGENE_BLAST.out.blast_results
    transgene_summary_json = PARSE_TRANSGENE_BLAST.out.json_results
    transgene_summary_csv = PARSE_TRANSGENE_BLAST.out.csv_results
    nanoplot_results = all_nanoplot_results
    nanostats_summary = PARSE_NANOSTATS.out.summary_csv
    assembly_summary = GATHER_ASSEMBLY_STATS.out.assembly_stats
    
    // Assembly evaluation outputs (only if enabled)
    assembly_alignments = params.run_assembly_evaluation && params.reference_genome ? 
        ALIGN_ASSEMBLY_TO_GENOME.out.alignment_bam : channel.empty()
    alignment_stats = params.run_assembly_evaluation && params.reference_genome ? 
        CALCULATE_ALIGNMENT_STATS.out.stats_json : channel.empty()
    structural_variants = params.run_assembly_evaluation && params.reference_genome ? 
        IDENTIFY_STRUCTURAL_VARIANTS.out.sv_vcf : channel.empty()
    evaluation_report = params.run_assembly_evaluation && params.reference_genome ? 
        GENERATE_EVALUATION_REPORT.out.report_html : channel.empty()
}
