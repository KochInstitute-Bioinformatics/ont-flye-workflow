# Assembly Evaluation Module - Implementation Summary

## ✅ Implementation Complete!

Successfully implemented the **Assembly Evaluation Module** for your ONT-Flye workflow on the `with_eval_v5` branch.

### 📁 Files Created/Modified

1. **`modules/local/assembly_evaluation.nf`** (282 lines)
   - Complete module with 9 processes implementing your evaluation pipeline

2. **`ASSEMBLY_EVALUATION.md`** (323 lines)
   - Comprehensive documentation with usage examples and integration guide

3. **`nextflow.config`** (modified)
   - Added process configurations for all 9 evaluation processes
   - Added parameters for assembly evaluation control

### 🔬 Implemented Processes

The module includes **9 sequential processes** that mirror your original bash script workflow:

| # | Process Name | Function | Tools Used |
|---|-------------|----------|------------|
| 1 | **ALIGN_ASSEMBLY_TO_GENOME** | Initial alignment to reference | minimap2 + samtools |
| 2 | **REPAIR_ASSEMBLY** | Rename & orient contigs | repair_assembly.py (BioPython) |
| 3 | **ALIGN_ANNOTATED_ASSEMBLY** | Re-align repaired assembly | minimap2 + samtools |
| 4 | **FINALIZE_ASSEMBLY** | Concatenate chromosomal contigs | final_assembly.py (BioPython) |
| 5 | **ALIGN_FINAL_ASSEMBLY** | Validate final assembly | minimap2 + samtools |
| 6 | **MAP_READS_TO_ASSEMBLY** | Map ONT reads to final assembly | minimap2 + samtools (BAM output) |
| 7 | **BLAST_TRANSGENE_TO_ASSEMBLY** | Detect transgene insertions | BLAST (makeblastdb + blastn) |
| 8 | **CONVERT_BLAST_TO_BED** | Convert BLAST to BED format | blast_to_bed.py |
| 9 | **MAP_TRANSCRIPTS_TO_ASSEMBLY** | Map transcripts (optional) | minimap2 + sam2bed |

### ✨ Key Features

✅ **Fully Containerized** - Uses biocontainers for reproducibility  
✅ **Lint-Validated** - Passes `nextflow lint` with no errors  
✅ **Published Outputs** - All intermediate and final results saved  
✅ **Version Tracking** - Each process emits software versions  
✅ **Python Scripts Integrated** - Your existing scripts in `bin/` are used  
✅ **Resource Optimized** - Sensible CPU/memory defaults configured  
✅ **Well Documented** - Complete usage guide with examples  

### 🐳 Containers Used

- **Minimap2 + Samtools + BioPython**: `quay.io/biocontainers/mulled-v2-66534bcbb7031a148b13e2ad42583020b9cd25c4:1679e915ddb9d6b4abda91880c4b48857d471bd8-0`
- **BLAST**: `ncbi/blast:latest`
- **Python**: `bumproo/python:3.11`
- **Bedops**: `quay.io/biocontainers/mulled-v2-4ec040a1b354301aafa7731db41a5f6e9ca25e70:8110a70be2bfe7f75cbd12b7a4a10c5d8a6ef337-0`

### 📊 Output Structure

```
results/assembly_evaluation/<sample_name>/
├── alignments/           # All alignment files (.txt)
├── repaired/            # Repaired assembly + logs
├── final/               # Final assembly + logs  
├── read_mapping/        # BAM files with ONT read mapping
├── transgene_blast/     # BLAST results + BED files
└── transcript_mapping/  # Transcript alignments (optional)
```

### 🚀 How to Use

**1. Enable in your config or command line:**
```bash
--run_assembly_evaluation true \
--reference_genome /path/to/WT_mito.fa
```

**2. (Optional) Add transcript mapping:**
```bash
--transcripts_fasta /path/to/phaffi_transcripts.fa
```

**3. The module integrates after the FLYE assembly step** - See `ASSEMBLY_EVALUATION.md` for complete integration code

### 📋 Parameters Added to nextflow.config

```groovy
params {
    // Assembly evaluation parameters
    run_assembly_evaluation = false  // Set to true to run assembly evaluation
    reference_genome = null  // Path to WT reference genome (e.g., WT_mito.fa)
    transcripts_fasta = null  // Optional: Path to transcripts file for mapping
    min_assembly_depth = 20  // Minimum depth for assembly evaluation
    max_assembly_depth = 100  // Maximum depth for assembly evaluation
}
```

### ⚙️ Process Resources Configured

| Process | CPUs | Memory | Time |
|---------|------|--------|------|
| ALIGN_ASSEMBLY_TO_GENOME | 4 | 8 GB | 2h |
| REPAIR_ASSEMBLY | 2 | 4 GB | 1h |
| ALIGN_ANNOTATED_ASSEMBLY | 4 | 8 GB | 2h |
| FINALIZE_ASSEMBLY | 2 | 4 GB | 1h |
| ALIGN_FINAL_ASSEMBLY | 4 | 8 GB | 2h |
| MAP_READS_TO_ASSEMBLY | 8 | 16 GB | 4h |
| BLAST_TRANSGENE_TO_ASSEMBLY | 2 | 4 GB | 1h |
| CONVERT_BLAST_TO_BED | 1 | 2 GB | 30m |
| MAP_TRANSCRIPTS_TO_ASSEMBLY | 4 | 8 GB | 2h |

### 📝 Next Steps

1. **Review the files** created in your repository
2. **Read `ASSEMBLY_EVALUATION.md`** for detailed integration instructions
3. **Integrate into your main workflow** using the provided code examples
4. **Commit and push** these changes to your `with_eval_v5` branch

### 🔍 Quality Checks Performed

- ✅ Nextflow lint validation passed
- ✅ All processes follow DSL2 conventions
- ✅ Version tracking implemented for all tools
- ✅ Proper publishDir settings for all outputs
- ✅ Container specifications for all processes
- ✅ Resource allocations configured

### 📚 Documentation Files

- **`ASSEMBLY_EVALUATION.md`** - Complete user guide with:
  - Detailed process descriptions
  - Integration examples
  - Configuration options
  - Troubleshooting guide
  - Resource requirements
  - Output structure
  - Citation information

### 🎯 What This Module Does

Based on your original bash script, this module:

1. **Aligns** your assembly to a reference genome to understand contig structure
2. **Repairs** the assembly by orienting contigs correctly (forward/reverse complement)
3. **Validates** the repaired assembly with a second alignment
4. **Finalizes** the assembly by concatenating multiple contigs per chromosome
5. **Confirms** the final assembly quality with a third alignment
6. **Maps** your original ONT reads back to the final assembly for validation
7. **Detects** transgene insertions using BLAST
8. **Converts** BLAST results to BED format for genome browsers
9. **Maps** transcripts to the assembly (optional step)

All of this is now available as reusable, containerized Nextflow processes!

---

## Generated: 2025-12-23 07:41:24
## Branch: with_eval_v5
## Repository: https://github.com/KochInstitute-Bioinformatics/ont-flye-workflow
