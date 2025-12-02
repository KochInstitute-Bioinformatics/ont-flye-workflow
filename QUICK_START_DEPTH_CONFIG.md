# Quick Start: Configuring Assembly Depth Thresholds

## TL;DR

Add these two lines to your `nextflow.config` under the `params` section:

```groovy
params {
    min_assembly_depth = 10    // Minimum coverage for assembly
    max_assembly_depth = 300   // Maximum coverage for assembly
}
```

## Three Ways to Use Custom Thresholds

### 1. Edit your config file (Permanent)

Add to `../nextflow.config`:

```groovy
params {
    min_assembly_depth = 5
    max_assembly_depth = 500
}
```

Run normally:

```bash
nextflow run main.nf -c ../nextflow.config -profile singularity,slurm
```

### 2. Command line override (One-time)

Keep your config as-is, override when running:

```bash
nextflow run main.nf \
    -c ../nextflow.config \
    -profile singularity,slurm \
    --min_assembly_depth 5 \
    --max_assembly_depth 500
```

### 3. Update your SBATCH script (Convenient)

Modify your submission script:

```bash
# Run the workflow with custom thresholds
nextflow run main.nf \
    -c ../nextflow.config \
    -profile singularity,slurm \
    --min_assembly_depth 5 \
    --max_assembly_depth 500
```

## Common Use Cases

**Include low-coverage samples:**

```bash
--min_assembly_depth 3 --max_assembly_depth 300
```

**Handle high-coverage samples:**

```bash
--min_assembly_depth 10 --max_assembly_depth 1000
```

**No upper limit:**

```bash
--min_assembly_depth 10 --max_assembly_depth 999999
```

**Only high-quality assemblies:**

```bash
--min_assembly_depth 20 --max_assembly_depth 200
```

## Output Files

Check these files to see which samples were selected/filtered:

- `assembly_candidates.csv` - Selected for assembly
- `assembly_filtered.csv` - Filtered out with reasons

## Options for Adjusting Acceptable Depths

1. Set defaults in config file
2. Override with `--min_assembly_depth` and `--max_assembly_depth` on command line
3. Or edit the config file for permanent changes
