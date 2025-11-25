# Running CellBender for Ambient RNA Removal

## Overview

CellBender is a Python tool developed by the Broad Institute for removing ambient RNA contamination from single-cell RNA-seq data. This script provides a convenient wrapper for running CellBender's `remove-background` command.

## What is Ambient RNA?

Ambient RNA refers to cell-free RNA molecules present in the solution surrounding cells during the scRNA-seq experiment. These molecules can be captured along with the RNA from actual cells, leading to contamination of the gene expression profiles. CellBender uses a machine learning approach to identify and remove this technical artifact.

## Installation

Before using this script, you need to install CellBender:

```bash
# Create a conda environment (recommended)
conda create -n cellbender python=3.9
conda activate cellbender

# Install CellBender
pip install cellbender
```

For GPU acceleration (optional but recommended for large datasets):
```bash
# Install PyTorch with CUDA support first
pip install torch torchvision torchaudio --index-url https://download.pytorch.org/whl/cu118
pip install cellbender
```

For more details, see the [CellBender documentation](https://cellbender.readthedocs.io/en/latest/installation/index.html).

## Basic Usage

```bash
./scripts/run_cellbender.sh \
    --input raw_feature_bc_matrix.h5 \
    --output-dir cellbender_output \
    --expected-cells 5000 \
    --total-droplets 15000
```

## Required Arguments

- `--input` or `-i`: Path to the raw h5 file from Cell Ranger (e.g., `raw_feature_bc_matrix.h5`)
- `--output-dir` or `-o`: Directory where output files will be saved
- `--expected-cells` or `-e`: Expected number of real cells in your sample
- `--total-droplets` or `-t`: Total number of droplets to include in the analysis

## Optional Arguments

- `--fpr` or `-f`: False positive rate (default: 0.01). Lower values are more conservative.
- `--epochs`: Number of training epochs (default: 150). More epochs may improve results but take longer.
- `--cuda`: Enable GPU acceleration (requires CUDA-compatible GPU)
- `--learning-rate`: Learning rate for training (default: 0.0001)
- `--z-dim`: Dimension of latent variable z (default: 100)
- `--z-layers`: Dimension of hidden layers (default: 500)
- `--low-count-threshold`: Threshold for excluding low count features (default: 15)
- `--extra-args`: Additional arguments to pass directly to cellbender

## How to Determine Parameters

### Expected Cells

Check your Cell Ranger summary metrics or estimate based on:
- Your experimental design (how many cells you loaded)
- The knee plot in Cell Ranger output
- The "Estimated Number of Cells" in the Cell Ranger web summary

### Total Droplets

Typically set to 2-3x the expected number of cells. This should include:
- Real cells
- Empty droplets
- Droplets with ambient RNA

A good starting point: `total_droplets = expected_cells × 2.5`

### False Positive Rate (FPR)

- Lower FPR (e.g., 0.001): More conservative, removes less ambient RNA
- Higher FPR (e.g., 0.05): More aggressive, removes more ambient RNA
- Default (0.01): Good balance for most applications

## Examples

### Basic run without GPU
```bash
./scripts/run_cellbender.sh \
    --input /path/to/raw_feature_bc_matrix.h5 \
    --output-dir /path/to/output \
    --expected-cells 5000 \
    --total-droplets 12500
```

### Run with GPU acceleration
```bash
./scripts/run_cellbender.sh \
    --input /path/to/raw_feature_bc_matrix.h5 \
    --output-dir /path/to/output \
    --expected-cells 5000 \
    --total-droplets 12500 \
    --cuda
```

### Custom parameters
```bash
./scripts/run_cellbender.sh \
    --input /path/to/raw_feature_bc_matrix.h5 \
    --output-dir /path/to/output \
    --expected-cells 10000 \
    --total-droplets 25000 \
    --fpr 0.005 \
    --epochs 200 \
    --cuda
```

### Large dataset with custom parameters
```bash
./scripts/run_cellbender.sh \
    --input /path/to/raw_feature_bc_matrix.h5 \
    --output-dir /path/to/output \
    --expected-cells 20000 \
    --total-droplets 50000 \
    --epochs 200 \
    --learning-rate 0.00005 \
    --cuda
```

## Output Files

CellBender creates several output files in the specified output directory:

- `output_filtered.h5`: Main output file with ambient RNA removed (use this for downstream analysis)
- `output.h5`: Contains both filtered and raw counts
- `output_cell_barcodes.csv`: List of cell barcodes that passed filtering
- `output.pdf`: Diagnostic plots
- `output.log`: Log file with training details

## Using CellBender Output in Seurat

After running CellBender, you can read the filtered data into Seurat:

```R
library(Seurat)

# Read CellBender output
counts <- Read10X_h5("cellbender_output/output_filtered.h5")

# Create Seurat object
seurat_obj <- CreateSeuratObject(
    counts = counts,
    project = "my_project"
)
```

## Troubleshooting

### "cellbender is not installed or not in PATH"
Make sure you've activated the conda environment where CellBender is installed:
```bash
conda activate cellbender
```

### Out of memory errors
- Reduce `--total-droplets` to include fewer droplets
- Use a machine with more RAM
- If using GPU, try reducing `--z-dim` or `--z-layers`

### Training takes too long
- Reduce `--epochs` (try 100 instead of 150)
- Use GPU acceleration with `--cuda`
- Reduce `--total-droplets`

### Results look over-corrected
- Decrease `--fpr` to be more conservative (e.g., 0.001)
- Increase `--epochs` for better convergence

### Results look under-corrected
- Increase `--fpr` to be more aggressive (e.g., 0.05)
- Check that `--total-droplets` includes enough empty droplets

## Best Practices

1. **Always use the raw h5 file**: CellBender needs both cell-containing and empty droplets to learn the ambient RNA profile.

2. **Start with default parameters**: The defaults work well for most datasets. Only adjust if you have specific reasons.

3. **Check the diagnostic plots**: Review the `output.pdf` file to ensure the model converged properly.

4. **Compare with original data**: It's good practice to compare your results before and after CellBender to understand what was removed.

5. **Use GPU if available**: GPU acceleration significantly speeds up processing, especially for large datasets.

## References

- CellBender paper: [Fleming et al., 2023, Nature Methods](https://www.nature.com/articles/s41592-023-01943-7)
- CellBender documentation: https://cellbender.readthedocs.io/
- CellBender GitHub: https://github.com/broadinstitute/CellBender

## Support

For issues with the script, please contact the DRESDEN-concept Genome Center bioinformatics team.

For CellBender-specific issues, refer to:
- [CellBender GitHub Issues](https://github.com/broadinstitute/CellBender/issues)
- [CellBender Documentation](https://cellbender.readthedocs.io/)
