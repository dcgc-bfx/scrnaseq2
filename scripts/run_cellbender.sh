#!/usr/bin/env bash

# Bash script for running CellBender remove-background
# CellBender is a software package for eliminating technical artifacts from 
# high-throughput single-cell RNA sequencing (scRNA-seq) data.
#
# CellBender documentation: https://cellbender.readthedocs.io/
# CellBender GitHub: https://github.com/broadinstitute/CellBender

set -euo pipefail

# Default values
INPUT_FILE=""
OUTPUT_DIR=""
EXPECTED_CELLS=""
TOTAL_DROPLETS=""
FPR="0.01"
EPOCHS="150"
CUDA=""
LEARNING_RATE="0.0001"
Z_DIM="100"
Z_LAYERS="500"
LOW_COUNT_THRESHOLD="15"
EXTRA_ARGS=""

# Function to display usage
usage() {
    cat << EOF
Usage: $(basename "$0") [OPTIONS]

Run CellBender remove-background to eliminate ambient RNA contamination from scRNA-seq data.

Required Arguments:
    -i, --input FILE            Input raw h5 file from Cell Ranger (e.g., raw_feature_bc_matrix.h5)
    -o, --output-dir DIR        Output directory for results
    -e, --expected-cells NUM    Expected number of real cells
    -t, --total-droplets NUM    Total droplets to include (usually 2-3x expected cells)

Optional Arguments:
    -f, --fpr FLOAT            Target false positive rate (default: 0.01)
    --epochs NUM               Number of training epochs (default: 150)
    --cuda                     Use CUDA for GPU acceleration
    --learning-rate FLOAT      Learning rate for training (default: 0.0001)
    --z-dim NUM                Dimension of latent variable z (default: 100)
    --z-layers NUM             Dimension of hidden layers (default: 500)
    --low-count-threshold NUM  Threshold for excluding low count features (default: 15)
    --extra-args "ARGS"        Additional arguments to pass to cellbender
    -h, --help                 Display this help message

Example:
    $(basename "$0") \\
        --input raw_feature_bc_matrix.h5 \\
        --output-dir cellbender_output \\
        --expected-cells 5000 \\
        --total-droplets 15000 \\
        --cuda

Notes:
    - CellBender must be installed and available in your environment
    - For GPU acceleration, use --cuda flag (requires CUDA-compatible GPU)
    - Total droplets should be 2-3x the expected number of cells
    - Output will be saved as: {output-dir}/output_filtered.h5

EOF
    exit 0
}

# Parse command line arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        -i|--input)
            INPUT_FILE="$2"
            shift 2
            ;;
        -o|--output-dir)
            OUTPUT_DIR="$2"
            shift 2
            ;;
        -e|--expected-cells)
            EXPECTED_CELLS="$2"
            shift 2
            ;;
        -t|--total-droplets)
            TOTAL_DROPLETS="$2"
            shift 2
            ;;
        -f|--fpr)
            FPR="$2"
            shift 2
            ;;
        --epochs)
            EPOCHS="$2"
            shift 2
            ;;
        --cuda)
            CUDA="--cuda"
            shift
            ;;
        --learning-rate)
            LEARNING_RATE="$2"
            shift 2
            ;;
        --z-dim)
            Z_DIM="$2"
            shift 2
            ;;
        --z-layers)
            Z_LAYERS="$2"
            shift 2
            ;;
        --low-count-threshold)
            LOW_COUNT_THRESHOLD="$2"
            shift 2
            ;;
        --extra-args)
            EXTRA_ARGS="$2"
            shift 2
            ;;
        -h|--help)
            usage
            ;;
        *)
            echo "Error: Unknown option: $1" >&2
            echo "Use -h or --help for usage information" >&2
            exit 1
            ;;
    esac
done

# Validate required arguments
if [[ -z "$INPUT_FILE" ]]; then
    echo "Error: Input file is required (-i or --input)" >&2
    exit 1
fi

if [[ -z "$OUTPUT_DIR" ]]; then
    echo "Error: Output directory is required (-o or --output-dir)" >&2
    exit 1
fi

if [[ -z "$EXPECTED_CELLS" ]]; then
    echo "Error: Expected number of cells is required (-e or --expected-cells)" >&2
    exit 1
fi

if [[ -z "$TOTAL_DROPLETS" ]]; then
    echo "Error: Total droplets is required (-t or --total-droplets)" >&2
    exit 1
fi

# Validate input file exists
if [[ ! -f "$INPUT_FILE" ]]; then
    echo "Error: Input file does not exist: $INPUT_FILE" >&2
    exit 1
fi

# Create output directory if it doesn't exist
mkdir -p "$OUTPUT_DIR"

# Check if cellbender is installed
if ! command -v cellbender &> /dev/null; then
    echo "Error: cellbender is not installed or not in PATH" >&2
    echo "Please install CellBender: pip install cellbender" >&2
    exit 1
fi

# Display configuration
echo "========================================"
echo "CellBender Remove-Background"
echo "========================================"
echo "Input file:         $INPUT_FILE"
echo "Output directory:   $OUTPUT_DIR"
echo "Expected cells:     $EXPECTED_CELLS"
echo "Total droplets:     $TOTAL_DROPLETS"
echo "False positive rate: $FPR"
echo "Epochs:             $EPOCHS"
echo "Learning rate:      $LEARNING_RATE"
echo "Z dimension:        $Z_DIM"
echo "Z layers:           $Z_LAYERS"
echo "Low count threshold: $LOW_COUNT_THRESHOLD"
if [[ -n "$CUDA" ]]; then
    echo "GPU acceleration:   Enabled"
else
    echo "GPU acceleration:   Disabled"
fi
echo "========================================"
echo ""

# Build the cellbender command
CMD="cellbender remove-background \
    --input \"$INPUT_FILE\" \
    --output \"$OUTPUT_DIR/output.h5\" \
    --expected-cells $EXPECTED_CELLS \
    --total-droplets-included $TOTAL_DROPLETS \
    --fpr $FPR \
    --epochs $EPOCHS \
    --learning-rate $LEARNING_RATE \
    --z-dim $Z_DIM \
    --z-layers $Z_LAYERS \
    --low-count-threshold $LOW_COUNT_THRESHOLD"

# Add CUDA flag if specified
if [[ -n "$CUDA" ]]; then
    CMD="$CMD $CUDA"
fi

# Add extra arguments if specified
if [[ -n "$EXTRA_ARGS" ]]; then
    CMD="$CMD $EXTRA_ARGS"
fi

# Run CellBender
echo "Running CellBender..."
echo "Command: $CMD"
echo ""

eval "$CMD"

# Check if output was created
OUTPUT_FILE="$OUTPUT_DIR/output_filtered.h5"
if [[ -f "$OUTPUT_FILE" ]]; then
    echo ""
    echo "========================================"
    echo "CellBender completed successfully!"
    echo "Output file: $OUTPUT_FILE"
    echo "========================================"
else
    echo ""
    echo "========================================"
    echo "Warning: Expected output file not found: $OUTPUT_FILE"
    echo "Check the output directory for results: $OUTPUT_DIR"
    echo "========================================"
fi

exit 0
