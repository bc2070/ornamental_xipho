#!/usr/bin/env bash

#SBATCH --job-name=LK_pipeline
#SBATCH --partition=long
#SBATCH --cpus-per-task=8
#SBATCH --mem=40G
#SBATCH --output=LK_pipeline_%j.out
#SBATCH --error=LK_pipeline_%j.err

set -euo pipefail

# ==========================================================
# Locate the pipeline directory
# ==========================================================

# SLURM_SUBMIT_DIR is the directory from which sbatch was called.
# This avoids using the temporary directory created by SLURM.
ROOT="${SLURM_SUBMIT_DIR:-$(pwd)}"

cd "$ROOT"

# ==========================================================
# Define directories
# ==========================================================

CONFIG="$ROOT/config.ini"
SCRIPT="$ROOT/scripts"
LOG="$ROOT/logs"
OUTPUT="$ROOT/output"

# ==========================================================
# Create output directories
# ==========================================================

mkdir -p "$LOG"
mkdir -p "$OUTPUT"

# ==========================================================
# Check required files and directories
# ==========================================================

if [[ ! -f "$CONFIG" ]]; then
    echo "ERROR: config.ini not found:"
    echo "$CONFIG"
    exit 1
fi

if [[ ! -d "$SCRIPT" ]]; then
    echo "ERROR: scripts directory not found:"
    echo "$SCRIPT"
    exit 1
fi

if [[ ! -f "$SCRIPT/run_pipeline.py" ]]; then
    echo "ERROR: run_pipeline.py not found:"
    echo "$SCRIPT/run_pipeline.py"
    exit 1
fi

# ==========================================================
# Load required modules
# ==========================================================

module load bedtools
module load R

# ==========================================================
# Print job information
# ==========================================================

echo "=================================================="
echo "LK pipeline started"
echo "=================================================="
echo "Job ID        : ${SLURM_JOB_ID:-N/A}"
echo "Job name      : ${SLURM_JOB_NAME:-N/A}"
echo "Submit dir    : ${SLURM_SUBMIT_DIR:-N/A}"
echo "Pipeline root : $ROOT"
echo "Config        : $CONFIG"
echo "Scripts       : $SCRIPT"
echo "Output        : $OUTPUT"
echo "=================================================="

# ==========================================================
# Run pipeline
# ==========================================================

python3 "$SCRIPT/run_pipeline.py" "$CONFIG"

# ==========================================================
# Finish
# ==========================================================

echo "=================================================="
echo "LK pipeline completed"
echo "=================================================="

