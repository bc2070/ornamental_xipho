#!/usr/bin/env bash

#SBATCH --job-name=LK_pipeline
#SBATCH --partition=long
#SBATCH --cpus-per-task=8
#SBATCH --mem=40G

# ==========================================================
# Locate the pipeline directory
# ==========================================================

ROOT="${SLURM_SUBMIT_DIR:-$(pwd)}"

# ==========================================================
# Define directories
# ==========================================================

CONFIG="$ROOT/config.ini"
SCRIPT="$ROOT/scripts"
LOG="$ROOT/logs"
OUTPUT="$ROOT/output"

# ==========================================================
# Create required directories
# ==========================================================

mkdir -p "$LOG"
mkdir -p "$OUTPUT"

# ==========================================================
# Configure SLURM output files
# ==========================================================

# The log directory must exist before sbatch starts the job.
# Therefore, create it before submitting the job:
#     mkdir -p logs
#
# SLURM output/error files are stored in the logs directory.

#SBATCH --output=logs/LK_pipeline_%j.out
#SBATCH --error=logs/LK_pipeline_%j.err

# ==========================================================
# Initialize shell
# ==========================================================

set -euo pipefail

cd "$ROOT"

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
echo "Logs          : $LOG"
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

