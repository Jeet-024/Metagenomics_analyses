#!/bin/bash

set -euo pipefail

INPUT="/mnt/faruq2/lab_data/raw_data/rRNA-Seq/241027_M05986_0170_LG3CY"
OUTPUT="/mnt/faruq2/lab_users/manojeet_docker/rRNA_Seq_analysis/241027_M05986_0170_LG3CY/raw_fastq_output"

mkdir -p "$OUTPUT"

bcl-convert \
    --bcl-input-directory "$INPUT" \
    --output-directory "$OUTPUT" \
    --sample-sheet /mnt/faruq2/lab_data/raw_data/rRNA-Seq/241027_M05986_0170_LG3CY/all_50_sample_sheet_rRNA_bcl2fastq_run.csv \
    --bcl-num-compression-threads 20 \
    --force
