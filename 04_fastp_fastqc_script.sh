#!/bin/bash

INPUT="/home/manojeet_docker/rRNA_Seq_analysis/241027_M05986_0170_LG3CY/raw_fastqs"
OUTPUT="/home/manojeet_docker/rRNA_Seq_analysis/241027_M05986_0170_LG3CY/trimmed_fastq_output"
OUT_QC_TRIM="/home/manojeet_docker/rRNA_Seq_analysis/241027_M05986_0170_LG3CY/trimmed_QC_report"

mkdir -p "$OUTPUT"
mkdir -p "$OUT_QC_TRIM"

for f1 in "$INPUT"/*_R1_*.fastq.gz;
do
   [[ "$f1" == *Undetermined* ]] && continue
   
       f2="${f1/_R1_/_R2_}"
   
       if [[ ! -f "$f2" ]]; then
           echo "Missing pair for $f1"
           continue
       fi
   
       sample=$(basename "${f1%%_L001_R*}")
   
    echo "Processing $sample"
   fastp \
    -i "$f1" \
    -I "$f2" \
    -o "$OUTPUT/${sample}_trimmed_L001_R1.fastq.gz" \
    -O "$OUTPUT/${sample}_trimmed_L001_R2.fastq.gz" \
    --detect_adapter_for_pe \
    --thread 20 \
    -h "$OUTPUT/${sample}_report.html" \
    -j "$OUTPUT/${sample}_report.json"
   
   echo "Running Trimmed QC..."
   
   fastqc \
    -t 10 \
    -o "$OUT_QC_TRIM" \
    "$OUTPUT/${sample}_trimmed_L001_R1.fastq.gz" \
    "$OUTPUT/${sample}_trimmed_L001_R2.fastq.gz"
   
done
   
echo "Workflow complete!"
