#!/bin/bash

INPUT="/home/manojeet_docker/rRNA_Seq_analysis/241027_M05986_0170_LG3CY/raw_fastq_output"
OUTPUT="/home/manojeet_docker/rRNA_Seq_analysis/241027_M05986_0170_LG3CY/fastqc_report_output"

mkdir -p "$OUTPUT"

for file in "$INPUT"/*.fastq.gz;
do
   if [[ "$file" == *"Undetermined"* ]]; then
   	continue
   fi

   fastqc -t 10 "$file" -o "$OUTPUT"
done
