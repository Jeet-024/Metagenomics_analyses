#!/bin/bash

#--fn_as_s_name Use the log filename as the sample name
#--module specifies to use a particular module only, here fastqc

INPUT="/home/manojeet_docker/rRNA_Seq_analysis/241027_M05986_0170_LG3CY/fastqc_report_output"
OUTPUT="/home/manojeet_docker/rRNA_Seq_analysis/241027_M05986_0170_LG3CY/MultiQC_report"

multiqc -o "$OUTPUT" \
--module fastqc \
--fn_as_s_name \
"$INPUT"
