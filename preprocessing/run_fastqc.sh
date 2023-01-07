#!/bin/bash

DIR=$(cd ../../ && pwd)
START=$(date)

$HOME/app/FastQC/fastqc -t 3 $DIR/fastq_samples/* -o $DIR/fastqc_original
#$HOME/app/FastQC/fastqc -t 3 $DIR/trimmed_samples/*P.fq.gz -o $DIR/fastqc_trimmed

echo "Start: $START"
echo "End: $(date)"

