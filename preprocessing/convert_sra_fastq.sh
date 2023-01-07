#!/bin/bash

DIR=$(cd ../../ && pwd)
FASTQDUMP=$HOME/app/sratoolkit.2.10.5-ubuntu64/bin/fastq-dump
OUTDIR=$DIR/fastq_samples
START=$(date)

fastq-dump --split-files $DIR/samples/* --outdir $OUTDIR

echo "Start: $START"
echo "End: $(date)"
