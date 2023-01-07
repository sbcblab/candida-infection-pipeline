#!/bin/bash

DIR=$(cd ../../ && pwd)

RSEM=$HOME/app/RSEM
STAR=$HOME/app/STAR/bin/Linux_x86_64_static

READS_DIR=$DIR/trimmed_samples
OUTDIR=$DIR/rsem_output
SAMPLE=$1

START=$(date)

$RSEM/rsem-calculate-expression -p 12 \
	--paired-end \
        --no-bam-output \
        --star \
        --star-path $STAR \
        --star-gzipped-read-file \
        $READS_DIR/${SAMPLE}_1P.fq.gz $READS_DIR/${SAMPLE}_2P.fq.gz \
        $OUTDIR/reference/mm_ensembl/mus_musculus \
        $OUTDIR/mm_${SAMPLE}_quals

echo "-------------------------------------------------------------------------"
echo "Start: $START"
echo "End: $(date)"
