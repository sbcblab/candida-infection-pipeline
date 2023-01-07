#!/bin/bash

DIR=$(cd ../../ && pwd)

RSEM=$HOME/Programas/RSEM-1.3.3
STAR=$HOME/Programas/STAR-2.7.10a/bin/Linux_x86_64_static

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
        $OUTDIR/reference/ca_ensembl/candida_albicans \
        $OUTDIR/expression/ca/ca_${SAMPLE}_quals

echo "-------------------------------------------------------------------------"
echo "Start: $START"
echo "End: $(date)"
