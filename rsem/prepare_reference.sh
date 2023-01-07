#!/bin/bash

DIR=$(cd ../../ && pwd)

RSEM=$HOME/app/RSEM
STAR=$HOME/app/STAR/bin/Linux_x86_64_static

REF_DIR=$DIR/reference
OUT_DIR=$DIR/rsem_output

START=$(date)

$RSEM/rsem-prepare-reference --gtf $REF_DIR/Mus_musculus.GRCm39.103.gtf \
	--star \
        --star-path $STAR \
        --p 8 \
        $REF_DIR/Mus_musculus.GRCm39.dna.primary_assembly.fa \
        $OUT_DIR/reference/mm_ensembl/mus_musculus

echo "Start: $START"
echo "End: $(date)"

