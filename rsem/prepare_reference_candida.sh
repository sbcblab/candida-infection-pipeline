#!/bin/bash

DIR=$(cd ../../ && pwd)

RSEM=$HOME/Programas/RSEM-1.3.3
STAR=$HOME/Programas/STAR-2.7.10a/bin/Linux_x86_64_static

REF_DIR=$DIR/reference
OUT_DIR=$DIR/rsem_output

START=$(date)

$RSEM/rsem-prepare-reference --gtf $REF_DIR/ca_GCF_000182965.3_ASM18296v3_genomic.gtf \
	--star \
        --star-path $STAR \
        --p 8 \
        $REF_DIR/ca_GCF_000182965.3_ASM18296v3_genomic.fna \
        $OUT_DIR/reference/ca_ensembl/candida_albicans

echo "Start: $START"
echo "End: $(date)"

