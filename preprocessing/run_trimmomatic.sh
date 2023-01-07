#!/bin/bash

DIR=$(cd ../../ && pwd)
TRIMMOMATIC=$HOME/app/Trimmomatic-0.39/trimmomatic-0.39.jar
ADAPTERS=$DIR/candida-infection-pipeline/preprocessing/adapters.fasta
SAMPLES=$DIR/fastq_samples
OUTDIR=$DIR/trimmed_samples

ID=(
	SRR7821862.1 SRR7821863.1 SRR7821864.1 SRR7821865.1 SRR7821866.1 SRR7821867.1 SRR7821868.1 SRR7821869.1 \
	SRR7821870.1 SRR7821871.1 SRR7821872.1 SRR7821873.1 SRR7821874.1 SRR7821875.1 SRR7821876.1 SRR7821877.1 \
	SRR7821878.1 SRR7821879.1 SRR7821880.1
)

START=$(date)

for i in ${ID[@]}; do
	java -jar $TRIMMOMATIC PE -threads 6 -phred33 \
	$SAMPLES/${i}_1.fastq $SAMPLES/${i}_2.fastq \
	-baseout $OUTDIR/${i}.fq.gz \
	ILLUMINACLIP:$ADAPTERS:2:30:10 \
	CROP:74 HEADCROP:1 \
	SLIDINGWINDOW:4:15 MINLEN:36 
done

echo "Start: $START"
echo "End: $(date)"
