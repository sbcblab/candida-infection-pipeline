#!/bin/bash

DIR=$(cd ../../ && pwd)
OUTPUT_DIR=$DIR/samples
START=$(date)

wget --input-file=$DIR/candida-infection-pipeline/preprocessing/datasets_url.txt --directory-prefix=$OUTPUT_DIR

echo "Start: $START"
echo "End: $(date)"
