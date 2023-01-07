#!/bin/bash

cd ..
DIR=$(pwd)

MODULES=(turquoise blue)

for MOD in ${MODULES[@]}; do
	FILE=$DIR/results/CytoscapeInput-edges-$MOD.txt
	LNCRNA=$DIR/results/degs_mmgs_ml_lncrna_$MOD.txt
	OUTPUT_FILE=$DIR/results/CytoscapeInput-edges-$MOD-lncrna.txt

	(head -1 $FILE; grep -F -f $LNCRNA $FILE) > $OUTPUT_FILE
done

