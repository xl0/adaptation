#!/bin/bash

OUT_DIR=out

mkdir -p ${OUT_DIR}

echo "### Extracting spacers ..."
for experiment in $(ls reads); do
	echo "Eextracting spacers from ${experiment}"
	mkdir -p ${OUT_DIR}/spacers/${experiment}
	time python2.7 -u ./extract-spacers-mp.py reads/${experiment}/*.fastq* ${OUT_DIR}/spacers/${experiment}/ || exit 1
done
