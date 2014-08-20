#!/bin/bash

OUT_DIR=/home/xl0/work/CRI/internships/Pasteur/nas/Alex/out2

rm -fr ${OUT_DIR}/analysis/
for experiment in $(ls ${OUT_DIR}/alignments/); do
	for input in $(ls ${OUT_DIR}/alignments/${experiment}/*phiNM4*.gb); do
		name=$(basename ${input} .gb)
		mkdir -p ${OUT_DIR}/analysis/${experiment}
		output=${OUT_DIR}/analysis/${experiment}/${name}

		./analyze-pam-frequency.py -o ${output} ${input} || exit 1
	done
	./pam_frequency_table.py ${OUT_DIR}/analysis/${experiment}/*.json -o \
		${OUT_DIR}/analysis/${experiment}-pam-hit-table.txt -j ${OUT_DIR}/analysis/${experiment}-pam-hit-table.json
done

./pam_frequency_table.py ${OUT_DIR}/analysis/*/*.json -o ${OUT_DIR}/analysis/pam-hit-table.txt -j ${OUT_DIR}/analysis/pam-hit-table.json

