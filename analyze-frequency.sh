#!/bin/bash

OUT_DIR=../nas/Alex/out4

rm -fr ${OUT_DIR}/analysis/

for experiment in $(ls ${OUT_DIR}/alignments/); do
	for input in $(ls ${OUT_DIR}/alignments/${experiment}/1phiNM4gamma4*.gb); do
		name=$(basename ${input} .gb)
		mkdir -p ${OUT_DIR}/analysis/${experiment}
		output=${OUT_DIR}/analysis/${experiment}/${name}

		echo ""
		./analyze-pam-frequency.py -o ${output}.json ${input} || exit 1

		./graph_pam_frequency.py -o ${output}.svg ${output}.json
		./graph_pam_motif.py -o ${output}_motif ${output}.json

	done
	./pam_frequency_table.py ${OUT_DIR}/analysis/${experiment}/*.json \
		-o ${OUT_DIR}/analysis/${experiment}-pam-hit-table.txt \
		-j ${OUT_DIR}/analysis/${experiment}-pam-hit-table.json \
		-c ${OUT_DIR}/analysis/${experiment}-pam-hit-table.csv
	echo ""
	./graph_pam_frequency.py -o ${OUT_DIR}/analysis/${experiment}-pam-hits.svg ${OUT_DIR}/analysis/${experiment}/*.json
	./graph_pam_frequency_pairs.py -o ${OUT_DIR}/analysis/${experiment}-pam-hits-correlations.svg ${OUT_DIR}/analysis/${experiment}/*.json

done

./pam_frequency_table.py ${OUT_DIR}/analysis/*/*.json -o ${OUT_DIR}/analysis/pam-hit-table.txt \
		-j ${OUT_DIR}/analysis/pam-hit-table.json -c ${OUT_DIR}/analysis/pam-hit-table.csv
echo ""
./graph_pam_frequency.py -o ${OUT_DIR}/analysis/pam-hits.svg ${OUT_DIR}/analysis/*/*.json

#./graph_pam_frequency_pairs.py -o ${OUT_DIR}/analysis/pam-hits-correlations.svg ${OUT_DIR}/analysis/*/*.json

./graph_pam_frequency_pairs.py -o ${OUT_DIR}/analysis/pam-hits-wt.svg ${OUT_DIR}/analysis/*/1phiNM4gamma4-PAMS-CAG.json ${OUT_DIR}/analysis/*/1phiNM4gamma4-PAMS-GAC.json 

./graph_pam_frequency_pairs.py -o ${OUT_DIR}/analysis/pam-hits-pRH135.svg ${OUT_DIR}/analysis/*/1phiNM4gamma4-PAMS-TCA.json ${OUT_DIR}/analysis/*/1phiNM4gamma4-PAMS-ACT.json

./graph_pam_frequency_pairs.py -o ${OUT_DIR}/analysis/pam-hits-pRH196.svg ${OUT_DIR}/analysis/*/1phiNM4gamma4-PAMS-GTC.json ${OUT_DIR}/analysis/*/1phiNM4gamma4-PAMS-CTG.json

./graph_pam_frequency_pairs.py -o ${OUT_DIR}/analysis/pam-hits-pRH163.svg ${OUT_DIR}/analysis/*/1phiNM4gamma4-PAMS-AGT.json ${OUT_DIR}/analysis/*/1phiNM4gamma4-PAMS-TGA.json

#for experiment in $(ls ${OUT_DIR}/alignments/); do
#	for input in $(ls ${OUT_DIR}/alignments/${experiment}/2pWJ40*.gb); do
#		name=$(basename ${input} .gb)
#		mkdir -p ${OUT_DIR}/analysis/${experiment}
#		output=${OUT_DIR}/analysis/${experiment}/${name}
#
#		echo ""
#		./analyze-pam-frequency.py -o ${output}.json ${input} || exit 1
#	done
#	./pam_frequency_table.py ${OUT_DIR}/analysis/${experiment}/*.json -o \
#		${OUT_DIR}/analysis/${experiment}-pam-hit-table.txt -j ${OUT_DIR}/analysis/${experiment}-pam-hit-table.json
#done


