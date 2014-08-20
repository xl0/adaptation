#!/bin/bash

OUT_DIR=out


rm -fr tmp
mkdir -p ${OUT_DIR} tmp


echo "### Annotating PAMs ..."
mkdir -p tmp/pam-annotated-templates/
for template in $(ls templates/*.gb); do
	./annotate-pams.py $template tmp/pam-annotated-templates/$(basename $template .gb)-PAMS.gb || exit 1
done

echo "### Extracting spacers ..."
for experiment in $(ls reads); do
	echo "Eextracting spacers from ${experiment}"
	mkdir -p ${OUT_DIR}/spacers/${experiment}
	time python2.7 -u ./extract-spacers-mp.py reads/${experiment}/*.fastq* ${OUT_DIR}/spacers/${experiment}/ || exit 1
done


echo "### Aligning spacers to templates ..."
rm -fr ${OUT_DIR}/alignments
for experiment in $(ls ${OUT_DIR}/spacers/); do
	for spacer_file in $(ls ${OUT_DIR}/spacers/${experiment}/*.fastq); do
		tag=$(basename ${spacer_file} .fastq)
		echo "Aligning spacers from ${experiment}:${tag}"
		mkdir -p ${OUT_DIR}/alignments/${experiment}/${tag}

		for template in $(ls tmp/pam-annotated-templates/*.gb); do
			cp $template ${OUT_DIR}/alignments/${experiment}/${tag}/$(basename $template .gb)-${tag}.gb
			echo "$template -> ${OUT_DIR}/alignments/${experiment}/${tag}/$(basename $template .gb)-${tag}.gb"
		done
       	
		time python2.7 -u ./align-spacers-mp.py ${OUT_DIR}/spacers/${experiment}/${spacer_file} $(ls ${OUT_DIR}/alignments/${experiment}/${tag}/*.gb) || exit 1

	done
done

for experiment in $(ls ${OUT_DIR}/alignments/); do
	for tag in $(ls ${OUT_DIR}/alignments/${experiment}); do
		for result in $(ls ${OUT_DIR}/alignments/${experiment}/${tag}/*phiNM4*.gb); do
			echo ${result}
			./analyze-pam-frequency.py ${result} || exit 1
		done
	done
done


rm -fr tmp
