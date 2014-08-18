#!/bin/bash

rm -fr out tmp
mkdir out tmp


echo "### Annotating PAMs ..."
mkdir -p tmp/pam-annotated-templates/
for template in $(ls templates/*.gb); do
	./annotate-pams.py $template tmp/pam-annotated-templates/$(basename $template .gb)-PAMS.gb
done

echo "### Extracting spacers ..."
for experiment in $(ls reads); do
	echo "Eextracting spacers from $experiment"
	mkdir -p out/spacers/$experiment
	./extract-spacers-mp.py reads/$experiment/*.fastq* out/spacers/$experiment/
done

echo "### Aligning spacers to templates ..."
for experiment in $(ls out/spacers/); do
	for spacer_file in $(ls out/spacers/$experiment); do
		tag=$(basename $spacer_file .fastq)
		echo "Aligning spacers from $experiment:$tag"
		mkdir -p out/alignments/$experiment/$tag

		for template in $(ls tmp/pam-annotated-templates/*.gb); do
			cp $template out/alignments/$experiment/$tag/$(basename $template .gb)-$tag.gb
			echo "$template -> out/alignments/$experiment/$tag/$(basename $template .gb)-$tag.gb"
#			echo $template
		done
		
		./align-spacers-mp.py out/spacers/$experiment/$spacer_file $(ls out/alignments/$experiment/$tag/*.gb)

	done
done

rm -fr tmp
