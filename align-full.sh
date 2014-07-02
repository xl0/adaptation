#!/bin/bash

set -x

rm -fr out
mkdir out

rm -fr tmp
mkdir tmp

# Pre-generate PAMs.

./annotate-pams.py phiNM4.gb out/phiNM4-PAMS.gb
./annotate-pams.py pWJ40.gb out/pWJ40-PAMS.gb
./annotate-pams.py NCTC8325.gb out/NCTC8325-PAMS.gb

for i in AGT CAG GTC TCA;
do
	cp out/phiNM4-PAMS.gb out/phiNM4-NoIndex_L001_${i}-PAMS.gb
	cp out/pWJ40-PAMS.gb out/pWJ40-NoIndex_L001_${i}-PAMS.gb
	cp out/NCTC8325-PAMS.gb out/NCTC8325-NoIndex_L001_${i}-PAMS.gb
	time python2.7 -u ./align-spacers-mp.py spacers_NoIndex_L001_${i}.fastq out/phiNM4-NoIndex_L001_${i}-PAMS.gb out/pWJ40-NoIndex_L001_${i}-PAMS.gb out/NCTC8325-NoIndex_L001_${i}-PAMS.gb
done

exit 0
for i in AGT CAG GTC TCA;
do
	cp out/phiNM4-PAMS.gb out/phiNM4-S1_L001_${i}-PAMS.gb
	cp out/pWJ40-PAMS.gb out/pWJ40-S1_L001_${i}-PAMS.gb
	cp out/NCTC8325-PAMS.gb out/NCTC8325-S1_L001_${i}-PAMS.gb
	time python2.7 -u ./align-spacers-mp.py spacers_S1_L001_${i}.fastq out/phiNM4-S1_L001_${i}-PAMS.gb out/pWJ40-S1_L001_${i}-PAMS.gb out/NCTC8325-S1_L001_${i}-PAMS.gb
done

