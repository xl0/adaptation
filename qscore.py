#!/usr/bin/env python2.7

import sys
import gzip
from Bio import SeqIO
from Bio.Seq import Seq
import numpy as np
import matplotlib.pyplot as plt


def usage():
	print "Usage:"
	print "\t" + sys.argv[0] + " <file.fastq>"

if (len(sys.argv) != 2) :
	usage();
	exit(0);

fh = gzip.open(sys.argv[1], "r");
seqh = SeqIO.parse(fh, 'fastq');


averages = [];
i = 0;

for seq in seqh:
	i = i + 1;

	qscore = seq.letter_annotations["phred_quality"];
	average = sum(qscore) / float(len(qscore))

	averages.append(average);
	if not (i % 100000):
		print i;


plt.hist(averages, range(1,1000))
plt.show();










