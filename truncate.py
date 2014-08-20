#!/usr/bin/env python2.7

import sys
from Bio import SeqIO
from Bio.Seq import Seq
from utils import *


def usage():
	print "Copy a number of sequences from a file"
	print "Usage:"
	print "\t" + sys.argv[0] + " <input.fastq> <output.fastq> <number>"

if (len(sys.argv) != 4) :
	usage();
	exit(0);

number = float(sys.argv[3])


ifh = open_maybe_gzip(sys.argv[1]);
iseqh = SeqIO.parse(ifh, 'fastq');

ofh = open(sys.argv[2], "wr+");

i = 0;
for seq in iseqh:
	i = i + 1;

	if not (i % 100000):
		print i;

	if (i >= number):
		break;

	SeqIO.write(seq, ofh, "fastq");

ofh.close();
iseqh.close();
ifh.close();

