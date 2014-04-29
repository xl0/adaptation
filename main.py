#!/usr/bin/env python2.7

import sys
from Bio import SeqIO
from Bio.Seq import Seq
import numpy as np
import matplotlib.pyplot as plt

from utils import *


def usage():
	print "Usage:"
	print "\t" + sys.argv[0] + " <file1.fastq> <file2.fastq>"

if (len(sys.argv) != 3) :
	usage();
	exit(0);

fh1 = open(sys.argv[1], "r");
seqh1 = SeqIO.parse(fh1, 'fastq');

fh2 = open(sys.argv[2], "r");
seqh2 = SeqIO.parse(fh2, "fastq");

primers = load_primers("primers.fasta");

while True:
	seq1 = seqh1.next();
	seq2 = seqh2.next();

	annotate_primers(seq1, primers);
	annotate_primers(seq2, primers)







