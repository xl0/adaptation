#!/usr/bin/env python2.7

import sys
import gzip
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import generic_dna

import numpy as np
import matplotlib.pyplot as plt

from utils import *

def usage():
	print "Usage:"
	print "\t" + sys.argv[0] + " <file1.fastq[.gz]> <file2.fastq[.gz]>"

def open_maybe_gzip(filename):
	fh = gzip.open(filename)
	try:
		tmp = fh.readline()
	except IOError:
		fh.close()
		fh = open(filename)

		return fh

	fh.close();
	fh = gzip.open(filename)

	return fh


if (len(sys.argv) != 3) :
	usage()
	exit(0)

fh1 = open_maybe_gzip(sys.argv[1])
seqh1 = SeqIO.parse(fh1, 'fastq', generic_dna)

fh2 = open_maybe_gzip(sys.argv[2])
seqh2 = SeqIO.parse(fh2, "fastq", generic_dna)

primers = load_sequences("primers.fasta")
tags = load_sequences("tags.fasta")

repeat = SeqRecord(Seq("GTTTTAGAGCTATGCTGTTTTGAATGGTCCCAAAAC", generic_dna), id="Repeat");

wfh = open("test.gb", "a");

print "Primers:"
for primer in primers:
	print primer.id, "\t", str(primer.seq)

print "Tags:"
for tag in tags:
	print str(tag.seq)


i = 0
nseq_primer = 0
for seq1 in seqh1:
	i = i + 1

	seq2 = seqh2.next()
	seq1.name = seq1.id = "R1_" + str(i)
	seq2.name = seq2.id = "R2_" + str(i)

	qs1 = seq1.letter_annotations["phred_quality"]
	avg1 = sum(qs1) / float(len(qs1))

	qs2 = seq2.letter_annotations["phred_quality"]
	avg2 = sum(qs2) / float(len(qs2))

#	print str(seq1.seq)
	annotate_primers(seq1, primers)
#	annotate_primers(seq2, primers)

	annotate_tags(seq1, tags)

	annotate_repeats(seq1, repeat)

	annotate_spacers(seq1)
#	show_features(seq1)

#	annotate_tags(seq2, tags)

#	print seq1.id
#	show_features(seq1)

#	SeqIO.write(seq1, wfh, "genbank");


#	if len(seq1.features) == 0:
#		print "!!!!!!!!!!!!!!!!"

	if i % 100 == 0:
		print i

wfh.close();



