#!/usr/bin/env python2.7

import sys
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import generic_dna


from utils import *

def usage():
	print "Usage:"
	print "\t" + sys.argv[0] + "<input.gb[.gz]> <output.gb>"

def main(seqh, out):

	# We assume a single sequence
	seq = seqh.next()


	positions = []
	# We miss a single PAM if GG is at the circle split.
	start = 0
	while True:
		start = seq.seq.find("GG", start)
		if start == -1:
			break

		positions.append((start - 1, 1))
		start = start + 1

	# And for the - strand
	start = 0
	while True:
		start = seq.seq.find("CC", start)
		if start == -1:
			break

		positions.append((start, -1))
		start = start + 1

	positions = sorted(positions, key = lambda position: position[0])
	for position in positions:
		seq.features.append(SeqFeature(FeatureLocation(position[0], position[0] + 3),
					type = "PAM", strand = position[1]))

	SeqIO.write(seq, out, "genbank")

if (len(sys.argv) != 3):
	usage()
	exit(0)

fh1 = open_maybe_gzip(sys.argv[1])
seqh = SeqIO.parse(fh1, 'genbank', generic_dna)

fh2 = open(sys.argv[2], "wr+")
fh2.truncate()


print "Annotating ", fh1.name, " to ", fh2.name
main(seqh, fh2)

fh2.close()
fh1.close()

