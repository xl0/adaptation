#!/usr/bin/env python2.7

import sys
import gzip
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import generic_dna


from utils import *
from sequences import *

import cProfile
import pstats


def usage():
	print "Usage:"
	print "\t" + sys.argv[0] + " <read1.fastq[.gz]> <read2.fastq[.gz]> <output.fastq>"

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

def main(seqh1, seqh2, out_fh):
	primers = get_primers()
	tags = get_tags()
	repeat = get_repeat()

	print "Primers:"
	for primer in primers:
		print primer.id, "\t", str(primer.seq)

	print "Tags:"
	for tag in tags:
		print str(tag.seq)

	s = i = 0
	nseq_primer = 0
	for seq1 in seqh1:
		i = i + 1

		seq2 = seqh2.next()
		seq1.name = seq1.id = "R1_" + str(i)
		seq2.name = seq2.id = "R2_" + str(i)

#		print str(seq1.seq)
		annotate_primers(seq1, primers)
	#	annotate_primers(seq2, primers)

		annotate_tags(seq1, tags)
		annotate_repeats(seq1, repeat)
		annotate_spacers(seq1)
#		show_features(seq1)

		spacers = extract_spacers(seq1)

		for spacer in spacers:
#			print str(spacer.seq)
#			print sum(spacer.letter_annotations["phred_quality"]) / float(len(spacer.letter_annotations["phred_quality"]))
			s = s + 1
			SeqIO.write(spacer, out_fh, "fastq")


		if i % 1000 == 0:
			print i

		if i > 10000:
			print s
			return


if (len(sys.argv) != 4):
	usage()
	exit(0)

fh1 = open_maybe_gzip(sys.argv[1])
seqh1 = SeqIO.parse(fh1, 'fastq', generic_dna)

fh2 = open_maybe_gzip(sys.argv[2])
seqh2 = SeqIO.parse(fh2, "fastq", generic_dna)

fh3 = open(sys.argv[3], "a+")
fh3.truncate()

main(seqh1, seqh2, fh3)

fh3.close()
fh2.close()
fh1.close()


#cProfile.run('main(seqh1, seqh2)', 'bzz')
#p = pstats.Stats('bzz')
#p.sort_stats('cumtime')
#p.print_stats()






