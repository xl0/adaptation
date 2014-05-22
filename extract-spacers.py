#!/usr/bin/env python2.7

import sys
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import generic_dna


from utils import *
from sequences import *

import cProfile
import pstats

import matplotlib as mpl
import matplotlib.pyplot as plt



def usage():
	print "Usage:"
	print "\t" + sys.argv[0] + " <read1.fastq[.gz]> <read2.fastq[.gz]> <out_file_base>"

def main(seqh1, seqh2, out_file_base):
	primers = get_primers()
	tags = get_tags()
	repeat = get_repeat()

	print "Primers:"
	for primer in primers:
		print primer.id, "\t", str(primer.seq)

	print "Tags:"
	for tag in tags:
		print str(tag.seq)

	out_files = {}
	for tag in tags:
		out_files[tag.seq] = open(out_file_base + "_" + tag.seq + ".fastq", "a+")
		out_files[tag.seq].truncate()
	out_file_all = open(out_file_base + ".fastq", "a+")
	out_file_all.truncate()

	q = []
	s = i = 0
	nseq_primer = 0

	multi_tag = 0
	mismatch = 0

	for seq1 in seqh1:
		i = i + 1

		seq2 = seqh2.next()
		seq1.name = seq1.id = "R1_" + str(i)
		seq2.name = seq2.id = "R2_" + str(i)

#		print str(seq1.seq)
		annotate_primers(seq1, primers)
	#	annotate_primers(seq2, primers)

		found_tags = annotate_tags(seq1, tags)
		if len(found_tags) > 1:
			multi_tag = multi_tag + 1
		if len(found_tags) > 2:
			print "Wtf? %d tags?" % len(found_tags)
			show_sequence(seq1)

		utags = list(set(found_tags))
		if len(utags) > 1:
			mismatch = mismatch + 1
			

		annotate_repeats(seq1, repeat)
		annotate_spacers(seq1)
#		show_features(seq1)

		spacers = extract_spacers(seq1)

		for spacer in spacers:
#			print str(spacer.seq)
			q.append(sum(spacer.letter_annotations["phred_quality"]) / float(len(spacer.letter_annotations["phred_quality"])))
			SeqIO.write(spacer, out_files[tag.seq], "fastq")
			SeqIO.write(spacer, out_file_all, "fastq")
		
#		show_sequence(seq1)

		if i % 1000 == 0:
			print i

		if i > 100000:
			break

	print "Multi-tag: %d, Tag mismatch: %d (%d %%)" % (multi_tag, mismatch, (mismatch * 100) / multi_tag)

	for fh in out_files.values():
		fh.close()

#	plt.hist(q, range(1,1000))
#	plt.show()



if (len(sys.argv) != 4):
	usage()
	exit(0)

fh1 = open_maybe_gzip(sys.argv[1])
seqh1 = SeqIO.parse(fh1, 'fastq', generic_dna)

fh2 = open_maybe_gzip(sys.argv[2])
seqh2 = SeqIO.parse(fh2, "fastq", generic_dna)

#fh3 = open(sys.argv[3], "a+")
#fh3.truncate()

main(seqh1, seqh2, sys.argv[3])

#fh3.close()
fh2.close()
fh1.close()


#cProfile.run('main(seqh1, seqh2)', 'bzz')
#p = pstats.Stats('bzz')
#p.sort_stats('cumtime')
#p.print_stats()

