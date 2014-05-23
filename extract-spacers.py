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

	print "Repeat:"
	print str(repeat.seq)


	out_files = {}
	for tag in tags:
		out_files[tag.seq] = open(out_file_base + "_" + tag.seq + ".fastq", "a+")
		out_files[tag.seq].truncate()
	out_file_all = open(out_file_base + ".fastq", "a+")
	out_file_all.truncate()

	out_file_none = open(out_file_base + "_none.fastq", "a+")
	out_file_none.truncate()

	# XXX Add stats
	stats = {}
	stats["tags"] = tags
	stats["primers"] = primers

	q = []
	s = i = 0
	nseq_primer = 0

	multi_tag1 = 0
	mismatch1 = 0
	multi_tag2 = 0
	mismatch2 = 0

	global_mismatch = 0
	total_mismatch = 0

	for seq1 in seqh1:
		i = i + 1
		if i % 1000 == 0:
			print i

		if i > 10000:
			break

		seq2 = seqh2.next()
		seq1.name = seq1.id = "R1_" + str(i)
		seq2.name = seq2.id = "R2_" + str(i)

#		print str(seq1.seq)
		annotate_primers(seq1, primers)
		annotate_primers(seq2, primers)

		found_tags1 = annotate_tags(seq1, tags)
		found_tags2 = annotate_tags(seq2, tags)
		if len(found_tags1) > 1:
			multi_tag1 += 1

		annotate_repeats(seq1, repeat)
		annotate_repeats(seq2, repeat)
		annotate_spacers(seq1)
		annotate_spacers(seq2)

		# Unique tags
		utags1 = list(set(found_tags1))
		if len(utags1) > 1:
			mismatch1 += 1

		if len(found_tags2) > 1:
			multi_tag2 += 1

		utags2 = list(set(found_tags2))
		if len(utags2) > 1:
			mismatch2 += 1

		if len(utags1) == 1 and len(utags2) == 1 and not utags1[0] == utags2[0]:
			global_mismatch += 1

		if not len(utags1) == 1 or not len(utags2) == 1 or not utags1[0] == utags2[0]:
			total_mismatch += 1
			continue

		spacers1 = extract_spacers(seq1, utags1[0])
		spacers2 = extract_spacers(seq2, utags2[0])

		spacers = spacers1 + spacers2
		if len(spacers) > 3:
			print "WTF? %d spacers?" % len(spacers)
			show_sequence(seq1)
			show_sequence(seq2)

		for spacer in spacers:
			q.append(sum(spacer.letter_annotations["phred_quality"]) / float(len(spacer.letter_annotations["phred_quality"])))
			SeqIO.write(spacer, out_files[utags1[0].seq], "fastq")

	print "R1: Multi-tag: %d, Tag mismatch: %d (%d %%)" % (multi_tag1, mismatch1, (mismatch1 * 100) / multi_tag1)
	print "R2: Multi-tag: %d, Tag mismatch: %d (%d %%)" % (multi_tag2, mismatch2, (mismatch2 * 100) / multi_tag2)

	print "Global mismatches: %d, Total mismatches: %d (%d %%)" % (global_mismatch, total_mismatch, (total_mismatch * 100) / i)

	for fh in out_files.values():
		fh.close()
	out_file_all.close()
	out_file_none.close()

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

