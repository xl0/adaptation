#!/usr/bin/env python2.7

import sys
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import generic_dna


from utils import *
from sequences import *

from multiprocessing import Pool


import cProfile
import pstats


def usage():
	print "Usage:"
	print "\t" + sys.argv[0] + " <read1.fastq[.gz]> <read2.fastq[.gz]> <out_file_base>"

class DoubleSeqIterable:
	def __init__(self, seqh1, seqh2, tags, primers, repeat):
		self.seqh1 = seqh1
		self.seqh2 = seqh2
		self.tags = tags
		self.primers = primers
		self.repeat = repeat
		self.n = 0

	def next(self):
		tmp = (self.seqh1.next(), self.seqh2.next(), self.tags, self.primers, self.repeat, self.n)
		self.n += 1

		if not self.n % 10000:
			print self.n
		if self.n > 1000000:
			raise StopIteration
		return tmp

	def __iter__(self):
		return self

def worker_function(arg):
	seq1 = arg[0]
	seq2 = arg[1]
	tags = arg[2]
	primers = arg[3]
	repeat = arg[4]
	i = arg[5]


	seq1.name = seq1.id = "R1_" + str(i)
	seq2.name = seq2.id = "R2_" + str(i)

	annotate_primers(seq1, primers)
	annotate_primers(seq2, primers)

	found_tags1 = annotate_tags(seq1, tags)
	found_tags2 = annotate_tags(seq2, tags)

	annotate_repeats(seq1, repeat)
	annotate_repeats(seq2, repeat)
	annotate_spacers(seq1)
	annotate_spacers(seq2)

	# Unique tags
	utags1 = list(set(found_tags1))

	utags2 = list(set(found_tags2))

	# If we get an abnormal number of tags, or the tags do not match, skip this pair.
	if not len(utags1) == 1 or not len(utags2) == 1 or not utags1[0] == utags2[0]:
		return None

	spacers = extract_spacers(seq1, utags1[0]) + extract_spacers(seq2, utags2[0])

	if len(spacers) > 2:
		print "WTF? %d spacers?" % len(spacers)
		show_sequence(seq1)
		show_sequence(seq2)

	# If we got exactly 2 spacers, see if they have the same sequence.
	if len(spacers) == 2:
		if len(spacers[0]) == len(spacers[1]):
			# Let's allow up to 3 mismatches, and no insertions/deletions.
			if(seq_mismatches(spacers[0], spacers[1]) <= 3):
				spacer = ambiguous_merge(spacers[0], spacers[1])
				return (spacer, utags1[0])

	return None


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
		out_files[tag.seq] = open(out_file_base + "_" + tag.seq + ".fastq", "wr+")
		out_files[tag.seq].truncate()
	out_file_all = open(out_file_base + ".fastq", "wr+")
	out_file_all.truncate()

	out_file_none = open(out_file_base + "_none.fastq", "wr+")
	out_file_none.truncate()


	s = i = 0
	nseq_primer = 0

	n_perfect_spacer = 0


	iterable = DoubleSeqIterable(seqh1, seqh2, tags, primers, repeat)

	pool = Pool(6)

	print "Hello!"



	for res in pool.imap_unordered(worker_function, iterable, 100):
		if res:
			spacer = res[0]
			tag = res[1]
			SeqIO.write(spacer, out_files[tag.seq], "fastq")
			n_perfect_spacer += 1


	print "Bye!"

	print "i = ", iterable.n

#	print "R1: Multi-tag: %d, Tag mismatch: %d (%d %%)" % (multi_tag1, mismatch1, (mismatch1 * 100) / multi_tag1)
#	print "R2: Multi-tag: %d, Tag mismatch: %d (%d %%)" % (multi_tag2, mismatch2, (mismatch2 * 100) / multi_tag2)

#	print "Global mismatches: %d, Total mismatches: %d (%d %%)" % (global_mismatch, total_mismatch, (total_mismatch * 100) / i)
	print "Perfect spacers found: %d" % (n_perfect_spacer)

	for fh in out_files.values():
		fh.close()
	out_file_all.close()
	out_file_none.close()

if (len(sys.argv) != 4):
	usage()
	exit(0)

fh1 = open_maybe_gzip(sys.argv[1])
seqh1 = SeqIO.parse(fh1, 'fastq', generic_dna)

fh2 = open_maybe_gzip(sys.argv[2])
seqh2 = SeqIO.parse(fh2, "fastq", generic_dna)


main(seqh1, seqh2, sys.argv[3])
#cProfile.run('main(seqh1, seqh2, sys.argv[3])', 'bzz')
#p = pstats.Stats('bzz')
#p.sort_stats('cumtime')
#p.print_stats()

fh2.close()
fh1.close()


