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

# We run into IO bottlenech at around 6
num_threads = 16

def usage():
	print "Usage:\t" + sys.argv[0] + " <read1.fastq[.gz]> <read2.fastq[.gz]> <out_file_base>"

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

	found_primers1 = annotate_primers(seq1, primers)
	found_primers2 = annotate_primers(seq2, primers)

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
				return (spacer,  utags1[0])

	return None


def main():
	if (len(sys.argv) != 4):
		usage()
		exit(0)

	fh1 = open_maybe_gzip(sys.argv[1])
	seqh1 = SeqIO.parse(fh1, 'fastq', generic_dna)

	fh2 = open_maybe_gzip(sys.argv[2])
	seqh2 = SeqIO.parse(fh2, "fastq", generic_dna)

	tags = get_tags()
	out_files = {}
	for tag in tags:
		out_files[tag.seq] = open(sys.argv[3] + tag.seq + ".fastq", "wr+")

	primers = get_primers()
	repeat = get_repeat()


	print "Extrating spacer information from %s and %s..."  % (fh1.name, fh2.name)
	print "Primers:"
	for primer in primers:
		print primer.id, "\t", str(primer.seq)

	print "Tags:"
	for tag in tags:
		print str(tag.seq)

	print "Repeat:"
	print str(repeat.seq)



	iterable = DoubleSeqIterable(seqh1, seqh2, tags, primers, repeat)

	pool = Pool(num_threads)

	print "Starting pool workers in %d threads" % num_threads

	num_spacers = 0
	num_processed = 0
	for res in pool.imap_unordered(worker_function, iterable, 100):
		if num_processed == 0:
			sys.stdout.write('Processing spacers... ')
			sys.stdout.flush()
		if not num_processed % 1000000:
			sys.stdout.write("%d ..." % num_processed)
			sys.stdout.flush()

		num_processed += 1
		if res:
			spacer = res[0]
			tag = res[1]
			SeqIO.write(spacer, out_files[tag.seq], "fastq")
			num_spacers += 1
	sys.stdout.write('%d\n' %  num_processed)

	print "Extracted %d spacers from %d sequences" % (num_spacers, iterable.n)
	fh1.close()
	fh2.close()
	for fh in out_files.values():
		pos = fh.tell()
		# Nothing in this file - remove it.
		if pos == 0:
			os.unlink(fh.name)
			
		fh.close()


if __name__ == "__main__":
    main()

#cProfile.run('main(seqh1, seqh2, sys.argv[3])', 'bzz')
#p = pstats.Stats('bzz')
#p.sort_stats('cumtime')
#p.print_stats()


