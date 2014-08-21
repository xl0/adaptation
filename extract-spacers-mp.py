#!/usr/bin/env python2.7

import sys
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import generic_dna
from collections import defaultdict

from utils import *
from sequences import *

from multiprocessing import Pool

import json

import cProfile
import pstats

# We run into IO bottlenech at around 6
num_threads = 16

def usage():
	print "Usage:\t" + sys.argv[0] + " <read1.fastq[.gz]> <read2.fastq[.gz]> <out_file_base>"

class DoubleSeqIterable:
	def __init__(self, seqh1, seqh2):
		self.seqh1 = seqh1
		self.seqh2 = seqh2
		self.n = 0

	def next(self):
		tmp = (self.seqh1.next(), self.seqh2.next(), self.n)
		self.n += 1

		return tmp

	def __iter__(self):
		return self

global g_tags
global g_repeat
global g_primers

def worker_function(arg):
	seq1 = arg[0]
	seq2 = arg[1]
	i = arg[2]
	tags = g_tags
	primers = g_primers
	repeat = g_repeat


	seq1.name = seq1.id = "R1_" + str(i)
	seq2.name = seq2.id = "R2_" + str(i)

	found_primers1 = annotate_primers(seq1, primers)
	found_primers2 = annotate_primers(seq2, primers)

	found_tags1 = annotate_tags(seq1, tags)
	found_tags2 = annotate_tags(seq2, tags)

	# Unique tags
	utags1 = list(set(found_tags1))
	utags2 = list(set(found_tags2))

	# If we get an abnormal number of tags, or the tags do not match, skip this pair.
	if not len(utags1) == 1 or not len(utags2) == 1 or not utags1[0] == utags2[0]:
		return None

	annotate_repeats(seq1, repeat)
	annotate_repeats(seq2, repeat)
	annotate_spacers(seq1)
	annotate_spacers(seq2)



	spacers = extract_spacers(seq1, utags1[0]) + extract_spacers(seq2, utags2[0])

#	if utags1[0] == 'GAC' and not len(spacers) == 2:
#		print '======='
#		show_sequence(seq1)
#		show_sequence(seq2)
		
	if len(spacers) > 2:
		print "\nWTF? %d spacers?" % len(spacers)
		show_sequence(seq1)
		show_sequence(seq2)


	# If we got exactly 2 spacers, see if they have the same sequence.
	if len(spacers) == 2:
		if len(spacers[0]) == len(spacers[1]):
			# Let's allow up to 3 mismatches, and no insertions/deletions.
			if(seq_mismatches(spacers[0], spacers[1]) <= 3):
				spacer = ambiguous_merge(spacers[0], spacers[1])
				return (utags1[0], spacer)

	# Good tag, but no spacer - still return tag for statistics
	return (utags1[0], None)

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
		out_files[tag.id] = open(sys.argv[3] + tag.id + ".fastq", "wr+")

	out_files['stats'] = open(sys.argv[3] + 'stats.json', 'wr+')
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




	global g_tags
	global g_repeat
	global g_primers
	def share_global_data(tags, primers, repeat):
		global g_tags
		global g_repeat
		global g_primers
		g_tags = tags
		g_primers = primers
		g_repeat = repeat

	iterable = DoubleSeqIterable(seqh1, seqh2)
	pool = Pool(initializer=share_global_data, initargs=(tags, primers, repeat,))

	stats = {	'num_sequences' : 0,
			'num_sequences_w_tag' : 0,
			'num_spacers' : 0,
			'spacer_lengths' : defaultdict(int),
		}
	for tag in tags:
		stats['num_sequences_' + tag.id] = 0
		stats['num_spacers_' + tag.id] = 0
		stats['spacer_lengths_' + tag.id] = defaultdict(int)

	stats['spacer_lengths'] = defaultdict(int)


	print "Starting pool workers in %d threads" % num_threads
	sys.stdout.write('Processing spacers... ')
	for res in pool.imap_unordered(worker_function, iterable, 1):
		if not stats['num_sequences'] % 1000000:
			sys.stdout.write("%d ..." % stats['num_sequences'])
	
		stats['num_sequences'] += 1
		if res:
			tag = res[0]
			stats['num_sequences_w_tag'] += 1
			stats['num_sequences_' + tag] += 1
			spacer = res[1]
			if spacer:
				length = len(spacer.seq)
				SeqIO.write(spacer, out_files[tag], "fastq")
				stats['num_spacers'] += 1
				stats['num_spacers_' + tag] += 1
				stats['spacer_lengths'][length] += 1
				stats['spacer_lengths_' + tag][length] += 1
	print stats['num_sequences']

	# Ignore sequences with tags that represent less than 1% of the tagged sequences.
	threshold = float(stats['num_sequences_w_tag']) / 100

	print 'Analyzed %d sequences, %d had valid tags (%.2f %%).' % (
		stats['num_sequences'],
		stats['num_sequences_w_tag'],
		float(stats['num_sequences_w_tag']) / stats['num_sequences'] * 100)
	print 'Breackdown by tag:'
	for tag in tags:
		num_seq = stats['num_sequences_' + tag.id]
		if num_seq > 0:
			if num_seq > threshold:
				print "%s: %d (%.2f %%)" % (tag.id, num_seq, float(num_seq) / stats['num_sequences_w_tag'] * 100)
			else:
				print "%s: %d (%.2f %%) (ignored)" % (tag.id, num_seq, float(num_seq) / stats['num_sequences_w_tag'] * 100)

	print 'Extracted %d spacers from %d tagged sequences (%.2f %%, %.2f %% of all sequences)' % (
		stats['num_spacers'],
		stats['num_sequences_w_tag'],
		float(stats['num_spacers']) / stats['num_sequences_w_tag'] * 100,
		float(stats['num_spacers']) / stats['num_sequences'] * 100)
	print 'Breackdown by tag:'
	for tag in tags:
		num_seq = stats['num_sequences_' + tag.id]
		if num_seq > threshold:
			print "%s: %d (%.2f %% of spacers, %.2f %% sequences with this tag had good spacers)" % (
				tag.id,
				stats['num_spacers_' + tag.id],
				float(stats['num_spacers_' + tag.id]) / stats['num_spacers'] * 100,
				float(stats['num_spacers_' + tag.id]) / stats['num_sequences_' + tag.id] * 100)

	json.dump(stats, out_files['stats'], sort_keys=True, indent=4)

	fh1.close()
	fh2.close()
	for tag in tags:
		num_seq = stats['num_sequences_' + tag.id]
		if num_seq < threshold:
			os.unlink(out_files[tag.id].name)
		out_files[tag.id].close()
	
	out_files['stats'].close()	

if __name__ == "__main__":
    main()

#cProfile.run('main(seqh1, seqh2, sys.argv[3])', 'bzz')
#p = pstats.Stats('bzz')
#p.sort_stats('cumtime')
#p.print_stats()


