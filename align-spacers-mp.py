#!/usr/bin/env python2.7

import sys
import os.path
import argparse

from Bio import SeqIO
from Bio.Seq import Seq, reverse_complement
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import generic_dna

from utils import *
from collections import defaultdict, OrderedDict
import traceback
from multiprocessing import Pool

import json

def usage():
	print "Usage:"
	print "\t" + sys.argv[0] + " <spacers.fastq[.gz]> <pamdb1.gb> [pamdb2.gb ...]"

# found_spacer dict:
# key: spacer sequence
# value: (template_name, (same as key in the spacer_db))
found_spacer_dict = {}

# spacer_db dictionary
#	key: (template, start, strand, length, mismatch, [seq str])
#	value: number of hits
spacer_db = defaultdict(int)

class DoubleSeqIterable:
	def __init__(self, seqh):
		self.seqh = seqh
		self.n = 0
		self.hits = 0
		self.misses = 0
		self.last_hits = 0
		self.last_misses = 0

	def next(self):
		# Check if the primer was already found.

		while True:
			seq = self.seqh.next()
			str_seq = str(seq.seq)
			self.n += 1
			if not self.n % 10000:
				print self.n, self.hits, self.misses, self.last_hits, self.last_misses
				self.last_hits = 0
				self.last_misses = 0

#			if self.n > 100000:
#				raise StopIteration

			if found_spacer_dict.has_key(str_seq):
				self.hits += 1
				self.last_hits += 1
				# Sequence does not align -> None
				if found_spacer_dict[str_seq]:
					for alignment in found_spacer_dict[str_seq]:
						spacer_db[tuple(alignment)] += 1

			else:
				self.misses += 1
				self.last_misses += 1
				break

		return (str_seq, self.n - 1)

	def __iter__(self):
		return self
	def get_hits(self):
		return self.hits
	def get_misses(self):
		return self.misses


global shared_template_db
def real_worker_function(arg):
	seq = arg[0]
	template_db = shared_template_db

	match = 0
	mismatch = 0

	# Ignore spacers that are obviously too short
	if len(seq) < 15:
		return (False, seq, None, None)

	# Align spacers to the DB sequence, counting the number of occurance

	found = False
	alignments = []
	for template in shared_template_db.keys():
		# Try a simple search - it's much faster than alignment.
		db_seq = shared_template_db[template][0]
		db_file_name = shared_template_db[template][1]
		location = 0
		while True:
			location = db_seq.seq.find(seq, location)
			if location == -1:
				break

			alignments.append((template, location, 1, len(seq), 0))
			found = True
			location += 1
			match += 1
		location = 0
		while True:
			location = db_seq.seq.find(reverse_complement(seq), location)
			if location == -1:
				break

			alignments.append((template, location, -1, len(seq), 0))
			found = True
			location += 1
			match += 1

		# Found something? No need to align. Also don't align large sequences - it's too slow.
		if found:
			break

		if len(db_seq) > 50000:
			continue

		score, pos_db, pos_seq = align(seq, db_file_name + ".seq")

		# Score is 5 * len, -5 for a mismatch, -10 for a gap
		# Let's allow up to 4 mismatches, or 2 gaps (20 points penalty)

		penalty = len(seq) * 5 - score

		if (penalty <= 20):
			# Acceptable match. Make sure we account for the firsr and last BP(s) in
			# case the mismatch was there.
			head_correction = -pos_seq[0]
			tail_correction = len(seq) - pos_seq[1]
			mismatch += 1

			found = True
			length = pos_db[1] + tail_correction - (pos_db[0] - head_correction)

			alignments.append((template, pos_db[0] - head_correction, 1, length , penalty / 5, seq))


		score, pos_db, pos_seq = align(reverse_complement(seq), db_file_name + ".seq")
		# Score is 5 * len, -5 for a mismatch, -10 for a gap
		# Let's allow up to 4 mismatches, or 2 gaps (20 points penalty)

		penalty = len(seq) * 5 - score

		if (penalty <= 20):
			# Acceptable match. Make sure we account for the firsr and last BP(s) in
			# case the mismatch was there.
			head_correction = -pos_seq[0]
			tail_correction = len(seq) - pos_seq[1]
			mismatch += 1

			found = True
			length = pos_db[1] + tail_correction - (pos_db[0] - head_correction)

			alignments.append((template, pos_db[0] - head_correction, -1, length , penalty / 5, seq))

		if found:
			break

	return (found, seq, alignments)

def worker_function(arg):

	try:
		return real_worker_function(arg)
	except:
		print "".join(traceback.format_exception(*sys.exc_info()))
		raise RuntimeError

def gen_tmp_sequence(template_seq, template_file):
	name = template_file + ".seq"

	fh = open(name, "wr+")
	fh.write(str(template_seq.seq))
	fh.write("\n")
	fh.close()

	return name

def main(argv):

	parser = argparse.ArgumentParser(description='Align spacers to templates')
	parser.add_argument('spacers', metavar='<spacer_file.fastq[.gz]>', type=str,
				help='Input file with PAMs and spacers annotated.')
	parser.add_argument('templates', metavar='<template.gb>', type=str, nargs='+',
				help='Template sequences')
	parser.add_argument('-j', metavar='<out_file.json>', type=str,
				help='Write statistics to a json file')
	parser.add_argument('-c', metavar='<cache.json>', type=str,
				help='Spacer alignment cache file')
	args = parser.parse_args()

	print 'Analyzing', args.spacers

	seqfh = open_maybe_gzip(args.spacers)
	seqh = SeqIO.parse(seqfh, 'fastq', generic_dna)

	# template_db
	# key: template file name
	# value: (db_seq, db_fd)
	template_db = OrderedDict()

	# List of temporary files, fed to water, and to be removed later.
	template_tmp_seq_list = []

	match = 0
	mismatch = 0
	for template_file in args.templates:
		template_name = os.path.basename(template_file)
		print "Parsing ", template_name

		fd = open(template_file)
		template_seqh = SeqIO.parse(fd, "genbank")
		template_seq = template_seqh.next()
		template_seqh.close()
		fd.close()

		if template_name in template_db.keys():
			print "Duplicate DB: ", template_name
			return 1

		# Remove spacer information from the DB to avoid duplicating it.
		features = []
		for feature in template_seq.features:
			if not feature.type == "Spacer":
				features.append(feature)

		template_seq.features = features

		template_db[template_name] = (template_seq, template_file)

		temp_name = gen_tmp_sequence(template_seq, template_file)
		template_tmp_seq_list.append(temp_name)

	if args.c:
		if os.path.isfile(args.c):
			cache_f = open(args.c)
			global found_spacer_dict
			found_spacer_dict = json.load(cache_f)
			cache_f.close()

			print 'Loaded %s entries from alignment cache %s' % (len(found_spacer_dict), args.c)

	print "Aligning spacer sequences from %s to %s." % (seqfh.name, ' '.join([x for x in template_db.keys()]))

	iterable = DoubleSeqIterable(seqh)

	global shared_template_db
	def share_template_db(db):
		global shared_template_db
		shared_template_db = db

	pool = Pool(initializer=share_template_db, initargs=(template_db,))

	for res in pool.imap(worker_function, iterable, 10):
		found = res[0]
		spacer_seq = res[1]
		alignments = res[2]

		if found:
			for alignment in alignments:
				spacer_db[tuple(alignment)] += 1
			found_spacer_dict[spacer_seq] = alignments
		else:
			found_spacer_dict[spacer_seq] = None


	print "Spacers analyzed: ", iterable.n
	print "Cache hits: ", iterable.get_hits()
	print "Cache misses: ", iterable.get_misses()

	print "Processing collected data ..."

	for alignment in spacer_db.keys():
		template = alignment[0]
		template_seq = template_db[template][0]

		hits = spacer_db[alignment]
		start = alignment[1]
		strand = alignment[2]
		length = alignment[3]

		qualifiers = {
			'hits' : hits,
			'len' : length
		}

		mismatches = alignment[4]
		if mismatches > 0:
			sequence = alignment[5]
			qualifiers['mismatches'] = mismatches
			qualifiers['sequence'] = sequence

		template_seq.features.append(SeqFeature(FeatureLocation(start, start+length),
			type="Spacer", strand=strand, qualifiers=qualifiers))

	for template in template_db.values():
		template_seq = template[0]
		template_filename = template[1]

		print "Updating %s ..." % template_filename

		fd = open(template_filename, 'wr+')
		SeqIO.write(template_seq, fd, 'genbank')
		fd.close()

	if args.c:
		print 'Updating alignment cache %s' % args.c
		cache_f = open(args.c, 'wr+')
		json.dump(found_spacer_dict, cache_f, sort_keys=True, indent=4)
		cache_f.close()

	for tmp in template_tmp_seq_list:
		os.unlink(tmp)

if __name__ == "__main__":
    main(sys.argv)


#cProfile.run('main(seqh, dbs)', 'bzz')
#p = pstats.Stats('bzz')
#p.sort_stats('cumtime')
#p.print_stats()

#for db in dbs:
#	db.close()

