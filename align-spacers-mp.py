#!/usr/bin/env python2.7

import sys
from Bio import SeqIO
from Bio.Seq import Seq, reverse_complement
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import generic_dna

from utils import *
from collections import defaultdict, OrderedDict

import os.path

import traceback

from multiprocessing import Pool, sharedctypes
#import matplotlib.pyplot as plt
import cProfile
import pstats

def usage():
	print "Usage:"
	print "\t" + sys.argv[0] + " <spacers.fastq[.gz]> <pamdb1.gb> [pamdb2.gb ...]"

# found_spacer dict:
# key: spacer sequence
# value: (template_name, (same as key in the spacer_db))
found_spacer_dict = defaultdict(list)

# spacer_db
# key: template file name
# value: dictionary
#	key: (start, strand, length, mismatch, seq str)
#	value: number of hits
spacer_db = {}


class DoubleSeqIterable:
	def __init__(self, seqh):
		self.seqh = seqh
		self.n = 0

	def next(self):
		# Check if the primer was already found.		

		while True:
			seq = self.seqh.next()
			str_seq = str(seq.seq)
			self.n += 1
			if not self.n % 100000:
				print self.n

#			if self.n > 100000:
#				raise StopIteration

			if found_spacer_dict.has_key(str_seq):
				spacer_dict_key_list = found_spacer_dict[str_seq]
				if spacer_dict_key_list[0]:
					for spacer_dict_key in spacer_dict_key_list[1]:
						spacer_db[spacer_dict_key[0]][spacer_dict_key[1]] += 1

			else:
				break

		return (str_seq, self.n - 1)

	def __iter__(self):
		return self

global shared_template_db
def real_worker_function(arg):
	seq = arg[0]
	template_db = shared_template_db
	n = arg[1]

	match = 0
	mismatch = 0


	# Align spacers to the DB sequence, counting the number of occurance

	spacer_db = defaultdict(int)
	for key in shared_template_db.keys():
		spacer_db[key] = defaultdict(int)

	found = False
	for key in shared_template_db.keys():
		# Try a simple search - it's much faster than alignment.
		db_seq = shared_template_db[key][0]
		spacers = spacer_db[key]
		location = 0
		while True:
			location = db_seq.seq.find(seq, location)
			if location == -1:
				break

			spacers[(location, 1, len(seq), 0)] += 1
			found = True
			location += 1
			match += 1
		location = 0
		while True:
			location = db_seq.seq.find(reverse_complement(seq), location)
			if location == -1:
				break

			spacers[(location, -1, len(seq), 0)] += 1
			found = True
			location += 1
			match += 1

		# Found something? No need to align. Also don't align large sequences - it's too slow.
		if found:
			break

		if len(db_seq) > 50000:
			continue

		score, pos_db, pos_seq = align(seq, str(key) + ".seq")

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

			spacers[(pos_db[0] - head_correction, 1, length , penalty / 5, seq)] += 1


		score, pos_db, pos_seq = align(reverse_complement(seq), str(key) + ".seq")
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

			spacers[(pos_db[0] - head_correction, -1, length , penalty / 5, seq)] += 1

	return (found, seq, spacer_db, n)

def worker_function(arg):

	try:
		return real_worker_function(arg)
	except:
		print "".join(traceback.format_exception(*sys.exc_info()))
		raise RuntimeError

def gen_tmp_sequence(db_seq, db_fd):
	name = db_fd.name + ".seq"

	fh = open(name, "wr+")
	fh.write(str(db_seq.seq))
	fh.write("\n")
	fh.close()

def main(seqh, dbs):


	# template_db
	# key: template file name
	# value: (db_seq, db_fd)
	template_db = OrderedDict()


	match = 0
	mismatch = 0
	for db_fd in dbs:
		print "Parsing ", db_fd.name
		db_seqh = SeqIO.parse(db_fd, "genbank")
		db_seq = db_seqh.next()
		db_seqh.close()

		if db_fd.name in template_db.keys():
			print "Duplicate DB: ", db_name
			return

		# Remove spacer information from the DB to avoid duplicating it.
		features = []
		for feature in db_seq.features:
			if not feature.type == "Spacer":
				features.append(feature)

		db_seq.features = features

		template_db[db_fd.name] = (db_seq, db_fd)
		spacer_db[db_fd.name] = defaultdict(int)

		gen_tmp_sequence(db_seq, db_fd)


	print "Aligning spacer sequences..."

	iterable = DoubleSeqIterable(seqh)

	global shared_template_db
	def share_template_db(db):
		global shared_template_db
		shared_template_db = db

	pool = Pool(initializer=share_template_db, initargs=(template_db,))

	for res in pool.imap(worker_function, iterable, 10):
		found = res[0]
		spacer_seq = res[1]
		spacer_db_update = res[2]
		n = res[3]

		if found:
			update_list = []

			# For each template
			for template in spacer_db_update.keys():
				spacer_dict_update = spacer_db_update[template]
				# copy the found spacers to the central DB
				for spacer_key in spacer_dict_update.keys():
					update_list.append((template, spacer_key))

			for update in update_list:
				spacer_db[update[0]][update[1]] += 1

			found_spacer_dict[spacer_seq] = (True, update_list)
		else:
			found_spacer_dict[spacer_seq] = (False, None)

	print "Updating the database files..."

	# Add spacer annotations to the db sequence
	for template_key  in template_db.keys():
		print "Updating %s ..." % template_key
		db_seq = template_db[template_key][0]
		db_fd = template_db[template_key][1]
		db_spacers = spacer_db[template_key]

		for spacer_key in db_spacers.keys():
			qualifiers = { 'hits' : db_spacers[spacer_key],
				       'len' : spacer_key[2]}

			if spacer_key[3]:
				qualifiers['mismatches'] = spacer_key[3]
				qualifiers['sequence'] = spacer_key[4]

			db_seq.features.append(SeqFeature(FeatureLocation(spacer_key[0], spacer_key[0] + spacer_key[2]),
				type="Spacer", strand=spacer_key[1], qualifiers=qualifiers))
		db_name = db_fd.name
		db_fd.close()

		db_fd = open(db_name, "wr+")
		SeqIO.write(db_seq, db_fd, "genbank")
		db_fd.close()


if (len(sys.argv) < 3):
	usage()
	exit(0)

seqfh = open_maybe_gzip(sys.argv[1])
seqh = SeqIO.parse(seqfh, 'fastq', generic_dna)

dbs = []

#print sys.argv
for db in sys.argv[2:]:
	print db
	fh = open(db, "rw")
	dbs.append(fh) #db_seqh.next())

main(seqh,dbs)
#cProfile.run('main(seqh, dbs)', 'bzz')
#p = pstats.Stats('bzz')
#p.sort_stats('cumtime')
#p.print_stats()

#for db in dbs:
#	db.close()

