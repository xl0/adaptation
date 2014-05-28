#!/usr/bin/env python2.7

import sys
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import generic_dna

from utils import *
from collections import defaultdict

#import matplotlib.pyplot as plt
import cProfile
import pstats

def usage():
	print "Usage:"
	print "\t" + sys.argv[0] + " <spacers.fastq[.gz]> <pamdb1.gb> [pamdb2.gb ...]"

def main(seqh, dbs):

	fast_db = {}

	match = 0
	mismatch = 0
	for db_fd in dbs:
		print "Parsing ", db_fd.name
		db_seqh = SeqIO.parse(db_fd, "genbank")
		db_seq = db_seqh.next()
		db_seqh.close()

		if db_fd.name in fast_db.keys():
			print "Duplicate DB: ", db_name
			return

		# Remove spacer information from the DB to avoid duplicating it.

		features = []
		for feature in db_seq.features:
			if not feature.type == "Spacer":
				features.append(feature)

		db_seq.features = features


		# 3rd argument - a dict with a tuple as keys:
		# [0] - start
		# [1] - strand
		# [2] - length
		# [3] - number of mismatches
		# [4] - spacer sequence - in case of mismatch
		fast_db[db_fd.name] = [db_seq, db_fd, defaultdict(int)]

	print "Aligning spacer sequences..."

	# Align spacers to the DB sequence, counting the number of occurance
	nseq = 0
	for seq in seqh:
		nseq += 1
		if not nseq % 100:
			print nseq

		found = 0
		for db in fast_db.values():
			# Try a simple search - it's much faster than alignment.
			db_seq = db[0]
			db_spacers = db[2]
			location = 0
			while True:
				location = db_seq.seq.find(str(seq.seq), location)
				if location == -1:
					break

				db_spacers[(location, 1, len(seq), 0)] += 1
				found = 1
				location += 1
				match += 1
			location = 0
			while True:
				location = db_seq.seq.find(str(seq.seq.reverse_complement()), location)
				if location == -1:
					break

				db_spacers[(location, -1, len(seq), 0)] += 1
				found = 1
				location += 1
				match += 1

			# Found something? No need to align.
			if found or len(db_seq) > 50000:
				break

			score, pos_db, pos_seq = align(seq.seq, db_seq.seq)

			# Score is 5 * len, -5 for a mismatch, -10 for a gap
			# Let's allow up to 4 mismatches, or 2 gaps (20 points penalty)

			penalty = len(seq) * 5 - score

			if (penalty <= 20):
				# Acceptable match. Make sure we account for the firsr and last BP(s) in
				# case the mismatch was there.


				head_correction = -pos_seq[0]
				tail_correction = len(seq) - pos_seq[1]
				mismatch += 1

				length = pos_db[1] + tail_correction - (pos_db[0] - head_correction)

				db_spacers[(pos_db[0] - head_correction, 1, length , penalty / 5, str(seq.seq))] += 1


			score, pos_db, pos_seq = align(seq.seq.reverse_complement(), db_seq.seq)
			# Score is 5 * len, -5 for a mismatch, -10 for a gap
			# Let's allow up to 4 mismatches, or 2 gaps (20 points penalty)

			penalty = len(seq) * 5 - score

			if (penalty <= 20):
				# Acceptable match. Make sure we account for the firsr and last BP(s) in
				# case the mismatch was there.


				head_correction = -pos_seq[0]
				tail_correction = len(seq) - pos_seq[1]
				mismatch += 1

				length = pos_db[1] + tail_correction - (pos_db[0] - head_correction)

				db_spacers[(pos_db[0] - head_correction, -1, length , penalty / 5, str(seq.seq))] += 1


	print "%d spacers found out of %d sequences (%d %%)" % (match + mismatch, nseq, (mismatch + match) * 100 / nseq)
	print "%d spacers matched perfectly (%d %%)" % (match, match * 100 / nseq)
	print "%d spacers matched partially (%d %%)" % (mismatch, mismatch * 100 / nseq)

	print "Updating the database files..."

	# Add spacer annotations to the db sequence
	for db_entry in fast_db.values():
		print "Updating %s ..." % db_fd.name
		db_seq = db_entry[0]
		db_fd = db_entry[1]
		db_spacers = db_entry[2]

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

		db_fd = open(db_name, "a+")
		db_fd.truncate()
		SeqIO.write(db_seq, db_fd, "genbank")


if (len(sys.argv) < 3):
	usage()
	exit(0)

seqfh = open_maybe_gzip(sys.argv[1])
seqh = SeqIO.parse(seqfh, 'fastq', generic_dna)

dbs = []

print sys.argv
for db in sys.argv[2:]:
	print db
	fh = open(db, "a+")
	dbs.append(fh) #db_seqh.next())

main(seqh,dbs)
#cProfile.run('main(seqh, dbs)', 'bzz')
#p = pstats.Stats('bzz')
#p.sort_stats('cumtime')
#p.print_stats()

for db in dbs:
	db.close()

