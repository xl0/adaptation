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

	n = 0
	nseq = 0
	for db_fd in dbs:
		print "Parsing ", db_fd.name
		db_seqh = SeqIO.parse(db_fd, "genbank")
		db_seq = db_seqh.next()
		db_seqh.close()

		if db_fd.name in fast_db.keys():
			print "Duplicate DB: ", db_name
			return

		fast_db[db_fd.name] = [db_seq, db_fd, defaultdict(int)]

	print "Aligning spacer sequences..."

	# Align spacers to the DB sequence, counting the number of occurance
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

				db_spacers[(location + len(seq), 1)] += 1
				found = 1
				n += 1
				location += 1
			location = 0
			while True:
				location = db_seq.seq.find(str(seq.seq.reverse_complement()), location)
				if location == -1:
					break

				db_spacers[(location, -1)] += 1
				found = 1
				location += 1
				n += 1

			# Don't look past the first hit.
			if found:
				break

		# Could not find the whole seq - try looking at the last 20 nt.
#		if not found:



#			# Don't align on large sequences - it's too slow.
#			if len(db_sequence <= 50000):

	print "Found %d spacers in %d sequences (%d %%)" % (n, nseq, n * 100 / nseq)

	print "Updating the database files..."

	# Add spacer annotations to the db sequence
	for db_entry in fast_db.values():
		db_seq = db_entry[0]
		db_fd = db_entry[1]
		db_spacers = db_entry[2]

		for spacer_key in db_spacers.keys():
#			print spacer_key[0]
			if spacer_key[1] == 1:
				db_seq.features.append(SeqFeature(FeatureLocation(spacer_key[0] - 30, spacer_key[0]),
					type="Spacer", strand=1, qualifiers={"hits" : str(db_spacers[spacer_key])}))
			else:
				db_seq.features.append(SeqFeature(FeatureLocation(spacer_key[0], spacer_key[0] + 30),
					type="Spacer", strand=-1, qualifiers={"hits" : str(db_spacers[spacer_key])}))
		db_name = db_fd.name
		db_fd.close()

		db_fd = open(db_name, "a+")
		db_fd.truncate()
		SeqIO.write(db_seq, db_fd, "genbank")

#		print db_spacers
#		for spacer in db_spacers:


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

