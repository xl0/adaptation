#!/usr/bin/env python2.7

import os
import sys
from Bio import SeqIO
from collections import defaultdict
import matplotlib
matplotlib.use("Agg")

from matplotlib import pyplot

import argparse
import json

from utils import *

def main():

	parser = argparse.ArgumentParser(description='Analyze spacer frequency')
	parser.add_argument('input', metavar='<file.json>', type=str,
				help='Input file with PAMs and spacers annotated.')
	parser.add_argument('-o', metavar='out_file.json', required=True)
	args = parser.parse_args()

	(direcotry, experiment, template, tag) = split_filename(args.input)

	print 'Analyzing', args.input

	fh = open(args.input)
	seqh = SeqIO.parse(fh, "genbank")
	seq = seqh.next()

	pam_dict = {}
	non_pam_dict = defaultdict(int)

	stats = { 'num_spacers' : 0,
		  'num_spacers_pos' : 0,
		  'num_spacers_neg' : 0,
		  'hit_spacers' : 0,
		  'hit_spacers_pos' : 0,
		  'hit_spacers_neg' : 0,
		  'hit_spacers_good_pam' : 0,
		  'hit_spacers_good_pam_pos' : 0,
		  'hit_spacers_good_pam_neg' : 0,
		  'hit_spacers_bad_pam' : 0,
		  'hit_spacers_bad_pam_pos' : 0,
		  'hit_spacers_bad_pam_neg' : 0,
		  'pams' : 0,
		  'pams_hits' : 0,
		  'pams_pos' : 0,
		  'pams_pos_hits' : 0,
		  'pams_neg' : 0,
		  'pams_neg_hits' : 0
		}

	for feature in seq.features:
		if feature.type == "PAM":
			stats['pams'] += 1
			pam_dict[(int(feature.location.start), feature.location.strand)] = 0
			if feature.location.strand == 1:
				stats['pams_pos'] += 1
			else:
				stats['pams_neg'] += 1

	for feature in seq.features:
		if feature.type == "Spacer":
			hits = int(feature.qualifiers['hits'][0])
			stats['num_spacers'] += 1
			stats['hit_spacers'] += hits
			if feature.location.strand == 1:
				stats['hit_spacers_pos'] += hits
				stats['num_spacers_pos'] += 1
				if pam_dict.has_key((int(feature.location.end), 1)):
					stats['hit_spacers_good_pam'] += hits
					stats['hit_spacers_good_pam_pos'] += hits
					pam_dict[(int(feature.location.end), 1)] += hits
				else:
					stats['hit_spacers_bad_pam'] += hits
					stats['hit_spacers_bad_pam_pos'] += hits
					non_pam_dict[int(feature.location.end), 1] += hits
			else:
				stats['hit_spacers_neg'] += hits
				stats['num_spacers_neg'] += 1
				if pam_dict.has_key((int(feature.location.start) - 3, -1)):
					stats['hit_spacers_good_pam'] += hits
					stats['hit_spacers_good_pam_neg'] += hits
					pam_dict[(int(feature.location.start) - 3, -1)] += hits
				else:
					stats['hit_spacers_bad_pam'] += hits
					stats['hit_spacers_bad_pam_neg'] += hits
					non_pam_dict[int(feature.location.start) - 3, -1] += hits 




	pam_dump_dict = {}
	for key in pam_dict.keys():
		assert(not pam_dump_dict.has_key(key[0] * key[1]))
		if pam_dict[key] > 0:
			if key[1] == 1:
				stats['pams_pos_hits'] += 1
			else:
				stats['pams_neg_hits'] += 1
		
		pam_dump_dict[key[0] * key[1]] = pam_dict[key]

	stats['pams_hits'] = stats['pams_pos_hits'] + stats['pams_neg_hits']

	non_pam_dump_dict = {}
	for key in non_pam_dict.keys():
		assert(not non_pam_dump_dict.has_key(key[0] * key[1]))
		non_pam_dump_dict[key[0] * key[1]] = non_pam_dict[key]

	print 'Sequence contains %d PAMs. %d (%.2f %%) in the positive, and %d (%.2f %%) on the negative strand.' % (
		stats['pams'],
		stats['pams_pos'], float(stats['pams_pos']) / stats['pams'] * 100,
		stats['pams_neg'], float(stats['pams_neg']) / stats['pams'] * 100)
	print 'Analyzed %d spacers. %d (%.2f %%) on positive strand, %d (%.2f %%) on negative.' % (
		stats['hit_spacers'],
		stats['hit_spacers_pos'], float(stats['hit_spacers_pos']) / stats['hit_spacers'] * 100,
		stats['hit_spacers_neg'], float(stats['hit_spacers_neg']) / stats['hit_spacers'] * 100)
	print '%d (%.2f %%) spacers had a good PAM, %d (%.2f %%) on positive, %d (%.2f %%) on negative.' % (
		stats['hit_spacers_good_pam'], float(stats['hit_spacers_good_pam']) / stats['hit_spacers'] * 100,
		stats['hit_spacers_good_pam_pos'], float(stats['hit_spacers_good_pam_pos']) / stats['hit_spacers_pos'] * 100,
		stats['hit_spacers_good_pam_neg'], float(stats['hit_spacers_good_pam_neg']) / stats['hit_spacers_neg'] * 100)
	print '%d (%.2f %%) spacers didn\'t have a good pam, %d (%.2f %%) on positive, %d (%.2f) on negative.' % (
		stats['hit_spacers_bad_pam'], float(stats['hit_spacers_bad_pam']) / stats['hit_spacers'] * 100,
		stats['hit_spacers_bad_pam_pos'], float(stats['hit_spacers_bad_pam_pos']) / stats['hit_spacers_pos'] * 100,
		stats['hit_spacers_bad_pam_neg'], float(stats['hit_spacers_bad_pam_neg']) / stats['hit_spacers_neg'] * 100)

	print "%d (%.2f %%) ouf of %d pams were hit. %d (%.2f %%) of %d on positive strand, %d (%.2f %%) of %d on the negative." % (
		stats['pams_hits'], float(stats['pams_hits']) / stats['pams'] * 100, stats['pams'],
		stats['pams_pos_hits'], float(stats['pams_pos_hits']) / stats['pams_pos'] * 100, stats['pams_pos'],
		stats['pams_neg_hits'], float(stats['pams_neg_hits']) / stats['pams_neg'] * 100, stats['pams_neg'])
	

	data_file = open(args.o, 'wr+')
	data = {
		'stats' : stats,
		'pams' : pam_dump_dict,
		'non_pams' : non_pam_dump_dict,
		'experiment' : experiment,
		'template' : template,
		'template_seq' : str(seq.seq),
		'template_len' : len(seq.seq),
		'tag' : tag
	}

	json.dump(data, data_file, sort_keys=True, indent=4)
	data_file.close()


if __name__ == "__main__":
    main()
