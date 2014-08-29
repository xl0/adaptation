#!/usr/bin/env python2.7

import os
import sys
from Bio import SeqIO
from collections import defaultdict, OrderedDict
import matplotlib
matplotlib.use("Agg")

from matplotlib import pyplot

import argparse
import json

from utils import *

def main():
	parser = argparse.ArgumentParser(description='Extract spacers that consistently hit top or bottom 10%')
	parser.add_argument('inputs', metavar='<frequency.json>', type=str, nargs='+',
				help='Input file with PAMs and spacers annotated.')
	parser.add_argument('-o', metavar='out_file.json', required=True)
	args = parser.parse_args()

	print 'Extracting consistent top/bottom 20%% %s -> %s' % ( ' '.join(args.inputs), args.o)


	data_dict = OrderedDict()

	top_pam_dict = defaultdict(list)
	bottom_pam_dict = defaultdict(list)

	for infile in args.inputs:
		data = json.load(open(infile))

		experiment = data['experiment']
		tag = data['tag']
		hits = data['stats']['hit_spacers_good_pam']

		pams = data['pams']
		# Sorted bottom to top
		pam_list = sorted(pams.keys(), key = lambda e: pams[e])

		# Top/bottom 20%
		bottom_pams = pam_list[0:len(pam_list) / 8]
		top_pams = pam_list[(len(pam_list) * 8) / 10:]

		for pam in bottom_pams:
			bottom_pam_dict[pam].append(float(pams[pam]) / hits)

		for pam in top_pams:
			top_pam_dict[pam].append(float(pams[pam]) / hits)

		template = data['template_seq']

	num_measurements = len(args.inputs)

	all_top_pams = {}
	all_bottom_pams = {}

	for pam, numbers in top_pam_dict.iteritems():
		all_top_pams[pam] = (len(numbers), float(sum(numbers)) / len(numbers))

	for pam, numbers in bottom_pam_dict.iteritems():
		all_bottom_pams[pam] = (len(numbers), float(sum(numbers)) / len(numbers))

	print 'Top 10:'
	print json.dumps(all_top_pams, sort_keys=True, indent=4)
	print 'Bottom 10:'
	print json.dumps(all_bottom_pams, sort_keys=True, indent=4)

	output = {
		'top_pams' : all_top_pams,
		'bottom_pams' : all_bottom_pams,
		'template'  : template,
		'num_measurements' : num_measurements
	}

	json.dump(output, open(args.o, 'wr+'), sort_keys=True, indent=4)


if __name__ == '__main__':
	main()
