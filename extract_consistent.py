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

	print 'Extracting consistent top/bottom 10%% %s -> %s' % ( ' '.join(args.inputs), args.o)


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

		bottom_pams_10 = pam_list[0:len(pam_list) / 10]
		top_pams_10 = pam_list[(len(pam_list) * 9) / 10:]

		for pam in bottom_pams_10:


#		data_dict[(experiment, tag)] = data

#	for key in data_dict.keys():
#		data = data_dict



if __name__ == '__main__':
	main()
