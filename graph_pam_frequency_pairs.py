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
	parser = argparse.ArgumentParser(description='Visualize adapration correlations')
	parser.add_argument('inputs', metavar='<frequency.json>', type=str, nargs='+',
				help='Input file with PAMs and spacers annotated.')
	parser.add_argument('-o', metavar='out_file.svg/png', required=True)
	args = parser.parse_args()

	print 'Visualizing %s -> %s' % ( ' '.join(args.inputs), args.o)

	data_dict = OrderedDict()

	for infile in args.inputs:
		fd = open(infile)
		data = json.load(fd)
		fd.close()

		experiment = data['experiment']
		tag = data['tag']

		data_dict[(experiment, tag)] = data

	entry = data_dict[data_dict.keys()[0]]
	
	pam_list = entry['pams'].keys()

	pam_hits_dict = {}
	for key in data_dict.keys():
		print key
		hits_list = []
		for pam_key in data_dict[key]['pams'].keys():
			hits_list.append(float(data_dict[key]['pams'][pam_key]) / data_dict[key]['stats']['hit_spacers_good_pam'])
		
		pam_hits_dict[key] = hits_list

	num_sequences = len(pam_hits_dict)

	pyplot.figure(figsize=(5 * num_sequences, 5 * num_sequences))	


	i = 0
	for key1 in data_dict.keys():
		j = 1
		for key2 in data_dict.keys():
			pyplot.subplot(num_sequences, num_sequences, i * num_sequences + j)
	
			pyplot.plot(pam_hits_dict[key2], pam_hits_dict[key1], 'bo')

			pyplot.yscale('log')
			pyplot.ylabel('%s:%s' % (data_dict[key1]['experiment'], data_dict[key1]['tag']))
			pyplot.xlabel('%s:%s' % (data_dict[key2]['experiment'], data_dict[key2]['tag']))
			pyplot.xscale('log')
			pyplot.title('%s:%s - %s:%s' % (data_dict[key1]['experiment'], data_dict[key1]['tag'],
								data_dict[key2]['experiment'], data_dict[key2]['tag']), fontsize=11)
			pyplot.axis([-1, 1, -1, 1])
			pyplot.tick_params(axis='both', which='major', labelsize='small')
			pyplot.subplots_adjust(left=0.1, right=0.9, top=0.9, bottom=0.1)
			j += 1
		i += 1

	pyplot.tight_layout()
	pyplot.savefig(args.o)
		



if __name__ == '__main__':
	main()
