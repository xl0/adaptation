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
	parser = argparse.ArgumentParser(description='Visualize spacer frequency')
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


	num_sequences = len(data_dict)
	f = pyplot.figure(figsize=(15, 5 * num_sequences))

	n = 0
	fontsize = 10
	for entry in data_dict.values():
		pams = entry['pams']

		pos_pam_list = [ int(key) for key in pams.keys() if int(key) > 0 and pams[key] > 0]
		neg_pam_list = [ -int(key) for key in pams.keys() if int(key) < 0 and pams[key] > 0]

		non_pams = entry['non_pams']

		pos_nonpam_list = [ int(key) for key in non_pams.keys() if int(key) > 0 and non_pams[key] > 0]
		neg_nonpam_list = [ -int(key) for key in non_pams.keys() if int(key) < 0 and non_pams[key] > 0]

		num_spacer_hits = entry['stats']['hit_spacers']

		pos_pam_list = sorted(pos_pam_list)
		neg_pam_list = sorted(neg_pam_list)

		pos_nonpam_list = sorted(pos_nonpam_list)
		neg_nonpam_list = sorted(neg_nonpam_list)

		pos_pam_score_list = []
		for key in pos_pam_list:
			pos_pam_score_list.append(float(pams[str(key)]) / num_spacer_hits)

		neg_pam_score_list = []
		for key in neg_pam_list:
			neg_pam_score_list.append(float(pams[str(-key)]) / num_spacer_hits)

		pos_nonpam_score_list = []
		for key in pos_nonpam_list:
			pos_nonpam_score_list.append(float(non_pams[str(key)]) / num_spacer_hits)

		neg_nonpam_score_list = []
		for key in neg_nonpam_list:
			neg_nonpam_score_list.append(float(non_pams[str(-key)]) / num_spacer_hits)

		pyplot.subplot(num_sequences, 2, n * 2 + 1)
		pyplot.plot(pos_pam_list, pos_pam_score_list, 'bo', pos_nonpam_list, pos_nonpam_score_list, 'go', )
		pyplot.yscale('log')
		pyplot.title('%s:%s +' % (entry['experiment'], entry['tag']), fontsize=fontsize+3)
		pyplot.axis([0, entry['template_len'], 0, 1])
		pyplot.tick_params(axis='both', which='major', labelsize=fontsize)
		pyplot.subplots_adjust(left=0.1, right=0.9, top=0.9, bottom=0.1)

		pyplot.subplot(num_sequences, 2, n * 2 + 2 )
		pyplot.plot(neg_pam_list, neg_pam_score_list, 'bo', neg_nonpam_list, neg_nonpam_score_list, 'go')
		pyplot.yscale('log')
		pyplot.title('%s:%s -' % (entry['experiment'], entry['tag']), fontsize=fontsize+3)
		pyplot.axis([0, entry['template_len'], 0, 1])
		pyplot.tick_params(axis='both', which='major', labelsize=fontsize)
		pyplot.subplots_adjust(left=0.1, right=0.9, top=0.9, bottom=0.1)
		n += 1

	pyplot.tight_layout()
	pyplot.savefig(args.o)

if __name__ == '__main__':
	main()


