#!/usr/bin/env python2.7

import os
import sys
from Bio import SeqIO
from Bio import Seq
from collections import defaultdict, OrderedDict
import itertools
import matplotlib
#matplotlib.use("Agg")

from matplotlib import pyplot

import argparse
import json

from utils import *

def get_last_spacer_bits(seq, pam, nbits):
	pam = int(pam)
	if abs(pam) < nbits:
		return None

	if pam > 0:
		return seq[pam - nbits:pam]
	else:
		return Seq.reverse_complement(seq[-pam + 3: -pam + 3 + nbits:])

def main():
	parser = argparse.ArgumentParser(description='Visualize spacer spacer distribution')
	parser.add_argument('input', metavar='<frequency.json>', type=str,
				help='Input file with PAMs and spacers annotated.')
	parser.add_argument('-o', metavar='out_file.svg/png', required=True)
	args = parser.parse_args()

	print 'Spacer distributions %s -> %s' % (args.input, args.o)

	fd = open(args.input)
	data = json.load(fd)
	fd.close()

	template = data['template_seq']
	template_len = data['template_len']
	template_seq = str(data['template_seq'])
	pams = data['pams']

	experiment = data['experiment']
	tag = data['tag']

	product_list = list(''.join(p) for p in itertools.product('ATCG', repeat=4))

	spacer_part_dict = {}
	spacer_hit_dict = {}
	for i in product_list:
		spacer_part_dict[i] = 0
		spacer_hit_dict[i] = 0

	for pam in pams:
		seq = get_last_spacer_bits(template_seq, pam, 4)
		spacer_part_dict[seq] += 1
		spacer_hit_dict[seq] += pams[pam]

	xlist = range(0, len(product_list))
	ylist = [float(spacer_hit_dict[key]) / spacer_part_dict[key] if spacer_part_dict[key] > 0 else 0 for key in product_list]

	pyplot.bar(xlist, ylist, color='b')
	pyplot.xlim([0, len(product_list)])
	#pyplot.yscale('log')
#	pyplot.show()

	pyplot.savefig(args.o)

	sorted_products = sorted(product_list, key=lambda key: float(spacer_hit_dict[key]) / spacer_part_dict[key] if spacer_part_dict[key] > 0 else 0, reverse=True)

	print 'Zeros:'
	for i in product_list:
		if spacer_part_dict[i] > 0 and spacer_hit_dict[i] == 0:
			print i

	print 'Top 10:'
	for i in range(0, 10):
		if spacer_part_dict[sorted_products[i]] > 0:
			print sorted_products[i], float(spacer_hit_dict[sorted_products[i]]) / spacer_part_dict[sorted_products[i]]

if __name__ == '__main__':
	main()
