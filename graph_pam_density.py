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

	parser = argparse.ArgumentParser(description='Plot pam density')
	parser.add_argument('input', metavar='<density.json>', type=str)
	parser.add_argument('-o', metavar='out_file.json', required=True)
	args = parser.parse_args()

	print 'Visualizing PAM density %s -> %s' % (args.input, args.o)

	data = json.load(open(args.input))

	density = data['density']

	pos_density = data['pos_density']
	neg_density = data['neg_density']


	position_list = [ int(i) for i in density.keys() ]

	position_list = sorted(position_list)

	density_values = [ density[str(i)] for i in position_list]
	pos_density_values = [ pos_density[str(i)] for i in position_list]
	neg_density_values = [ neg_density[str(i)] for i in position_list]

	pyplot.plot(position_list, density_values, 'r', label="Both")
	pyplot.plot(position_list, pos_density_values, 'g', label='Positive')
	pyplot.plot(position_list, neg_density_values, 'b', label='Negative')
	pyplot.legend(loc='upper left')
	pyplot.yscale('linear')
	pyplot.title('PAM density')
#	pyplot.axis([0, len(template), 1, 11])
	pyplot.tick_params(axis='both', which='major')
	pyplot.subplots_adjust(left=0.1, right=0.9, top=0.9, bottom=0.1)
	pyplot.tight_layout()
	pyplot.savefig(args.o)

if __name__ == '__main__':
	main()

