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
def sign(number):
	if number == abs(number):
		return 1
	else:
		return -1

def main():
	parser = argparse.ArgumentParser(description='PAM hit frequency between experiments')
	parser.add_argument('input', metavar='<file.gb>', type=str, nargs='+',
				help='Input file with PAMs and spacers annotated.')
	parser.add_argument('-o', metavar='text-output-file.txt', required=True)
	parser.add_argument('-j', metavar='json-output-file.json')
	args = parser.parse_args()

	ofile = open(args.o, 'wr+')
	if args.j:
		jfile = open(args.j, 'wr+')

#	(direcotry, template, tag) = split_filename(args.input)

	print 'Analyzing', args.input
	
	data = OrderedDict()
	for i in args.input:
		print 'Reading', i
		f = open(i)
		data[i] = json.load(f)
		f.close

	accum_pam_dict = {}
	entry = data.values()[0]
	for pam in entry['pams'].keys():
		accum_pam_dict[int(pam)] = True

	accum_non_pam_dict = {}
	for entry in data.values():
		for non_pam in entry['non_pams'].keys():
			accum_non_pam_dict[int(non_pam)] = True

#	print '#', len(accum_pam_dict)

	accum_pam_list = sorted(accum_pam_dict.keys(), key=lambda entry: abs(entry))
#	print json.dumps(accum_pam_dict, indent=4) 

	accum_non_pam_list = sorted(accum_non_pam_dict.keys(), key=lambda entry: abs(entry))

	ofile.write('PAM hits per tag:\n')
	for entry in data.values():
		stats = entry['stats']
		ofile.write ("%s: %d (%d on positive, %d on negative strand)\n" % (
			entry['tag'],
			stats['hit_spacers_good_pam'],
			stats['hit_spacers_good_pam_pos'],
			stats['hit_spacers_good_pam_neg']))

	ofile.write ('\nPAM hits\n')
	header = 'Position Strand\t'
	for entry in data.values():
		header += entry['experiment'] + ':' + entry['tag'] + '\t'
	header += '\n'

	for i in range(0, len(accum_pam_dict)):
		if not i % 50:
			ofile.write(header)
		line = '%8d %2d     ' % (abs(accum_pam_list[i]), sign(accum_pam_list[i]))
		for entry in data.values():
			line += '%d\t' % entry['pams'][str(accum_pam_list[i])]
		line += '\n'
		ofile.write(line)

	ofile.write('Non-PAM hits per tag:\n')
	for entry in data.values():
		stats = entry['stats']
		ofile.write("%s: %d (%d on positive, %d on negative strand)\n" % (
			entry['tag'],
			stats['hit_spacers_bad_pam'],
			stats['hit_spacers_bad_pam_pos'],
			stats['hit_spacers_bad_pam_neg']))
	
	for i in range(0, len(accum_non_pam_dict)):
		if not i % 50:
			ofile.write(header)
		line = '%8d %2d     ' % (abs(accum_non_pam_list[i]), sign(accum_non_pam_list[i]))
		for entry in data.values():
			if entry['non_pams'].has_key(str(accum_non_pam_list[i])):
				hits = entry['non_pams'][str(accum_non_pam_list[i])]
			else:
				hits = 0
			line += '%d\t' % hits
		line += '\n'
		ofile.write(line)

	if args.j:
		json.dump(data, jfile, sort_keys=True, indent=4)
		jfile.close()
	ofile.close()


if __name__ == "__main__":
    main()
