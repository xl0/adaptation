#!/usr/bin/env python2.7

import os
import sys
from Bio import SeqIO
from Bio import motifs

import argparse
import json

def main():

	parser = argparse.ArgumentParser(description='Extract the top and bottom 10% of the spacers')
	parser.add_argument('input', metavar='<frequency.json>', type=str,
				help='Input file with PAMs and spacers annotated.')
	parser.add_argument('-o', metavar='out_base', required=True)
	args = parser.parse_args()

	print 'Visualizing %s -> %s' % (args.input, args.o)

	fd = open(args.input)
	data = json.load(fd)
	fd.close()

	before = 75
	after = 50

	template = data['template_seq']
	template_len = data['template_len']
	template_seq = data['template_seq']
	pams = data['pams']

	experiment = data['experiment']
	tag = data['tag']

	pos_pam_list = [key for key in pams.keys() if int(key) > 0 and int(key) > before and template_len - int(key) > after]
	pos_pam_list = sorted(pos_pam_list, key = lambda e: pams[e])

	pos_bottom_pams_10 = pos_pam_list[0:len(pos_pam_list) / 10]
	pos_top_pams_10 = pos_pam_list[(len(pos_pam_list) * 9) / 10:]

#	pos_mid_pams_10 = pos_pam_list[len(pos_pam_list) * 4.5 / 10:len(pos_pam_list) * 5.5 / 10 ]
#	pos_top_pams_10 = pos_pam_list[(len(pos_pam_list) * 9) / 10:]

	pos_bottom_pams_1 = pos_pam_list[0:len(pos_pam_list) / 100]
	pos_top_pams_1 = pos_pam_list[(len(pos_pam_list) * 99) / 100:]

#	print "bottom 1%"
	pos_bottom_1_seqs = []
	for pam in pos_bottom_pams_1:
		pos_bottom_1_seqs.append(template_seq[int(pam) - before:int(pam) + after])
#		print "%8d %s" % (pams[pam], template_seq[int(pam) - before:int(pam) + after])
#	print "bottom 10%"
	pos_bottom_10_seqs = []
	for pam in pos_bottom_pams_10:
#		print "%8d %s" % (pams[pam], template_seq[int(pam) - before:int(pam) + after])
		pos_bottom_10_seqs.append(template_seq[int(pam) - before:int(pam) + after])

#	print "top 1%"
	pos_top_1_seqs = []
	for pam in pos_top_pams_1:
		pos_top_1_seqs.append(template_seq[int(pam) - before:int(pam) + after])
#		print "%8d %s" % (pams[pam], template_seq[int(pam) - before:int(pam) + after])
#	print 'top 10%'
	pos_top_10_seqs = []
	for pam in pos_top_pams_10:
#		print '%8d %s' % (pams[pam], template_seq[int(pam) - before:int(pam) + after])
		pos_top_10_seqs.append(template_seq[int(pam) - before:int(pam) + after])

	pos_seqs = []
	for pam in pos_pam_list:
		pos_seqs.append(template_seq[int(pam) - before:int(pam) + after])

	pos_bottom1_motif = motifs.create(pos_bottom_1_seqs)
	pos_bottom10_motif = motifs.create(pos_bottom_10_seqs)

	pos_top1_motif = motifs.create(pos_top_1_seqs)
	pos_top10_motif = motifs.create(pos_top_10_seqs)

	pos_motif = motifs.create(pos_seqs)
#	print pos_top1_motif.counts
#	print pos_top10_motif.counts

	pos_bottom1_motif.weblogo(args.o + '_bottom1.png', logo_title='%s:%s bottom 1%%' % (experiment, tag))
	pos_bottom10_motif.weblogo(args.o + '_bottom10.png', logo_title='%s:%s bottom 10 %%' % (experiment, tag))

	pos_top1_motif.weblogo(args.o + '_top1.png', logo_title='%s:%s top 1%%' % (experiment, tag))
	pos_top10_motif.weblogo(args.o + '_top10.png', logo_title='%s:%s top 10%%' % (experiment, tag))

	pos_motif.weblogo(args.o + '_all.png', logo_title='%s:%s All sequences' % (experiment, tag))

if __name__ == '__main__':
	main()
