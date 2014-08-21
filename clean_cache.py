#!/usr/bin/env python2.7

import json
import argparse

def main():
	parser = argparse.ArgumentParser(description='Remove non-aligned spacers from alignment cache dictionary')
	parser.add_argument('dictionary', metavar='<dictionary.json>', type=str,
				help='Alignment cache dictionary')

	args = parser.parse_args()
	print 'cleaning', args.dictionary

	fh = open(args.dictionary)

	cache = json.load(fh)
	fh.close()

	clean_cache = {}

	for key in cache.keys():
		if cache[key]:
			clean_cache[key] = cache[key]

	fh = open(args.dictionary, 'wr+')
	json.dump(clean_cache, fh)
	fh.close()		

if __name__ == "__main__":
    main()
