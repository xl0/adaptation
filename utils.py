
import re
import gzip
import subprocess
import os

from StringIO import StringIO
from Bio import AlignIO
from Bio import SeqIO
from Bio import Seq
from Bio import SeqRecord
from Bio import SeqFeature
from Bio.SeqFeature import FeatureLocation

devnull = open(os.devnull, "w")

def open_maybe_gzip(filename):
	fh = gzip.open(filename)
	try:
		tmp = fh.readline()
	except IOError:
		fh.close()
		fh = open(filename)

		return fh

	fh.close();
	fh = gzip.open(filename)

	return fh

# Filename is '/dir/dir/experiment/4template-PAMS-TAG.gb'. Split into '/dir/dir/dir/', 'experiment', 'template', 'TAG'.
def split_filename(filename):
	(d, f) = os.path.split(filename)

	experiment = os.path.basename(d)

	name = f.split('.')[0]

	name_parts = name.split('-')
	template = name_parts[0]
	# Cut out the number used to order templates.
	template = template[1:]
	tag = name_parts[2]

	return (d, experiment, template, tag)

def sort_features(seq):
	seq.features = sorted(seq.features, key = lambda feature: int(feature.location.start))

def show_feature_type(seq, type,):
	feature_idx = 0
	cursor = 0
	string = ""

	string += type[0:3] + ": "
	for feature in seq.features:
		if feature.type == type:
			string += " " * (feature.location.start - cursor)
			string += str(seq.seq[feature.location.start:feature.location.end])
			cursor = feature.location.end
	if cursor > 0:
		print string

def show_features(seq):
	for feature in seq.features:
		print "\t", str(feature.location), ":", feature.type, feature.ref, seq.seq[feature.location.start:feature.location.end]

def show_sequence(seq):
	sort_features(seq)
	show_feature_type(seq, "Tag")
	show_feature_type(seq, "Primer")
	print "Seq: " + str(seq.seq)
	show_feature_type(seq, "Spacer")
	show_feature_type(seq, "Repeat")
	show_features(seq)


# Calculate Hamming distance. Source: http://en.wikipedia.org/wiki/Hamming_distance
def hamming_distance(s1, s2):
	if len(s1) != len(s2):
		raise ValueError("Undefined for sequences of unequal length")
	return sum(ch1 != ch2 for ch1, ch2 in zip(s1, s2))


def seq_mismatches(subseq, seq, pos = 0):
	return hamming_distance(str(subseq), str(seq[pos:pos + len(subseq)]))

def annotate_primers(seq, primers):
	for primer in primers:
		# Look for the primer on plus strand.
		positions = []
		start = 0

		# Make a list of possible start positions - matching the first 4 letters
		# XXX Rewrite in python-y way.
		while True:
			start = seq.seq.find(primer.seq[0:4], start)
			if start == -1:
				break
			positions.append(start)
			start = start + 1

		# Now actually compare the whole primer with the candidate
		for pos in positions:
			if pos + len(primer.seq) - 1 < len(seq) and seq_mismatches(primer.seq, seq.seq, pos) < 3:
				seq.features.append(SeqFeature(
					FeatureLocation(pos, pos + len(primer.seq)),
					type = "Primer", ref = primer.id, strand = +1))

		# Look for the primer on the minus strand.
		positions = []
		start = 0
		while True:
			start = seq.seq.find(primer.reverse_complement[0:4], start)
			if start == -1:
				break
			positions.append(start)
			start = start + 1

		for pos in positions:
			if pos + len(primer.seq) - 1 < len(seq) and seq_mismatches(primer.reverse_complement, seq.seq, pos) < 3:
				seq.features.append(SeqFeature(
					FeatureLocation(pos, pos + len(primer.seq)),
					type = "Primer", ref = primer.id, strand = -1))

def annotate_tag(seq, tags, start, end, strand):
	if strand == 1:
		end = start
		start = start - 3;
		if start < 0:
			return None

		seq_str = str(seq.seq[start:end])
	else:
		start = end
		end = end + 3
		if end > len(seq.seq):
			return None

		seq_str = str(seq.seq[start:end].reverse_complement())

	for tag in tags:
		if tag.seq == seq_str:
			seq.features.append(SeqFeature(
				FeatureLocation(start, end), type="Tag",
				ref=tag.id, strand = strand))
			return tag.id
	return None

def annotate_tags(seq, tags):
	found_tags = []
	# We need primers annotated before this
	for feature in seq.features:
		if feature.type == "Primer":
			new_tag = None
			new_tag = annotate_tag(seq, tags,
					int(feature.location.start),
					int(feature.location.end),
					feature.strand)
			if new_tag:
				found_tags.append(new_tag)
	return found_tags

def annotate_repeats(seq, repeat):
	primer_end = -1
	for feature in seq.features:
		if feature.type == "Primer" and feature.strand == 1:
			primer_end = feature.location.end;

	if primer_end == -1:
		return


	# Repeats are polyndromic on the ends - no need to chech RC.
	positions = []
	start = 0
	while True:
		start = seq.seq.find(repeat.seq[0:4], start)
		if start == -1:
			break
		positions.append(start)
		start = start + 1

	for position in positions:
		if position + len(repeat.seq) - 1 < len(seq) and seq_mismatches(repeat.seq, seq.seq, position) < 4:
			seq.features.append(SeqFeature(
				FeatureLocation(position, position + len(repeat.seq)),
				type = "Repeat", id = " ", strand = 1))

		if position + len(repeat.seq) - 1 < len(seq) and seq_mismatches(repeat.reverse_complement, seq.seq, position) < 4:
			seq.features.append(SeqFeature(
				FeatureLocation(position, position + len(repeat.seq)),
				type = "Repeat", id = " ", strand = -1))

def annotate_spacers(seq):

	repeats = [feature for feature in seq.features if feature.type == "Repeat"]

	# Actually should be sorted already, but sorting a sorted list is cheap, and just in case.
	repeats = sorted(repeats, key = lambda repeat:repeat.location.start)

	while len(repeats) > 1:
		seq.features.append(SeqFeature(
			FeatureLocation(repeats[0].location.end, repeats[1].location.start),
			type = "Spacer", id = " ", strand = repeats[0].strand));
		repeats = repeats[1:]

def sum_stats(stats, new):
	for ket in new.keys():
		stats[key] += new[key]

def extract_spacers(seq, tag):

	primer = tag = ""

	for feature in seq.features:
		if feature.type == "Primer" and feature.strand == 1:
			primer = feature.ref

	if not primer:
		return []

	spacers = []

	for feature in seq.features:
		if feature.type == "Spacer":
			spacer_seq = seq.seq[feature.location.start:feature.location.end]
			spacer_q = seq.letter_annotations["phred_quality"][feature.location.start:feature.location.end]
			if primer[0] == "R":
				spacer_seq = spacer_seq.reverse_complement()
				spacer_q.reverse()

			spacers.append(SeqRecord(spacer_seq, id = seq.id, name = tag, description=tag, letter_annotations = {'phred_quality' : spacer_q}))
	return spacers


def parse_water(stdout):

#	print stdout
#	fh = open(filename)
#	lines = fh.readlines()
#	fh.close()

	lines = stdout.splitlines()

	score_str = lines[27]
	score = int(float((score_str.split()[2])))

#	i = 0
#	m = re.search('[1-9]+',lines[27])
#	score = int(m.group(0))

	#parse the position
	query_position = map(lambda x: int(x),re.findall('[0-9]+',lines[32]))
	sequence_position = map(lambda x: int(x),re.findall('[0-9]+',lines[34]))

	sequence_position[0] -= 1
	query_position[0] -= 1

	return score, query_position, sequence_position


def align(spacer, refseq):


#	print len(spacer), str(spacer.seq)

	args = [
			"water",
			str("-aseq=" + refseq),
			"asis:" + spacer,
			"-gapopen=10",
			"-gapextend=10",
			"-datafile=water_score_matrix.txt",
			"stdout"
	]

#	print " ".join(args)

#	subprocess.call(args, stderr = devnull)

	child = subprocess.Popen(args, stdin = None, stdout = subprocess.PIPE, stderr = subprocess.PIPE)
#	stdout, stderr = child.communicate(str(refseq))
	stdout, stderr = child.communicate(None)

	

#	print stdout
	if not child.returncode == 0:
		print " ".join(args)
		print "stderr"
		print stderr
		raise subprocess.CalledProcessError

#	score = 0
#	pos_query = None
#	position = None
	score, pos_query, position = parse_water(stdout)
#	print score, pos_query, position

	return score, pos_query, position


def ambiguous_merge(seq1, seq2):
	if not len(seq1) == len(seq2):
		raise ValueError('Sequences should have same length')

	s_array = list(str(seq1.seq))
	q_array = seq1.letter_annotations['phred_quality']

	n = 0
	while n < len(seq1):
		if (seq1.seq[n] == seq2.seq[n]):
			# letters match - just pick the highest score here.
			q_array = max(seq1.letter_annotations["phred_quality"],
				    seq2.letter_annotations["phred_quality"])
		else:
			# No match - 'N' and average score.
			s_array[n] = "N"
			q_array[n] = seq1.letter_annotations["phred_quality"][n] + seq2.letter_annotations["phred_quality"][n] / float(2)
		n += 1

	seq_str = "".join(s_array)

	seq = SeqRecord(seq_str, id = seq1.id + seq2.id, name = seq1.name, description=seq1.description,
				letter_annotations = {'phred_quality' : q_array})

	return seq
