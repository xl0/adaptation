
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation

def show_features(seq):
	for feature in seq.features:
		print "\t", str(feature.location), ":", feature.type, feature.ref, seq.seq[feature.location.start:feature.location.end]


# Calculate Hamming distance. Source: http://en.wikipedia.org/wiki/Hamming_distance
def hamming_distance(s1, s2):
	if len(s1) != len(s2):
		raise ValueError("Undefined for sequences of unequal length")
	return sum(ch1 != ch2 for ch1, ch2 in zip(s1, s2))


def seq_mismatches(subseq, seq, pos):
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
					type = "Primer", ref = primer.id, strand = +1)
							)

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
					type = "Primer", ref = primer.id, strand = -1)
							)

def annotate_tag(seq, tags, start, end, strand):
	if strand == 1:
		end = start
		start = start - 3;
		if start < 0:
			return

		seq_str = str(seq.seq[start:end])
	else:
		start = end
		end = end + 3
		if end > len(seq.seq):
			return

		seq_str = str(seq.seq[start:end].reverse_complement())

	for tag in tags:
		if str(tag.seq) == seq_str:
			seq.features.append(SeqFeature(
				FeatureLocation(start, end), type="Tag",
				ref=tag.id, strand = strand))

def annotate_tags(seq, tags):
	# We need primers annotated before this
	for feature in seq.features:
		if feature.type == "Primer":
			annotate_tag(seq, tags,
				int(feature.location.start),
				int(feature.location.end),
				feature.strand)
       			
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
	

def extract_spacers(seq):

	primer = tag = ""

	for feature in seq.features:
		if feature.type == "Primer" and feature.strand == 1:
			primer = feature.ref
		if feature.type == "Tag" and feature.strand == 1:
			tag = feature.ref

	if not primer or not tag:
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
