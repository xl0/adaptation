
from Bio import SeqIO
from Bio.SeqFeature import SeqFeature, FeatureLocation

def load_sequences(filename):
	fh = open(filename, "r")
	
	sequences = []
	for seq in SeqIO.parse(fh, "fasta"):
		sequences.append(seq)
	
	return sequences

def show_features(seq):
	for feature in seq.features:
		print "\t", str(feature.location), ":", feature.type, feature.ref, seq.seq[feature.location.start:feature.location.end]


# Calculate Hammind distance. Source: http://en.wikipedia.org/wiki/Hamming_distance
def hamming_distance(s1, s2):
	if len(s1) != len(s2):
		raise ValueError("Undefined for sequences of unequal length")
	return sum(ch1 != ch2 for ch1, ch2 in zip(s1, s2))


def seq_mismatches(subseq, seq, pos):
	return hamming_distance(str(subseq), str(seq[pos:pos + len(subseq)]))

def annotate_primers(seq, primers):

	for primer in primers:
		# Look for the primer on plus strand.
		positions = [i for i in range(len(seq)) if seq.seq.startswith(primer.seq[0:4], i)]
		for pos in positions:
			if pos + len(primer) - 1 < len(seq) and seq_mismatches(primer.seq, seq.seq, pos) < 3:
				seq.features.append(SeqFeature(
					FeatureLocation(pos, pos + len(primer)),
					type = "Primer", ref = primer.id, strand = +1)
							)
				

		# Look for the primer on the minus strand.
		positions = [i for i in range(len(seq)) if seq.seq.startswith(primer.seq.reverse_complement()[0:4], i)]
		for pos in positions:
			if pos + len(primer) - 1 < len(seq) and seq_mismatches(primer.seq.reverse_complement(), seq.seq, pos) < 3:
				seq.features.append(SeqFeature(
					FeatureLocation(pos, pos + len(primer)),
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

	positions = [i for i in range(primer_end + 1, len(seq)) if seq.seq.startswith(repeat.seq[0:4], i)]
	for position in positions:
		if position + len(repeat) - 1 < len(seq) and seq_mismatches(repeat.seq, seq.seq, position) < 4:
			seq.features.append(SeqFeature(
				FeatureLocation(position, position + len(repeat)),
				type = "Repeat", id = "", strand = 1))

		if position + len(repeat) - 1 < len(seq) and seq_mismatches(repeat.seq.reverse_complement(), seq.seq, position) < 4:
			seq.features.append(SeqFeature(
				FeatureLocation(position, position + len(repeat)),
				type = "Repeat", id = "", strand = -1))

