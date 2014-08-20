
from Bio.Seq import Seq, reverse_complement
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import generic_dna

class SimpleSeq:
	def __init__(self, seq, id):
		self.seq = seq
		self.id = id
		self.reverse_complement = reverse_complement(seq)

def get_primers():
	return [
		SimpleSeq("CAAAATTTTTTAGACAAAAATAGTC", "D"),
		SimpleSeq("AAGAAGAAATCAACCAGCGC", "R"),
		SimpleSeq("TAACCCTCTTTCTCAAGTTATC", "R2")
	]

["CAG", "TCA", "GTC", "AGT"]

def get_tags():
	return [
		SimpleSeq("CAG", "CAG"),
		SimpleSeq("TCA", "TCA"),
		SimpleSeq("GTC", "GTC"),
		SimpleSeq("AGT", "AGT"),
		SimpleSeq("GAC", "GAC"),
		SimpleSeq("ACT", "ACT"),
		SimpleSeq("CTG", "CTG"),
		SimpleSeq("TGA", "TGA")
	]


def get_repeat():
	return SimpleSeq("GTTTTAGAGCTATGCTGTTTTGAATGGTCCCAAAAC", "Repeat")

