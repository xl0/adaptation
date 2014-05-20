#!/usr/bin/env python2.7

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import generic_dna

fh = open("tags.fasta", "a+")

tags = []

for seq in ["CAG", "TCA", "GTC", "AGT"]:
	SeqIO.write(SeqRecord(Seq(seq, generic_dna), id=seq), fh, "fasta")

fh.close()

