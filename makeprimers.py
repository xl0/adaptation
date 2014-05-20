#!/usr/bin/env python2.7

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import generic_dna

fh = open("primers.fasta", "a+")



primer1 = SeqRecord(Seq("CAAAATTTTTTAGACAAAAATAGTC", generic_dna), id="D")
primer2 = SeqRecord(Seq("AAGAAGAAATCAACCAGCGC", generic_dna), id="R")
primer3 = SeqRecord(Seq("TAACCCTCTTTCTCAAGTTATC", generic_dna), id="R2")

SeqIO.write(primer1, fh, "fasta")
SeqIO.write(primer2, fh, "fasta")
SeqIO.write(primer3, fh, "fasta")


fh.close()

