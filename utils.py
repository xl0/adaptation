
from Bio import SeqIO

def load_primers(filename):
	primer_fh = open(filename, "r");
	
	primers = [];
	for primer in SeqIO.parse(primer_fh, "fasta"):
		primers.append(primer);
	
	return primers;
			
	
