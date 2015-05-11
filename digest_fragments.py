#!/usr/bin/env python
from Bio import Restriction, SeqRecord
from Bio.Seq import Seq
from StringIO import StringIO
from write import write_seq, check_format
from Bio import SeqIO
from Bio import Restriction

genbank_local_dir = "/home/chymera/genbank/"
seqence_subdir = "genomes/Eukaryotes/vertebrates_mammals/Mus_musculus/GRCm38.p3/"
assembly_subdir = "Primary_Assembly/assembled_chromosomes/FASTA/"

def select_digest(sequence, rb, max_length=10000, min_length=100):

	if isinstance(sequence, str):
		sequence = Seq(sequence)
	elif isinstance(sequence, SeqRecord.SeqRecord):
		print("MUEEE")
		sequence = sequence.seq

	digest = getattr(Restriction, rb).catalyze(sequence)

	selected_fragments = []
	for fragment in digest:
		if min_length <= len(fragment) <= max_length:
			selected_fragments += [fragment]

	return selected_fragments

if __name__ == '__main__':
	template = SeqIO.read("/home/chymera/genbank/genomes/Eukaryotes/vertebrates_mammals/Mus_musculus/GRCm38.p3/Primary_Assembly/assembled_chromosomes/FASTA/chr1.fa", "fasta")
	results = select_digest(template, "XceI", min_length=200, max_length=5000)
	print[str(i) for i in results]
