#!/usr/bin/env python
from Bio import Restriction, SeqRecord
from Bio.Seq import Seq
from StringIO import StringIO
from write import write_seq, check_format
from Bio import SeqIO
from Bio import Restriction

def select_digest(sequence, rb, max_length=10000, min_length=100, save_to=False):

	if isinstance(sequence, str):
		sequence = Seq(sequence)
	elif isinstance(sequence, SeqRecord.SeqRecord):
		sequence = sequence.seq

	digest = getattr(Restriction, rb).catalyze(sequence)

	selected_fragments = []
	for ix, fragment in enumerate(digest):
		if min_length <= len(fragment) <= max_length:
			selected_fragments += [fragment]
			if save_to:
				write_seq(sequence_write_path=save_to, entrez_id="", sequence=fragment, sequence_id=str(len(fragment))+"-"+str(ix), formats=["fasta"])

	return selected_fragments

if __name__ == '__main__':
	template = SeqIO.read("/home/chymera/genbank/genomes/Eukaryotes/vertebrates_mammals/Mus_musculus/GRCm38.p3/Primary_Assembly/assembled_chromosomes/FASTA/chr1.fa", "fasta")
	results = select_digest(template, "XceI", min_length=800, max_length=2500, save_to="digests/")
