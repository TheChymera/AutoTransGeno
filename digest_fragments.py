#!/usr/bin/env python
from Bio import Restriction, SeqRecord
from Bio.Seq import Seq
from StringIO import StringIO
from write import write_seq, check_format
from Bio import SeqIO
from Bio.SeqUtils import GC
from Bio import Restriction

def select_digest(sequence, rb, max_length=10000, min_length=100, max_gc=100, min_gc=0, save_to=False):

	if isinstance(sequence, str):
		sequence = Seq(sequence)
	elif isinstance(sequence, SeqRecord.SeqRecord):
		sequence = sequence.seq

	digest = getattr(Restriction, rb).catalyze(sequence)

	selected_fragments = []
	for ix, fragment in enumerate(digest):
		if min_length <= len(fragment) <= max_length and min_gc <= GC(fragment) <= max_gc:
			selected_fragments += [fragment]
			if save_to:
				import time
				from os import makedirs, path
				append_date = "/"+time.strftime("%Y%m%d_%H%M")+"/"
				save_path = save_to + append_date
				if not path.exists(save_path):
					makedirs(save_path)
				write_seq(sequence_write_path=save_path, entrez_id="", sequence=fragment, sequence_id=str(len(fragment))+"-"+str(ix), formats=["fasta"])

	return selected_fragments

if __name__ == '__main__':
	template = SeqIO.read("/home/chymera/genbank/genomes/Eukaryotes/vertebrates_mammals/Mus_musculus/GRCm38.p3/Primary_Assembly/assembled_chromosomes/FASTA/chr1.fa", "fasta")
	results = select_digest(template, "XceI", min_length=1500, max_length=1500, max_gc=55, min_gc=30, save_to="digests/")
