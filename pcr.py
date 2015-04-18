#!/usr/bin/env python

from Bio.Emboss.Applications import PrimerSearchCommandline
from Bio.Seq import Seq

genbank_local_dir = "/home/chymera/genbank/"
seqence_subdir = "genomes/Eukaryotes/vertebrates_mammals/Mus_musculus/GRCm38.p3/"
assembly_subdir = "Primary_Assembly/assembled_chromosomes/FASTA/"

sequence_path = genbank_local_dir + seqence_subdir + assembly_subdir

def mock_pcr(template, primer_files="", primer_seqs="", primer_directory="", max_length=3000):
	"""
	Performs a rudimentary in-silico PCR, which returns the sequence and length of template fragments flanked by primers shother than max_length.
	"""
	if primer_files == "" and primer_seqs == "":
		raise Exception("Please specify at least two entries for `primer_seqs`, `primer_seqs`, or one for each.")
	if len(primer_files)+len(primer_seqs) < 2:
		raise Exception("Please specify at least two primers.")

	#converts files from whatever format they have into FASTA (required for ncbi blast)
	if primer_files:
		for ix, primer_file in enumerate(primer_files):
			new_primer_file = check_format(subject, "fasta")
			primer_files[ix] = new_primer_file
	else:
		primer_files=[]

	#writes sequences to files
	if primer_seqs:
		from write import write_seq
		for ix, primer_seq in enumerate(primer_seqs):
			if type(primer_seq) == str:
				dest = write_seq(sequence=primer_seq, sequence_id="primer"+str(ix))
				primer_files += dest

	for 

if __name__ == '__main__':
	template = "GAGGTGCATACCTTTGCGCATCATCTCCCAGGCTTGGGCTTCCTTTAGGGTAACTGGCCGCCGCCATGTTGCAAACGGGAAGGAAATGAATGAACCGCCGTTATGAAATCTTGCTTAGGCCTTCCTTCTTCCTAGCTTGTGACTAACCTCATTCCTCTCGGCTGGGTGGAGTGTCCTTTATCCTGTAGGCCAGGTGATGCAAGGCTTCCGTGCTCTCGAGAGAGTTCTACCTCACAATCTGTCTCACCTTATTAGCCTTAAAAGCCCTTGAGCCTTATTGTCCTCGGGCATAATGCGTATTCTAGATTATTCTCTGAAAATCAAAGCGGACTTACAGAGGTCCGCTTGACCTCCCAACCCCAGAGGTAGTTATGGCGTAGTGCAGAGCCGTGGGATGGGGAGCTGAGTCATGGTGGTTCTGAAAAGAAATTTTCCACCACAAAATGGCTCCTGTAGTAGCAGCCCCTTCCATCCCCTGCACTTCCCATCACAGCCTCGCACTGACCCAGGCCCTATAGGCCAGGATGTAAAGGTCATTAAGAGGATTGGGTGTCCCTGCGCCTCAGAATCCTGCCCTTCTCCCCGTTCCATCCTCCAGAAACCAGATCTCTCCACTCCGCCCTGATCTGAGGTTAAATTTAGCCGTGTGACCTTTCTGGATCTGGGGTCTGAGCGGGCTCTCCACCCTGCTCCCCCTACACACATCTGTTGCTCCGGCTCTCATTTTTGCCCGAGAAGAACAGGTGTTTCGCGAACGAGCCCTGGGATTAGGGTTGGAAACCCCCCACATGTTTTCTCAGTCTTTCCCCTTAGTTCGAGGGACTTGGAGGACACAGGTGGGCCCGCCCTGTGCTGCTCACGCTGACCTTTAGCCTTGCCCTTTGAGCTTGCTGATGAATGAGTTCACAGGTCTGCCCTGTCCAGGGGGTGTAGCCTGAAGTCCAGCCATGCTGGAACAAACTTCCCAGGGCATGAGTGAT"
	mock_pcr(template, primer_seqs=["CCTTTAGGGTAA", "GGCTTCCGTGCT"])
