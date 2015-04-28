#!/usr/bin/env python

from Bio.Blast.Applications import NcbiblastnCommandline
from Bio.Seq import Seq
from Bio.Blast import NCBIXML
from StringIO import StringIO
from write import write_seq

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

	#convert template format if necessary, if template not formatte a path write chached file for seuence
	if "." in template or "/" in template:
		template = check_format(template, "fasta")
	elif type(template) == str:
		dest = write_seq(sequence=template, sequence_id="template")
		template = dest

	#writes sequences to files
	if primer_seqs:
		for ix, primer_seq in enumerate(primer_seqs):
			if type(primer_seq) == str:
				dest = write_seq(sequence=primer_seq, sequence_id="primer"+str(ix))
				primer_files += [dest]

	print primer_files

	primer_start_positions=[]
	reverse_complement_start_positions=[]
	for primer_file in primer_files:
		output = NcbiblastnCommandline(query=template, subject=primer_file, outfmt=5, task="blastn")()[0]
		blast_result_record = NCBIXML.read(StringIO(output))

		for alignment in blast_result_record.alignments:
			hsp_list=alignment.hsps
			for hsp in hsp_list:
				print "muie:", hsp.sbjct_start, hsp.sbjct_end, hsp.query_start, hsp.query_end
				if hsp.sbjct_start > hsp.sbjct_end:
					strand = -1
					amplicon_5_end = hsp.query_end
				else:
					strand = +1
					amplicon_5_end = hsp.query_start
				print strand

if __name__ == '__main__':
	template = "AAAAAAAAAAAAAAAGGGGGGCCCCCCCCAAAAAAAAAAAAGAGGTGCATACCTTTGCGCATCATCTCCCAGGCTTGGGCTTCCTTTAGGGTAACTGGCCGCCGCCATGTTGCAAACGGGAAGGAAATGAATGAACCGCCGTTATGAAATCTTGCTTAGGCCTTCCTTCTTCCTAGCTTGTGACTAACCTCATTCCTCTCGGCTGGGTGGAGTGTCCTTTATCCTGTAGGCCAGGTGATGCAAAGCACGGAAGCCGGCTTCCGTGCTCTCGAGAGAGTTCTACCTCACAATCTGTCTCACCTTATTAGCCTTAAAAGCCCTTGAGCCTTATTGTCCTCGGGCATAATGCGTATTCTAGATTATTCTCTGAAAATCAAAGCGGACTTACAGAGGTCCGCTTGACCTCCCAACCCCAGAGGTAGTTATGGCGTAGTGCAGAGCCGTGGGATGGGGAGCTGAGTCATGGTGGTTCTGAAAAGAAATTTTCCACCACAAAATGGCTCCTGTAGTAGCAGCCCCTTCCATCCCCTGCACTTCCCATCACAGCCTCGCACTGACCCAGGCCCTATAGGCCAGGATGTAAAGGTCATTAAGAGGATTGGGTGTCCCTGCGCCTCAGAATCCTGCCCTTCTCCCCGTTCCATCCTCCAGAAACCAGATCTCTCCACTCCGCCCTGATCTGAGGTTAAATTTAGCCGTGTGACCTTTCTGGATCTGGGGTCTGAGCGGGCTCTCCACCCTGCTCCCCCTACACACATCTGTTGCTCCGGCTCTCATTTTTGCCCGAGAAGAACAGGTGTTTCGCGAACGAGCCCTGGGATTAGGGTTGGAAACCCCCCACATGTTTTCTCAGTCTTTCCCCTTAGTTCGAGGGACTTGGAGGACACAGGTGGGCCCGCCCTGTGCTGCTCACGCTGACCTTTAGCCTTGCCCTTTGAGCTTGCTGATGAATGAGTTCACAGGTCTGCCCTGTCCAGGGGGTGTAGCCTGAAGTCCAGCCATGCTGGAACAAACTTCCCAGGGCATGAGTGAT"
	mock_pcr(template, primer_seqs=["CCTTTAGGGTAA", "GGCTTCCGTGCT"])
