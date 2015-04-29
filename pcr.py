#!/usr/bin/env python

from Bio.Blast.Applications import NcbiblastnCommandline
from Bio.Seq import Seq
from Bio.Blast import NCBIXML
from StringIO import StringIO
from write import write_seq
import itertools
from Bio import SeqIO
import os

genbank_local_dir = "/home/chymera/genbank/"
seqence_subdir = "genomes/Eukaryotes/vertebrates_mammals/Mus_musculus/GRCm38.p3/"
assembly_subdir = "Primary_Assembly/assembled_chromosomes/FASTA/"

sequence_path = genbank_local_dir + seqence_subdir + assembly_subdir

def mock_pcr(template, primer_files="", primer_seqs="", primer_directory="", max_length=3000):
	"""
	Performs a rudimentary in-silico PCR, which returns the sequence and length of template fragments flanked by primers shother than max_length.
	"""
	if primer_files == "" and primer_seqs == "":
		raise Exception("Please specify at least two entries for `primer_seqs`, `primer_files`, or one for each.")
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
	if "." in template:
		template_list = [check_format(template, "fasta")]
	elif type(template) == str:
		dest = write_seq(sequence=template, sequence_id="template")
		templat_list = [dest]
	elif "/" in template:
		os.listdir(template)


	#writes sequences to files
	if primer_seqs:
		for ix, primer_seq in enumerate(primer_seqs):
			if type(primer_seq) == str:
				dest = write_seq(sequence=primer_seq, sequence_id="primer"+str(ix))
				primer_files += [dest]

	print primer_files

	for template in template_list:
		fw_primer_start_positions=[]
		rv_primer_start_positions=[]
		template_seq = SeqIO.read(template, "fasta")
		for primer_file in primer_files:
			output = NcbiblastnCommandline(query=template, subject=primer_file, outfmt=5, task="blastn")()[0]
			blast_result_record = NCBIXML.read(StringIO(output))

			for alignment in blast_result_record.alignments:
				hsp_list=alignment.hsps
				for hsp in hsp_list:
					print "muie:", hsp.sbjct_start, hsp.sbjct_end, hsp.query_start, hsp.query_end, hsp.align_length, hsp.sbjct
					if hsp.sbjct_start > hsp.sbjct_end:
						strand = -1
						amplicon_5_end = hsp.query_end
						rv_primer_start_positions += [hsp.query_end]
					else:
						strand = +1
						amplicon_5_end = hsp.query_start
						fw_primer_start_positions += [hsp.query_start]

		for fw_primer_start in fw_primer_start_positions:
			for rv_primer_start in rv_primer_start_positions:
				amplicon_length = rv_primer_start - fw_primer_start
				if amplicon_length <= max_length:
					return amplicon_length, template_seq[fw_primer_start:rv_primer_start].seq


def genome_mock_pcr(genome_path, unspecific_insert, primer_files="", primer_seqs="", primer_directory="", max_length=3000):


if __name__ == '__main__':
	template = "AAAAAAAAAAAAAAAGGGGGGCCCCCCCCAAAAAAAAAAAAGAGGTGCATACCTTTGCGCATCATCTCCCAGGCTTGGGCTTCCTTTAGGGTAACTGGCCGCCGCCATGTTGCAAACGGGAAGGAAATGAATGAACCGCCGTTATGAAATCTTGCTTAGGCCTTCCTTCTTCCTAGCTTGTGACTAACCTCATTCCTCTCGGCTGGGTGGAGTGTCCTTTATCCTGTAGGCCAGGTGATGCAAAGCACGGAAGCCGGCGAGAGTTCTACCTCACAATCTGTCTCACCTTATTAGCCTTAAAAGCCCTTGAGCCTTATTGTCCTCGGGCATAATGCGTATTCTAGATTATTCTCTGAAAATCAAAGCGGACTTACAGAGGTCCGCTTGACCTCCCAACCCCAGAGGTAGTTATGGCGTAGTGCAGAGCCGTGGGATGGGGAGCTGAGTCATGGTGGTTCTGAAAAGAAATTTTCCACCACAAAATGGCTCCTGTAGTAGCAGCCCCTTCCATCCCCTGCACTTCCCATCACAGCCTCGCACTGACCCAGGCCCTATAGGCCAGGATGTAAAGGTCATTAAGAGGATTGGGTGTCCCTGCGCCTCAGAATCCTGCCCTTCTCCCCGTTCCATCCTCCAGAAACCAGATCTCTCCACTCCGCCCTGATCTGAGGTTAAATTTAGCAGCACGGAAGCCGGCCGTGTGACCTTTCTGGATCTGGGGTCTGAGCGGGCTCTCCACCCTGCTCCCCCTACACACATCTGGGAAATCCCATTGACTTGCTCCGGCTCTCATTTTTGCCCGAGAAGAACAGGTGTTTCGCGAACGAGCCCTGGGATTAGGGTCGTGCCTTCGGCCGTTGGAAACCCCCCACATGTTTTCTCAGTCTTTCCCCTTAGTTCGAGGGACTTGGAGGACACAGGTGGGCCCGCCCTGTGCTGCTCACGCTGACCTTTAGCCTTGCCCTTTGAGCTTGCTGATGAATGAGTTCACAGGTCTGCCCTGTCCAGGGGGTGTAGCCTGAAGTCCAGCCATGCTGGAACAAACTTCCCAGGGCATGAGTGAT"
	mock_pcr(template, primer_seqs=["CCTTTAGGGTAACTG", "GCCGGCTTCCGTGCT"])
