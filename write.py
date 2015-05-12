#!/usr/bin/env python
__author__ = 'Horea Christian'

from Bio import SeqIO

format_extensions = {
	"genbank": [".gbk", ".gb", ".genbank"],
	"fasta": [".fasta"]
}

def check_format(sequence_path, format):
	"""
	Checks if a file is in the given format, and if not converts its equivalent in a different format to the defined format
	"""
	from os.path import splitext

	_ , file_extension = splitext(sequence_path)
	if any(format_extension == file_extension for format_extension in format_extensions[format]):
		return sequence_path
	else:
		for try_format in format_extensions:
			if any(format_extension == file_extension for format_extension in format_extensions[try_format]):
				sequence_path = convert_seq(sequence_path, try_format, format)
				break
	return sequence_path

def write_seq(sequence_write_path=".cache/", entrez_id="", sequence="", sequence_id="", formats=["fasta"]):
	# CAREFUL! in case of multiple exports will only return the last exported (genbank) destination!

	# elif isinstance(sequence, SeqRecord.SeqRecord):
	# 	sequence = sequence.seq

	if sequence and not sequence_id:
		raise Exception("Please specify an ID for your sequence (whatever string you choose)")

	if entrez_id:
		from Bio import Entrez
		Entrez.email = "A.N.Other@example.com"     # Always tell NCBI who you are
		handle = Entrez.efetch(db="nucleotide", id=entrez_id, rettype="gb", retmode="text")
		record = SeqIO.read(handle, "genbank")
		if "fasta" in formats:
			destination = sequence_write_path+entrez_id+".fasta"
			SeqIO.write(record, destination, "fasta")
		if "genbank" in formats:
			destination = sequence_write_path+entrez_id+".gbk"
			SeqIO.write(record, destination, "genbank")

	if sequence:
		from Bio.Seq import Seq
		from Bio.Alphabet import generic_dna
		from Bio.SeqRecord import SeqRecord

		if isinstance(sequence, str):
			sequence = Seq(sequence, generic_dna)

		record=SeqRecord(sequence, id=sequence_id)
		if "fasta" in formats:
			destination = sequence_write_path+sequence_id+".fasta"
			SeqIO.write(record, destination, "fasta")
		if "genbank" in formats:
			destination = sequence_write_path+sequence_id+".gb"
			SeqIO.write(record, destination, "genbank")

	return destination

def convert_seq(file_path, from_format, to_format):
	"""
	Reads the sequence saved at `file_path` as `from_format` format, and creates a new file in thesame directory in a `to_format` format.
	Returns the path to the new file.
	"""

	from os.path import splitext

	record = SeqIO.read(file_path, from_format)
	file_path_base, _ = splitext(file_path)
	new_file_path = file_path_base+format_extensions[to_format][0]
	SeqIO.write(record, new_file_path, to_format)
	return new_file_path
