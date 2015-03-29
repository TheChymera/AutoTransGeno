#!/usr/bin/env python
__author__ = 'Horea Christian'

from Bio import SeqIO

format_extensions = {
	"genbank": ".gbk",
	"fasta": ".fasta"
}

def write_seq(sequence_write_path=".cache/", entrez_id="", sequence="", sequence_id="", formats=["fasta"]):

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
		record=SeqRecord(Seq(sequence, generic_dna), id=sequence_id)
		if "fasta" in formats:
			destination = sequence_write_path+entrez_id+".fasta"
			SeqIO.write(record, destination, "fasta")
		if "genbank" in formats:
			destination = sequence_write_path+entrez_id+".gb"
			SeqIO.write(record, destination, "genbank")

	return destination

def convert_seq(file_path, from_format, to_format):
	from os.path import splitext

	record = SeqIO.read(file_path, from_format)
	file_path_base, _ = splitext(file_path)
	new_file_path = file_path_base+format_extensions[to_format]
	SeqIO.write(record, new_file_path, to_format)
	return new_file_path
