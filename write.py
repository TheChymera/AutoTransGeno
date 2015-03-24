#!/usr/bin/env python
__author__ = 'Horea Christian'

from Bio import SeqIO

def write_seq(sequence_write_path=".cache/", entrez_id="", sequence="", sequence_id=""):
	sequence_write_path="/home/chymera/data/CreBLseq/Cre_aj627603/preliminary_sequencing/"

	if sequence and not sequence_id:
		raise Exception("Please specify an ID for your sequence (whatever string you choose)")

	if entrez_id:
		from Bio import Entrez
		Entrez.email = "A.N.Other@example.com"     # Always tell NCBI who you are
		handle = Entrez.efetch(db="nucleotide", id=entrez_id, rettype="gb", retmode="text")
		record = SeqIO.read(handle, "genbank")
		destination = sequence_write_path+entrez_id+".fasta"
		SeqIO.write(record, destination, "fasta")

	if sequence_list:
		from Bio.Seq import Seq
		from Bio.Alphabet import generic_dna
		from Bio.SeqRecord import SeqRecord
		record=SeqRecord(Seq(sequence, generic_dna), id=sequence_id)
		destination = sequence_write_path+sequence_id+".fasta"
		SeqIO.write(record, destination, "fasta")

	return destination
