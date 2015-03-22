#!/usr/bin/env python
__author__ = 'Horea Christian'

from Bio import SeqIO

entrez_id=""
# Please enter a list of ids for your sequences to use as filenames (you can pick whatever names you like)
sequence_list=[]
sequence_list_ids=[]

if entrez_id:
	from Bio import Entrez
	Entrez.email = "A.N.Other@example.com"     # Always tell NCBI who you are
	handle = Entrez.efetch(db="nucleotide", id=entrez_id, rettype="gb", retmode="text")
	record = SeqIO.read(handle, "genbank")
	SeqIO.write(record, ".cache/"+entrez_id+".fasta", "fasta")

if sequence_list:
	from Bio.Seq import Seq
	from Bio.Alphabet import generic_dna
	from Bio.SeqRecord import SeqRecord
	for sequence_nr in enumerate(sequence_list):
		cre_fw=SeqRecord(Seq(sequence_list[nr], generic_dna), id=sequence_list_ids[nr])
		SeqIO.write(cre_fw, ".cache/"+sequence_list_ids[nr]+".fasta", "fasta")
