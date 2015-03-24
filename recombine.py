#!/usr/bin/env python
__author__ = 'Horea Christian'

from Bio import SeqIO
from Bio.Restriction import *

sequence_write_path="/home/chymera/data/CreBLseq/Cre_aj627603/preliminary_sequencing/"
sequence_read_path=""
entrez_id=""
# Please enter a list of ids for your sequences to use as filenames (you can pick whatever names you like)
sequence_list=[]
sequence_list_ids=[]
file_paths=["hit_full.fasta"]

if not sequence_write_path:
	sequence_write_path=".cache/"
if not sequence_read_path:
	sequence_read_path = sequence_write_path

if sequence_list:
	from Bio.Seq import Seq
	from Bio.Alphabet import generic_dna
	from Bio.SeqRecord import SeqRecord
	from numpy import arange
	for nr in arange(len(sequence_list)):
		sequence = Seq(sequence_list[nr], generic_dna)
		Ana = Analysis(CommOnly, sequence, linear=False)
		Ana.print_that()

if file_paths:
	from Bio.Seq import Seq
	from Bio.SeqRecord import SeqRecord
	from numpy import arange
	for nr in arange(len(file_paths)):
		record = SeqIO.read(sequence_read_path+file_paths[nr], "fasta")
		sequence = record.seq
		Ana = Analysis(CommOnly, sequence, linear=False)
		Ana.print_that()
