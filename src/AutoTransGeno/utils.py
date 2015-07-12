__author__ = 'Horea Christian'

from Bio.Restriction import *
from Bio.pairwise2 import format_alignment
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO, pairwise2
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
from os.path import isfile
from Bio.Alphabet import IUPAC
import numpy as np
import random
import string


format_extensions = {
	"genbank": [".gbk", ".gb", ".genbank"],
	"fasta": [".fasta"]
}

def simple_sequence_merge(sequences=[]):

	#make sure this works with string and sequence object inputs:
	if isinstance(sequences[0], str):
		pass
	else:
		sequences = [str(sequence.seq) for sequence in sequences]

	sequences = sorted(sequences, key=len, reverse=True)
	sequences_number = len(sequences)
	sequences_index = range(sequences_number)
	merged_sequence = ""
	for idx, position in enumerate(sequences[0]):
		added = False
		for character in ["-", "N"]:
			merged_sequence, added = check_position(position, character, merged_sequence, idx, sequences_index, sequences)
			if added:
				break
		if not added:
			merged_sequence += position

	return merged_sequence

def check_position(position, character, merged_sequence, idx, sequences_index, sequences, added=False):
	if position == character:
		for i in sequences_index[1:]:
			try:
				if sequences[i][idx] != character:
					merged_sequence += sequences[i][idx]
					added = True
			except IndexError:
				merged_sequence += position
				added = True
	return merged_sequence, added

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

def write_seq(sequence_write_path="/tmp/", entrez_id="", sequence="", ID="", record="", formats=["fasta"]):
	# CAREFUL! in case of multiple exports will only return the last exported (genbank) destination!

	# giving the sequence a random ID if not otherwise set
	if not ID:
		ID = ''.join(random.SystemRandom().choice(string.ascii_uppercase + string.digits) for _ in range(8))

	if sequence:
		if "fasta" in formats:
			destination = sequence_write_path+ID+".fasta"
			record = SeqRecord(Seq(sequence,IUPAC.ambiguous_dna), id="Unnamed sequence, outputted by AutoTransGeno")
			SeqIO.write(record, destination, "fasta")
		if "genbank" in formats:
			destination = sequence_write_path+ID+".gbk"
			record = SeqRecord(Seq(sequence,IUPAC.ambiguous_dna), id="Unnamed sequence, outputted by AutoTransGeno")
			SeqIO.write(record, destination, "genbank")

	if record:
		if "fasta" in formats:
			destination = sequence_write_path+ID+".fasta"
			SeqIO.write(record, destination, "fasta")
		if "genbank" in formats:
			destination = sequence_write_path+ID+".gbk"
			SeqIO.write(record, destination, "genbank")

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

def check_fetch_record(data_dir, construct_name, formats=["genbank"]):
	if "genbank" in formats:
		main_record_file = data_dir+construct_name+".gbk"

		if not isfile(main_record_file):
			main_record_file = write_seq(sequence_write_path=data_dir, entrez_id=construct_name, formats=["genbank"])
	return main_record_file

def overhangs(enzyme):

	if "?" in enzyme.elucidate() or "cut twice, not yet implemented sorry" in enzyme.elucidate():
		return("!!!NON-FATAL ERROR: UNDETERMINED OVERHANGS")
	else:
		split_5 = enzyme.elucidate().split("^")
		if "_" in split_5[0]:
			overhang = split_5[0].split("_")[1]
		elif "_" in split_5[1]:
			overhang = split_5[1].split("_")[0]
		return overhang

def standard_template(template):
	if "." in template:
		template = check_format(template, "fasta")
	elif type(template) == str:
		dest = write_seq(sequence=template, sequence_id="template")
		template = dest
	else:
		raise Exception("Pleae check your specified `template` value.")
	return template

def concatenate_overlapping_sequences(sequence1, sequence2):
	# add high penalties for opening gaps (we want no more than two, and low penalties for extending them - they should both be long)
	# e.g:
	#ATGATAGATAGCCAGATA------------------
	#||||||||||||||||||||||||||||||||||||
	#--------------GATATAACGTTANATTAWSCAT

	alignment=pairwise2.align.globalms(sequence1, sequence2, 1, -4, -90, 0)

	# if len(alignment) > 1:
		# raise Exception("There are multiple global alignment variants for your sequences. This should not happen and indicates that there is a flaw in our algorithm. Contact h.chr@mai.ru")

	for a in alignment:
		aligned_sequences = list(a)[:2]

	initial_gap_length = 0
	for a,b in zip(*aligned_sequences):
		if a == "-":
			initial_gap_length += 1
		else:
			break
	aligned_sequences = [i[initial_gap_length:] for i in aligned_sequences]
	concatenated_sequences = simple_sequence_merge(aligned_sequences)
	return concatenated_sequences

def extract_feature(sequence_id, data_dir, feature_names, write_file=False):
	# CAREFUL! only returns last detected sequence!

	sequence_path = check_fetch_record(data_dir+"/genbank/", sequence_id, formats=["genbank"])
	main_record = SeqIO.read(sequence_path, 'gb')

	for feature in main_record.features:
		if feature.type != "gene" and feature.type != "regulatory":
			# Exclude nongene features
			continue
		hit = 0
		for feature_name in feature_names:
			if feature_name in str(feature.qualifiers):
				hit += 1
		if hit > 0:
			extracted_feature = feature.extract(main_record)

	if write_file:
		extracted_feature_file = write_seq(sequence_write_path=data_dir, record=extracted_feature, ID=feature_names[0]+"_from_"+sequence_id)
	else:
		extracted_feature_file = ""

	return extracted_feature, extracted_feature_file
