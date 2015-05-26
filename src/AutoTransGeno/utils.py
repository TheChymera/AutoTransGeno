__author__ = 'Horea Christian'

from Bio import SeqIO
from os.path import isfile

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
		record=sequence
		#
		# record=SeqRecord(sequence, id=sequence_id)
		# print record, type(record)
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


def check_fetch_record(data_dir, construct_name, formats=["genbank"]):
	if "genbank" in formats:
		main_record_file = data_dir+construct_name+".gbk"

		if not isfile(main_record_file):
			write_seq(sequence_write_path=data_dir, entrez_id=construct_name, formats=["genbank"])
	return main_record_file

def overhangs(enzyme):
	from Bio.Restriction import *

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

def extract_feature(sequence_id, data_dir, feature_names, write_file=False):
	# CAREFUL! only returns last detected sequence!

	sequence_path = check_fetch_record(data_dir, sequence_id, formats=["genbank"])
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
		extracted_feature_file = write_seq(sequence_write_path=data_dir, sequence=extracted_feature, sequence_id=feature_names[0])
	else:
		extracted_feature_file = ""

	return extracted_feature, extracted_feature_file

if __name__ == '__main__':
	# outp = find_primers(sequence_path="/home/chymera/data/CreBLseq/Cre_aj627603/preliminary_sequencing/hit_full.fasta")
	outp = extract_feature(sequence_id="AJ627603", data_dir="/home/chymera/data/CreBLseq/Cre_aj627603/", feature_names=["cre", "CRE", "Cre"])
