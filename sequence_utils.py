from Bio import SeqIO
from os.path import isfile

def check_fetch_record(data_dir, construct_name, formats=["genbank"]):
	if "genbank" in formats:
		main_record_file = data_dir+construct_name+".gbk"

		if not isfile(main_record_file):
			from write import write_seq
			write_seq(sequence_write_path=datadir, entrez_id=construct_name, formats=["genbank"])
	return main_record_file

def check_blast_format(sequence_path):
	from os.path import splitext

	_ , sequence_format = splitext(sequence_path)
	if sequence_format == ".gbk" or sequence_format == ".gb" or sequence_format == ".genbank":
		from write import convert_seq
		sequence_path = convert_seq(sequence_path, "genbank", "fasta")
	return sequence_path

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
		from write import write_seq
		extracted_feature_file = write_seq(sequence_write_path=data_dir, sequence=extracted_feature, sequence_id=feature_names[0])
	else:
		extracted_feature_file = ""

	return extracted_feature, extracted_feature_file

if __name__ == '__main__':
	# outp = find_primers(sequence_path="/home/chymera/data/CreBLseq/Cre_aj627603/preliminary_sequencing/hit_full.fasta")
	outp = extract_feature(sequence_id="AJ627603", data_dir="/home/chymera/data/CreBLseq/Cre_aj627603/", feature_names=["cre", "CRE", "Cre"])
