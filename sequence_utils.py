from Bio import SeqIO
from os.path import isfile

def check_fetch_record(data_dir, construct_name, formats=["genbank"]):
	if "genbank" in formats:
		main_record_file = data_dir+construct_name+".gbk"

		if not isfile(main_record_file):
			from write import write_seq
			write_seq(sequence_write_path=datadir, entrez_id=construct_name, formats=["genbank"])
	return main_record_file


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