__author__ = "Horea Christian"

from Bio import SeqIO
from Bio.Emboss.Applications import Primer3Commandline
from Bio.Seq import Seq
from utils import write_seq

def find_primers(sequence="", sequence_path="", sequence_id="", output_path="", datadir=""):

	if not output_path:
		output_path = "output/last_primer_out.pr3"

	if sequence_id and not sequence and not sequence_path:
		sequence_path = write_seq(sequence_write_path=datadir, entrez_id=sequence_id, formats=["genbank", "fasta"])

	if sequence:
		sequence_path = write_seq(sequence=sequence, sequence_id=sequence_id)

	main_record = SeqIO.read(sequence_path, 'gb')

	for feature in main_record.features:
		if feature.type != "gene" and feature.type != "regulatory":
			# Exclude nongene features
			continue
		print feature.qualifiers
		if "cre" in str(feature.qualifiers) or "Cre" in str(feature.qualifiers):
			print feature.extract(main_record).seq

	primer_cl = Primer3Commandline(sequence=sequence_path,auto=True)
	primer_cl.outfile = output_path
	primer_cl.numreturn = 5
	primer_cl.osize = 20
	primer_cl.maxsize = 24
	primer_cl.opttm = 60
	primer_cl.mintm = 52
	primer_cl.mingc = 35
	primer_cl.maxgc = 75
	primer_cl.psizeopt = 200

	primer_cl()
	print(primer_cl)

	f = open(output_path, 'r')
	return f.read()

if __name__ == '__main__':
	outp = find_primers(sequence_path="/home/chymera/data/CreBLseq/Cre_aj627603/aj627603.gbk")
	# outp = extract_feature(sequence_id="AJ627603", data_dir="/home/chymera/data/CreBLseq/Cre_aj627603/", feature_names=["Cre", "cre", "CRE"])
	# print(outp)

	# from Bio.SeqUtils import MeltingTemp

	# print('%0.2f' % MeltingTemp.Tm_GC(Seq("GCGTTTTCCCTACATATGCTACCAG")))
