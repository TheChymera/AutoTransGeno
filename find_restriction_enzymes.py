from Bio import SeqIO

def enzyme_selector(sequence, restriction_interval, genome_frequency=False, deterministic_overhangs=False):
	from Bio.Restriction import *

	basic_analysis = Analysis(AllEnzymes, sequence.seq)
	respect_target = basic_analysis.only_between(restriction_interval[0],restriction_interval[1])
	# print respect_target

	if genome_frequency:
		respect_frequency = respect_target
	 	for enzyme, item in respect_target.items():
			if enzyme.frequency() < genome_frequency[0] or enzyme.frequency() > genome_frequency[1]:
				del respect_frequency[enzyme]
			else:
				if deterministic_overhangs:
					if "N" in enzyme.elucidate() or "R" in enzyme.elucidate() or "Y" in enzyme.elucidate():
						del respect_frequency[enzyme]

	return respect_frequency

if __name__ == '__main__':
	from sequence_utils import extract_feature

	sequence,_ = extract_feature(sequence_id="AJ627603", data_dir="/home/chymera/data/CreBLseq/Cre_aj627603/", feature_names=["Cre", "cre", "CRE"])
	outp = enzyme_selector(sequence=sequence, restriction_interval=[0,400], genome_frequency=[500,2000], deterministic_overhangs=True)