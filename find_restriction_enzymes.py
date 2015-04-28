from Bio import SeqIO

def enzyme_selector(sequence, restriction_interval, genome_frequency=False, deterministic_overhangs=False, rb=False):
	from Bio.Restriction import Analysis, AllEnzymes, RestrictionBatch

	if not rb:
		basic_analysis = Analysis(AllEnzymes, sequence.seq)
	else:
		basic_analysis = Analysis(rb, sequence.seq)
	respect_target = basic_analysis.only_between(restriction_interval[0],restriction_interval[1])
	# print respect_target

	if genome_frequency:
		respect_frequency = respect_target
	 	for enzyme, item in respect_target.items():
			if enzyme.frequency() < genome_frequency[0] or enzyme.frequency() > genome_frequency[1]:
				del respect_frequency[enzyme]
			else:
				if deterministic_overhangs:
					from sequence_utils import overhangs
					if any(bp_ID in overhangs(enzyme) for bp_ID in ["N", "R", "Y", "!!!", "S", "W", "M", "K", "B", "D", "H", "V"]) or overhangs(enzyme) == "":
						del respect_frequency[enzyme]

	return respect_frequency

if __name__ == '__main__':
	from sequence_utils import extract_feature

	sequence,_ = extract_feature(sequence_id="AJ627603", data_dir="/home/chymera/data2/gt.ep/sequences", feature_names=["Cre", "cre", "CRE"])
	outp = enzyme_selector(sequence=sequence, restriction_interval=[0,690], genome_frequency=[700,2000], deterministic_overhangs=True)
	print outp
