__author__ = "Horea Christian"

from Bio import Restriction, SeqRecord
from Bio.Seq import Seq
from StringIO import StringIO
from utils import write_seq, check_format
from Bio import SeqIO
from Bio.SeqUtils import GC
from Bio import Restriction

def select_digest(sequence, rb, max_length=10000, min_length=100, max_gc=100, min_gc=0, save_to=False):
	"""
	Usage Example:
		template = SeqIO.read("/home/chymera/genbank/genomes/Eukaryotes/vertebrates_mammals/Mus_musculus/GRCm38.p3/Primary_Assembly/assembled_chromosomes/FASTA/chr1.fa", "fasta")
		results = select_digest(template, "XceI", min_length=1500, max_length=1500, max_gc=55, min_gc=30, save_to="digests/")
	"""

	if isinstance(sequence, str):
		sequence = Seq(sequence)
	elif isinstance(sequence, SeqRecord.SeqRecord):
		sequence = sequence.seq

	digest = getattr(Restriction, rb).catalyze(sequence)

	selected_fragments = []
	for ix, fragment in enumerate(digest):
		if min_length <= len(fragment) <= max_length and min_gc <= GC(fragment) <= max_gc:
			selected_fragments += [fragment]
			if save_to:
				import time
				from os import makedirs, path
				append_date = "/"+time.strftime("%Y%m%d_%H%M")+"/"
				save_path = save_to + append_date
				if not path.exists(save_path):
					makedirs(save_path)
				write_seq(sequence_write_path=save_path, entrez_id="", sequence=fragment, sequence_id=str(len(fragment))+"-"+str(ix), formats=["fasta"])

	return selected_fragments

def enzyme_selector(sequence, restriction_interval, genome_frequency=False, deterministic_overhangs=False, rb=False):
	"""
	Usage Example:
		from sequence_utils import extract_feature
		sequence,_ = extract_feature(sequence_id="AJ627603", data_dir="/home/chymera/data2/gt.ep/sequences/", feature_names=["Cre", "cre", "CRE"])
		outp = enzyme_selector(sequence=sequence, restriction_interval=[0,690], genome_frequency=[700,2000], deterministic_overhangs=True)
		print outp
	"""


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
