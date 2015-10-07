__author__ = "Horea Christian"

from Bio.Blast.Applications import NcbiblastnCommandline
from Bio.Seq import Seq
from Bio.Blast import NCBIXML
from StringIO import StringIO
from utils import check_format
from utils import standard_template
import itertools
from Bio import SeqIO
import os

genbank_local_dir = "/home/chymera/genbank/"
seqence_subdir = "genomes/Eukaryotes/vertebrates_mammals/Mus_musculus/GRCm38.p3/"
assembly_subdir = "Primary_Assembly/assembled_chromosomes/FASTA/"

sequence_path = genbank_local_dir + seqence_subdir + assembly_subdir

def mock_pcr(template, primer_files=[], primer_seqs="", primer_directory="", max_length=3000, evalue=0.1):
	"""
	Performs a rudimentary in-silico PCR, which returns the sequence and length of template fragments flanked by primers shother than max_length.
	"""
	if primer_files == "" and primer_seqs == "":
		raise Exception("Please specify at least two entries for `primer_seqs`, `primer_files`, or one for each.")
	if len(primer_files)+len(primer_seqs) < 2:
		raise Exception("Please specify at least two primers.")

	#converts files from whatever format they have into FASTA (required for ncbi blast)
	if primer_files:
		for ix, primer_file in enumerate(primer_files):
			new_primer_file = check_format(primer_directory+primer_file, "fasta")
			primer_files[ix] = new_primer_file

	#convert template format if necessary, if template not formatte a path write chached file for seuence
	template = standard_template(template)

	#writes sequences to files
	if primer_seqs:
		for ix, primer_seq in enumerate(primer_seqs):
			if type(primer_seq) == str:
				dest = write_seq(sequence=primer_seq, sequence_id="primer"+str(ix))
				primer_files += [dest]

	fw_primer_start_positions=[]
	rv_primer_start_positions=[]
	template_seq = SeqIO.read(template, "fasta")
	for primer_file in primer_files:
		output = NcbiblastnCommandline(query=template, subject=primer_file, outfmt=5, evalue=evalue, task="blastn")()[0]
		blast_result_record = NCBIXML.read(StringIO(output))

		for alignment in blast_result_record.alignments:
			hsp_list=alignment.hsps
			for hsp in hsp_list:
				# print "hsp infos:", hsp.sbjct_start, hsp.sbjct_end, hsp.query_start, hsp.query_end, hsp.align_length, hsp.sbjct
				if hsp.sbjct_start > hsp.sbjct_end:
					strand = -1
					amplicon_5_end = hsp.query_end
					rv_primer_start_positions += [hsp.query_end]
				else:
					strand = +1
					amplicon_5_end = hsp.query_start
					fw_primer_start_positions += [hsp.query_start]

	for fw_primer_start in fw_primer_start_positions:
		for rv_primer_start in rv_primer_start_positions:
			amplicon_length = rv_primer_start - fw_primer_start
			if amplicon_length <= max_length and amplicon_length >= 0:
				return amplicon_length, template_seq[fw_primer_start:rv_primer_start].seq


def mock_pcr_star(a_b_c_d_e):
	"""Convert `mock_pcr([1,2,3,4,5,6])` to `mock_pcr(1,2,3,4,5,6)` call."""
	return mock_pcr(*a_b_c_d_e)

def genome_mock_pcr(genome_path, unspecific_insert=False, primer_files="", primer_seqs="", primer_directory="", max_length=3000, evalue=0.1, poolsize=8):
	from multiprocessing.dummy import Pool as ThreadPool
	pool = ThreadPool(poolsize)

	fasta_files = [genome_path+fasta_file for fasta_file in os.listdir(genome_path) if os.path.splitext(fasta_file)[-1] in [".fa",".fna","fasta"]]

	# for i in (zip(fasta_files,
	# 	itertools.repeat(primer_files),
	# 	itertools.repeat(primer_seqs),
	# 	itertools.repeat(primer_directory),
	# 	itertools.repeat(max_length))):
		# print i

	amplicons = pool.map(mock_pcr_star,
						zip(fasta_files,
							itertools.repeat(primer_files),
							itertools.repeat(primer_seqs),
							itertools.repeat(primer_directory),
							itertools.repeat(max_length),
							itertools.repeat(evalue)))

	#close the pool and wait for the work to finish
	pool.close()
	pool.join()

	if unspecific_insert:
		amplicons += [mock_pcr(unspecific_insert, primer_files=primer_files, primer_seqs=primer_seqs, primer_directory=primer_directory, max_length=max_length)]

	return amplicons
