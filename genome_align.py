#!/usr/bin/env python
from __future__ import division
__author__ = 'Horea Christian'

from Bio import Entrez
Entrez.email = "h.chr@mail.ru"
from Bio import pairwise2
from Bio.pairwise2 import format_alignment

from Bio.Blast import NCBIWWW
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
from Bio.Blast import NCBIXML

my_seq = Seq("AATATGGCGAGGAACACTGAAAAATTGGAAAATTTANATGTCCACTGTANGACATGGAATATGGCAAGAAAACTGAAAATCATGGAAAATGANAAACATCCACGTGACGACTTGAAAAATGACGAAATCACTGAAAAACGTGAAAAATGANAAATGCACACTCTCGGACCTGGAATATGGCGAGAAAACTGAAAATCACGGAAAATGAGAAATACACACTTTANGACGTGAAATATGGANAGGAAAACTGAAAAAGTGGAAAATTTANAAATGTCCACTGTAGGACTGTTGGGAGCCGCGCGAGATATGGCA", generic_dna)

result_handle = NCBIWWW.qblast("blastn", "nt", my_seq)
blast_records = NCBIXML.read(result_handle)

E_VALUE_THRESH = 1e-100

for alignment in blast_records.alignments:
	for hsp in alignment.hsps:
		if hsp.expect < E_VALUE_THRESH:
			print('****Alignment****')
			print('sequence:', alignment.title)
			print('length:', alignment.length)
			print('e value:', hsp.expect)
			print(hsp.query[0:75] + '...')
			print(hsp.match[0:75] + '...')
			print(hsp.sbjct[0:75] + '...')
