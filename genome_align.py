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

# with open("/home/chymera/data/CreBLseq/Cre_aj627603/hit_fw.txt") as handle:
# 	my_seq = Seq(handle.read().strip(), generic_dna)

my_seq = Seq("AATATGGCGAGGAACACTGAAAAATTGGAAAATTTANATGTCCACTGTANGACATGGAATATGGCAAGAAAACTGAAAATCATGGAAAATGANAAACATCCACGTGACGACTTGAAAAATGACGAAATCACTGAAAAACGTGAAAAATGANAAATGCACACTCTCGGACCTGGAATATGGCGAGAAAACTGAAAATCACGGAAAATGAGAAATACACACTTTANGACGTGAAATATGGANAGGAAAACTGAAAAAGTGGAAAATTTANAAATGTCCACTGTAGGACTGTTGGGAGCCGCGCGAGATATGGCA", generic_dna)

result_handle = NCBIWWW.qblast("blastn", "nt", my_seq)
blast_records = NCBIXML.read(result_handle)
# for i in result_handle:
# 	print(i)
# print(result_handle.getvalue())
# print(blast_records)

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

# hit_fw="AATATGGCGAGGAACACTGAAAAATTGGAAAATTTANATGTCCACTGTANGACATGGAATATGGCAAGAAAACTGAAAATCATGGAAAATGANAAACATCCACGTGACGACTTGAAAAATGACGAAATCACTGAAAAACGTGAAAAATGANAAATGCACACTCTCGGACCTGGAATATGGCGAGAAAACTGAAAATCACGGAAAATGAGAAATACACACTTTANGACGTGAAATATGGANAGGAAAACTGAAAAAGTGGAAAATTTANAAATGTCCACTGTAGGACTGTTGGGAGCCGCGCGAGATATGGCA"
#
# hit_rv="TGCGCCTGCTGGAAGATGGCGAGAAAACTGAAAATCACGGAAAATGAGAAATACACACTTTAGGACGTGAAATATGGCGAGGAACACTGAAAAATTGGAAAATTTAGATGTCCACTGTAGGACATGGAATATGGCAAGAAAACTGAAAATCATGGAAAATGAGAAACATCCACGTGACGACTTGAAAAATGACGAAATCACTGAAAAACGTGAAAAATGAGAAATGCACACTCTCGGACCTGGAATATGGCGAGAAAACTGAAAATCACGGAAAANGAGAAATACACACTTTAGGACGTGAAATATGGAGAGGAAAAC"
#
# for a in pairwise2.align.localxx(hit_fw, hit_rv):
# 	print(format_alignment(*a))
