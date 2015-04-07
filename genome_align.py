#!/usr/bin/env python
__author__ = "Horea Christian"

import pickle
from Bio import Entrez
Entrez.email = "h.chr@mail.ru"
from Bio import pairwise2
from Bio.pairwise2 import format_alignment

from Bio.Blast import NCBIWWW
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
from Bio.Blast import NCBIXML
from Bio import pairwise2

new=False 	 # fetch new alignment from ncbi (takes a lot of time)
printalign=False
conditional=False # print only hits above a certain e-value
database="GPIPE/10090/105/ref_top_level"
datadir="/home/chymera/data/CreBLseq/Cre_aj627603/"
cre_construct="AJ627603.txt"

my_fw_seq = Seq("AATATGGCGAGGAACACTGAAAAATTGGAAAATTTANATGTCCACTGTANGACATGGAATATGGCAAGAAAACTGAAAATCATGGAAAATGANAAACATCCACGTGACGACTTGAAAAATGACGAAATCACTGAAAAACGTGAAAAATGANAAATGCACACTCTCGGACCTGGAATATGGCGAGAAAACTGAAAATCACGGAAAATGAGAAATACACACTTTANGACGTGAAATATGGANAGGAAAACTGAAAAAGTGGAAAATTTANAAATGTCCACTGTAGGACTGTTGGGAGCCGCGCGAGATATGGCA", generic_dna) #len 312
my_rv_seq = Seq("TGCGCCTGCTGGAAGATGGCGAGAAAACTGAAAATCACGGAAAATGAGAAATACACACTTTAGGACGTGAAATATGGCGAGGAACACTGAAAAATTGGAAAATTTAGATGTCCACTGTAGGACATGGAATATGGCAAGAAAACTGAAAATCATGGAAAATGAGAAACATCCACGTGACGACTTGAAAAATGACGAAATCACTGAAAAACGTGAAAAATGAGAAATGCACACTCTCGGACCTGGAATATGGCGAGAAAACTGAAAATCACGGAAAANGAGAAATACACACTTTAGGACGTGAAATATGGAGAGGAAAAC", generic_dna) #len 318
my_seq = Seq("TGCGCCTGCTGGAAGATGGCGAGAAAACTGAAAATCACGGAAAATGAGAAATACACACTTTAGGACGTGAAATATGGCGAGGAACACTGAAAAATTGGAAAATTTANATGTCCACTGTANGACATGGAATATGGCAAGAAAACTGAAAATCATGGAAAATGANAAACATCCACGTGACGACTTGAAAAATGACGAAATCACTGAAAAACGTGAAAAATGANAAATGCACACTCTCGGACCTGGAATATGGCGAGAAAACTGAAAATCACGGAAAATGAGAAATACACACTTTANGACGTGAAATATGGANAGGAAAACTGAAAAAGTGGAAAATTTANAAATGTCCACTGTAGGACTGTTGGGAGCCGCGCGAGATATGGCA", generic_dna) #len 388

# print my_seq

primer_fw = SeqIO.parse(datadir+"preliminary_sequencing/primer_fw.txt", "gb")

f=open(datadir+"preliminary_sequencing/primer_fw.txt", 'r')
primer_fw = Seq(f.read(), generic_dna)
f.close()

f=open(datadir+"AJ627603.txt", 'r')
cre = Seq(f.read(), generic_dna)
f.close()

alignments = pairwise2.align.globalms(cre, primer_fw, 2, -1, -.5, -.1)
print("a")

for alignment in alignments:
	print("a")
	print(format_alignment(*alignment))

if new:
	result_handle = NCBIWWW.qblast("blastn", database, my_seq, megablast=False, alignments=1000, hitlist_size=100, expect=10.0)
	save_file = open("my_blast.xml", "w")
	save_file.write(result_handle.read())
	save_file.close()
	result_handle.close()

result_handle = open("my_blast.xml")

blast_records = NCBIXML.read(result_handle)


E_VALUE_THRESH = 1e-115

alignment_nr = 0
hsp_nr = 0

if printalign:
	print(blast_records)
	for alignment in blast_records.alignments:
		alignment_nr += 1
		for hsp in alignment.hsps:
			if conditional:
				if hsp.expect < E_VALUE_THRESH:
					if "musculus" in alignment.title or "Mus" in alignment.title or "mouse" in alignment.title:
						if not "chromosome Y" in alignment.title and not "chromosome y" in alignment.title:
							hsp_nr += 1
							print(hsp)
							print(alignment.title)
		 					print("\n")
			else:
				hsp_nr += 1
				print(hsp)
				print(alignment.title)
				print("\n")
	print alignment_nr
	print hsp_nr
