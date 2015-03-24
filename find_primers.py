#!/usr/bin/env python
# -*- coding: UTF-8 -*-

from Bio import SeqIO
from Bio.Emboss.Applications import Primer3Commandline
from Bio.Seq import Seq

def find_primers(sequence="", sequence_path="", sequence_id=""):
	if not sequence and not sequence_path:
		raise Exception("Please specify a sequence or the path to a sequence!")
	if sequence and sequence_path:
		raise Exception("Please specify ONLY a sequence or the path to a sequence!")

	if sequence:
		from write import write_seq
		sequence_path = write_seq(sequence=sequence, sequence_id=sequence_id)
	primer_cl = Primer3Commandline(sequence=sequence_path,auto=True)
	primer_cl.outfile = "out.pr3"
	primer_cl.numreturn = 3
	primer_cl.osize = 20
	primer_cl.maxsize = 26
	primer_cl.opttm = 58
	primer_cl.mintm = 52
	primer_cl.mingc = 35
	primer_cl.maxgc = 75
	primer_cl.psizeopt = 200

	primer_cl()
