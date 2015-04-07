#!/usr/bin/env python

from Bio.Emboss.Applications import PrimerSearchCommandline


genbank_local_dir = "/home/chymera/genbank/"
seqence_subdir = "genomes/Eukaryotes/vertebrates_mammals/Mus_musculus/GRCm38.p3/"
assembly_subdir = "Primary_Assembly/assembled_chromosomes/FASTA/"

sequence_path = genbank_local_dir + seqence_subdir + assembly_subdir
