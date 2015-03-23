#!/usr/bin/env python

import numpy as np
from Bio.Blast.Applications import NcbiblastnCommandline
from Bio.Blast import NCBIXML
from StringIO import StringIO
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio.Graphics import GenomeDiagram
from reportlab.lib.units import cm
from Bio.Seq import Seq
from Bio import SeqIO
from Bio.Alphabet import generic_dna
from reportlab.lib import colors

def new_track(diagram, track_name):
	construct_track = diagram.new_track(
						1,
						name=track_name,
						greytrack=True,
						scale=True,
						scale_format="SInt",
						scale_smalltick_interval=100,
						scale_largetick_interval=1000,
						scale_fontsize=3,
						scale_ticks=True,
						scale_color=colors.black,
						scale_largeticks=0.4,
						scale_smallticks=0.2,
						axis_labels=True,
						greytrack_labels=1)
	track_features = construct_track.new_set()
	return construct_track, track_features

def add_to_track(track_features, query, subject, annotation="", feature_color=""):
	from Bio.Blast.Applications import NcbiblastnCommandline

	output = NcbiblastnCommandline(query=query, subject=subject, outfmt=5, task="blastn")()[0]
	blast_result_record = NCBIXML.read(StringIO(output))

	for alignment in blast_result_record.alignments:
		for hsp in alignment.hsps:
			if hsp.sbjct_start > hsp.sbjct_end:
				strand = -1
			else:
				strand = +1
			feature = SeqFeature(FeatureLocation(hsp.sbjct_start, hsp.sbjct_end), strand=strand)
			track_features.add_feature(feature,
				name=annotation+str(hsp.identities)+"/"+str(hsp.align_length),
				label=True, color=feature_color,
				label_color=feature_color,
				label_angle=0,
				sigil="ARROW")

datadir="/home/chymera/data/CreBLseq/Cre_aj627603/"
cre_record = SeqIO.read(datadir+"aj627603.gbk", 'gb')

gdd = GenomeDiagram.Diagram('Test Diagram', x=0.05, fragment_size=0.2)
cre_construct_track = gdd.new_track(1, name="aj627603", greytrack=True, greytrack_labels=1)

primer_track, primer_features = new_track(gdd, "aj627603")
add_to_track(primer_features, datadir+"preliminary_sequencing/primer_fw.fasta", datadir+"aj627603.fasta", annotation="fw_primer", feature_color=colors.green)
add_to_track(primer_features, datadir+"cre_fw.fasta", datadir+"aj627603.fasta", annotation="cre_fw", feature_color=colors.green)
add_to_track(primer_features, datadir+"cre_rv.fasta", datadir+"aj627603.fasta", annotation="cre_rv", feature_color=colors.green)

hit_track, hit_features = nre_track(gdd, "aj627603")
gds_features = cre_construct_track.new_set()

# ########## Run BLAST and add features for hit alignment
output = NcbiblastnCommandline(query=datadir+"preliminary_sequencing/hit_full.fasta", subject=datadir+"aj627603.fasta", outfmt=5)()[0]
blast_result_record = NCBIXML.read(StringIO(output))

for alignment in blast_result_record.alignments:
	print("a")
	for hsp in alignment.hsps:
		if hsp.sbjct_start > hsp.sbjct_end:
			strand = -1
		else:
			strand = +1
		feature = SeqFeature(FeatureLocation(hsp.sbjct_start, hsp.sbjct_end), strand=strand)
		hit_features.add_feature(feature, name="hit_full_match "+str(hsp.identities)+"/"+str(hsp.align_length), label=True, color=colors.red, height=5)
##########

for feature in cre_record.features:
	if feature.type != "gene":
		#Exclude this feature
		continue
	if len(gds_features) % 2 == 0:
		color = colors.thistle
	else:
		color = colors.lightblue
	gds_features.add_feature(feature, sigil="ARROW", color=color, label_color=color, label=True, label_size = 14, label_angle=0, arrowshaft_height=1)

#I want to include some strandless features, so for an example
#will use EcoRI recognition sites etc.
restriction_list=[
	("GAATTC","EcoRI",colors.chartreuse),
	("CCCGGG","SmaI",colors.orange),
	("GTTAAC","HindI",colors.lemonchiffon),
	("GTCGAC","HindI",colors.cornsilk),
	("CAATC","HinfI",colors.olive),
	("CATTC","HinfI",colors.olive),
	("CAGTC","HinfI",colors.olive),
	("CACTC","HinfI",colors.olive),
	("AGCGCT","HaeII",colors.navy),
	("GGCGCC","HaeII",colors.navy),
	("GGATCC","BamHI",colors.purple)
	]
for site, name, color in restriction_list:
	index = 0
	while True:
		index  = cre_record.seq.find(site, start=index)
		if index == -1 : break
		feature = SeqFeature(FeatureLocation(index, index+len(site)))
		gds_features.add_feature(feature, color=color, name=name, label=True, label_size=3, label_position="end", label_color=color, label_angle=90)
		index += len(site)

gdd.draw(format='linear', pagesize="A4", fragments=1, start=0, end=len(cre_record))
gdd.write("output/plasmid_linear_nice.pdf", "PDF")
