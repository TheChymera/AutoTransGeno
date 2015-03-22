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

datadir="/home/chymera/data/CreBLseq/Cre_aj627603/"

cre_record = SeqIO.read(datadir+"aj627603.gbk", 'gb')

gdd = GenomeDiagram.Diagram('Test Diagram', x=0.05, fragment_size=0.2)
cre_construct_track = gdd.new_track(1, name="aj627603", greytrack=True, greytrack_labels=1)
primer_track = gdd.new_track(1, name="aj627603", greytrack=True, scale=True, scale_format="SInt", scale_smalltick_interval=100, scale_largetick_interval=1000, scale_fontsize=3, scale_ticks=True, scale_color=colors.black, scale_largeticks=0.4, scale_smallticks=0.2, axis_labels=True, greytrack_labels=1)
hit_track = gdd.new_track(1, name="aj627603", greytrack=True, greytrack_labels=1)
gds_features = cre_construct_track.new_set()
gds_features1 = primer_track.new_set(height=5)
gds_features2 = hit_track.new_set(height=5)

########## Run BLAST and add features for fw_primer
output = NcbiblastnCommandline(query=datadir+"preliminary_sequencing/primer_fw.fasta", subject=datadir+"aj627603.fasta", outfmt=5, task="blastn")()[0]
blast_result_record = NCBIXML.read(StringIO(output))

for alignment in blast_result_record.alignments:
	for hsp in alignment.hsps:
		if hsp.sbjct_start > hsp.sbjct_end:
			strand = -1
		else:
			strand = +1
		feature = SeqFeature(FeatureLocation(hsp.sbjct_start, hsp.sbjct_end), strand=strand)
		gds_features1.add_feature(feature, name="fw_primer "+str(hsp.identities)+"/"+str(hsp.align_length), label=True, color=colors.green, label_color=colors.green, label_angle=0, sigil="ARROW")
##########

########## Run BLAST and add features for cre_fw
output = NcbiblastnCommandline(query=datadir+"cre_fw.fasta", subject=datadir+"aj627603.fasta", outfmt=5, task="blastn")()[0]
blast_result_record = NCBIXML.read(StringIO(output))

print(len(blast_result_record.alignments))
for alignment in blast_result_record.alignments:
	for hsp in alignment.hsps:
		if hsp.sbjct_start > hsp.sbjct_end:
			strand = -1
		else:
			strand = +1
		feature = SeqFeature(FeatureLocation(hsp.sbjct_start, hsp.sbjct_end), strand=strand)
		gds_features1.add_feature(feature, name="cre_fw_primer "+str(hsp.identities)+"/"+str(hsp.align_length), label=True, color=colors.red, label_color=colors.red, label_angle=0, sigil="ARROW")
##########

########## Run BLAST and add features for cre_rv
output = NcbiblastnCommandline(query=datadir+"cre_rv.fasta", subject=datadir+"aj627603.fasta", outfmt=5, task="blastn")()[0]
blast_result_record = NCBIXML.read(StringIO(output))

print(len(blast_result_record.alignments))
for alignment in blast_result_record.alignments:
	for hsp in alignment.hsps:
		if hsp.sbjct_start > hsp.sbjct_end:
			strand = -1
		else:
			strand = +1
		feature = SeqFeature(FeatureLocation(hsp.sbjct_start, hsp.sbjct_end), strand=strand)
		gds_features1.add_feature(feature, name="cre_rv_primer "+str(hsp.identities)+"/"+str(hsp.align_length), label=True, color=colors.red, label_color=colors.red, label_angle=0, sigil="ARROW")
##########

# ########## Run BLAST and add features for hit alignment
output = NcbiblastnCommandline(query=datadir+"preliminary_sequencing/hit_full.fasta", subject=datadir+"aj627603.fasta", outfmt=5)()[0]
blast_result_record = NCBIXML.read(StringIO(output))

# Print some information on the result
for alignment in blast_result_record.alignments:
	print("a")
	for hsp in alignment.hsps:
		if hsp.sbjct_start > hsp.sbjct_end:
			strand = -1
		else:
			strand = +1
		feature = SeqFeature(FeatureLocation(hsp.sbjct_start, hsp.sbjct_end), strand=strand)
		gds_features2.add_feature(feature, name="hit_full_match "+str(hsp.identities)+"/"+str(hsp.align_length), label=True, color=colors.red, height=5)
		# feature = SeqFeature(FeatureLocation(hsp.sbjct_start-hsp.query_start, hsp.sbjct_end+hsp.query_start), strand=+1)
		# gds_features1.add_feature(feature, name="hit_full ", label=True, color=colors.pink, height=5)
##########

#Add three features to show the strand options,
# feature = SeqFeature(FeatureLocation(25, 1025), strand=+1)
# gds_features.add_feature(feature, name="Forward", label=True, color=colors.green)
# feature = SeqFeature(FeatureLocation(150, 250), strand=None)
# gds_features.add_feature(feature, name="Strandless", label=True)
# feature = SeqFeature(FeatureLocation(275, 375), strand=-1)
# gds_features.add_feature(feature, name="Reverse", label=True)


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
