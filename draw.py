#!/usr/bin/env python

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

gdd = GenomeDiagram.Diagram('Test Diagram', x=0.01, fragment_size=0.5)
cre_construct_track = gdd.new_track(1, greytrack=True)
alignment_track = gdd.new_track(1, greytrack=True)
gds_features = cre_construct_track.new_set()
gds_features1 = alignment_track.new_set(height=5)

########## Run BLAST and add features for fw_primer
output = NcbiblastnCommandline(query=datadir+"preliminary_sequencing/primer_fw.fasta", subject=datadir+"aj627603.fasta", outfmt=5, out=".cache/last_BLAST.xml")
output = open(".cache/last_BLAST.xml")
blast_result_record = NCBIXML.read(output)

for alignment in blast_result_record.alignments:
	for hsp in alignment.hsps:
		feature = SeqFeature(FeatureLocation(hsp.sbjct_start, hsp.sbjct_end), strand=+1)
		gds_features.add_feature(feature, name="fw_primer "+str(hsp.identities)+"/"+str(hsp.align_length), label=True, color=colors.green)
##########

########## Run BLAST and add features for fw_primer
output = NcbiblastnCommandline(query=datadir+"cre_fw.fasta", subject=datadir+"aj627603.fasta", outfmt=5, out=".cache/last_BLAST.xml")
output = open(".cache/last_BLAST.xml")
blast_result_record = NCBIXML.read(output)

print(len(blast_result_record.alignments))
for alignment in blast_result_record.alignments:
	for hsp in alignment.hsps:
		feature = SeqFeature(FeatureLocation(hsp.sbjct_start, hsp.sbjct_end), strand=+1)
		gds_features.add_feature(feature, name="cre_fw_primer "+str(hsp.identities)+"/"+str(hsp.align_length), label=True, color=colors.green)
##########

########## Run BLAST and add features for hit alignment
output = NcbiblastnCommandline(query=datadir+"preliminary_sequencing/hit_full.fasta", subject=datadir+"aj627603.fasta", outfmt=5)()[0]
blast_result_record = NCBIXML.read(StringIO(output))

# Print some information on the result
for alignment in blast_result_record.alignments:
	print("a")
	for hsp in alignment.hsps:
		print(hsp)
		print(hsp.query_end)
		print(hsp.query_start)
		print(len(hsp.query))
		feature = SeqFeature(FeatureLocation(hsp.sbjct_start, hsp.sbjct_end), strand=+1)
		gds_features1.add_feature(feature, name="hit_full_match "+str(hsp.identities)+"/"+str(hsp.align_length), label=True, color=colors.red, height=5)
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
		color = colors.lemonchiffon
	else:
		color = colors.lavenderblush
	gds_features.add_feature(feature, sigil="ARROW", color=color, label=True, label_size = 14, label_angle=0,arrowshaft_height=1)

#I want to include some strandless features, so for an example
#will use EcoRI recognition sites etc.
restriction_list=[
	("GAATTC","EcoRI",colors.green),
	("CCCGGG","SmaI",colors.orange),
	("GTTAAC","HindI",colors.brown),
	("GTCGAC","HindI",colors.brown),
	("CAATC","HinfI",colors.pink),
	("CATTC","HinfI",colors.pink),
	("CAGTC","HinfI",colors.pink),
	("CACTC","HinfI",colors.pink),
	("AGCGCT","HaeII",colors.blue),
	("GGCGCC","HaeII",colors.blue),
	("GGATCC","BamHI",colors.purple)
	]
for site, name, color in restriction_list:
	index = 0
	while True:
		index  = cre_record.seq.find(site, start=index)
		if index == -1 : break
		feature = SeqFeature(FeatureLocation(index, index+len(site)))
		gds_features.add_feature(feature, color=color, name=name, label=True, label_size = 10, label_color=color)
		index += len(site)

gdd.draw(format='linear', pagesize=(len(cre_record)/100*cm, 10	*cm), fragments=1, start=0, end=len(cre_record))
gdd.write("output/plasmid_linear_nice.pdf", "PDF")
