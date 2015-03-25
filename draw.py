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

def new_track(diagram, track_name, smalltick=100, largetick=1000):
	construct_track = diagram.new_track(
						1,
						name=track_name,
						greytrack=True,
						scale=True,
						scale_format="SInt",
						scale_smalltick_interval=smalltick,
						scale_largetick_interval=largetick,
						scale_fontsize=3,
						scale_ticks=True,
						scale_color=colors.black,
						scale_largeticks=0.4,
						scale_smallticks=0.2,
						axis_labels=True,
						greytrack_labels=1)
	track_features = construct_track.new_set()
	return construct_track, track_features

def add_to_track(track_features, query, subject, annotation="", feature_color="", forceone=False, sigil="ARROW", label_angle=0):
	from Bio.Blast.Applications import NcbiblastnCommandline

	subject = check_blast_format(subject)

	if forceone:
		output = NcbiblastnCommandline(query=query, subject=subject, outfmt=5, task="blastn", num_alignments=1)()[0]
	else:
		output = NcbiblastnCommandline(query=query, subject=subject, outfmt=5, task="blastn")()[0]
	blast_result_record = NCBIXML.read(StringIO(output))

	for alignment in blast_result_record.alignments:
		if forceone:
			hsp_list=alignment.hsps[:1]
		else:
			hsp_list=alignment.hsps
		for hsp in hsp_list:
			if hsp.sbjct_start > hsp.sbjct_end:
				strand = -1
			else:
				strand = +1
			feature = SeqFeature(FeatureLocation(hsp.sbjct_start, hsp.sbjct_end), strand=strand)
			track_features.add_feature(feature,
				name=annotation+" ("+str(hsp.identities)+"/"+str(hsp.align_length)+")",
				label=True, color=feature_color,
				label_color=feature_color,
				label_angle=label_angle,
				sigil=sigil,
				arrowshaft_height=1)

def multitrack(diagram, track_name, query, subject, smalltick=100, largetick=1000, annotation="", feature_color=""):
	from Bio.Blast.Applications import NcbiblastnCommandline

	output = NcbiblastnCommandline(query=query, subject=subject, outfmt=5, task="blastn")()[0]
	blast_result_record = NCBIXML.read(StringIO(output))

	for alignment in blast_result_record.alignments:
		for hsp in alignment.hsps:
			if hsp.sbjct_start > hsp.sbjct_end:
				strand = -1
			else:
				strand = +1
			construct_track = diagram.new_track(
								1,
								name=track_name,
								greytrack=True,
								scale=True,
								scale_format="SInt",
								scale_smalltick_interval=smalltick,
								scale_largetick_interval=largetick,
								scale_fontsize=3,
								scale_ticks=True,
								scale_color=colors.black,
								scale_largeticks=0.4,
								scale_smallticks=0.2,
								axis_labels=True,
								greytrack_labels=1)
			track_features = construct_track.new_set()
			feature = SeqFeature(FeatureLocation(hsp.sbjct_start, hsp.sbjct_end), strand=strand)
			track_features.add_feature(feature,
				name=annotation+" ("+str(hsp.identities)+"/"+str(hsp.align_length)+")",
				label=True, color=feature_color,
				label_color=feature_color,
				label_angle=0,
				sigil="ARROW",
				arrowshaft_height=1)

def draw_digest(track_features, main_record, restriction_list):
	for site, name, color in restriction_list:
		index = 0
		while True:
			index  = main_record.seq.find(site, start=index)
			if index == -1 : break
			feature = SeqFeature(FeatureLocation(index, index+len(site)))
			track_features.add_feature(feature, color=color, name=name, label=True, label_size=3, label_position="end", label_color=color, label_angle=90)
			index += len(site)

def check_blast_format(sequence_path):
	from os.path import splitext

	_ , sequence_format = splitext(sequence_path)
	if sequence_format == ".gbk" or sequence_format == ".gb" or sequence_format == ".genbank":
		from write import convert_seq
		sequence_path = convert_seq(sequence_path, "genbank", "fasta")
	return sequence_path

def cre(datadir, construct_name):
	from os.path import isfile

	main_record_file = datadir+construct_name+".gbk"
	# main_record_file

	if not isfile(main_record_file):
		from write import write_seq
		write_seq(sequence_write_path=datadir, entrez_id=construct_name, formats=["genbank", "fasta"])

	main_record = SeqIO.read(main_record_file, 'gb')

	gdd = GenomeDiagram.Diagram(construct_name+' Construct Diagram', x=0.05, track_size=0.2)

	genbank_track, genbank_features = new_track(gdd, construct_name+" features")
	for feature in main_record.features:
		if feature.type != "gene":
			# Exclude this feature
			continue
		if len(genbank_features) % 2 == 0:
			color = colors.thistle
		else:
			color = colors.lightblue
		genbank_features.add_feature(feature, sigil="ARROW", color=color, label_color=color, label=True, label_size = 14, label_angle=0, arrowshaft_height=1)

	restriction_track, restriction_features = new_track(gdd, construct_name+" restriction sites")
	draw_digest(restriction_features, main_record, restriction_list)

	primer_track, primer_features = new_track(gdd, construct_name+" primers")
	add_to_track(primer_features, datadir+"preliminary_sequencing/primer_fw.fasta", main_record_file, annotation="seq_fw_primer", feature_color=colors.green)
	add_to_track(primer_features, datadir+"cre_fw.fasta", main_record_file, annotation="cre_fw_primer", feature_color=colors.red)
	add_to_track(primer_features, datadir+"cre_rv.fasta", main_record_file, annotation="cre_rv_primer", feature_color=colors.red)


	hit_track, hit_features = new_track(gdd, construct_name+" hit alignments")
	add_to_track(hit_features, datadir+"preliminary_sequencing/hit_fw.fasta", main_record_file, annotation="alignement with hit_fw", feature_color=colors.cornflower, sigil="BOX", label_angle=30)
	add_to_track(hit_features, datadir+"preliminary_sequencing/hit_rv.fasta", main_record_file, annotation="alignement with hit_rv", feature_color=colors.salmon, sigil="BOX", label_angle=30)

	gdd.draw(format='linear', pagesize="A4", fragments=1, start=0, end=len(main_record))
	gdd.write("output/"+construct_name+"_diagram.pdf", "PDF")

def my_hit(datadir, construct_name):
	main_record = SeqIO.read(datadir+"preliminary_sequencing/hit_full.fasta", 'fasta')
	gdd = GenomeDiagram.Diagram(construct_name+' Construct Diagram', x=0.05, track_size=0.2)

	# multitrack(gdd, construct_name, datadir+"preliminary_sequencing/hit_fw.fasta", datadir+"preliminary_sequencing/hit_full.fasta", annotation="hit_fw", feature_color=colors.cornflower)
	# multitrack(gdd, construct_name, datadir+"preliminary_sequencing/hit_rv.fasta", datadir+"preliminary_sequencing/hit_full.fasta", annotation="ht_rv", feature_color=colors.salmon)

	hit_fw_track, hit_fw_features = new_track(gdd, construct_name+" forward hit", smalltick=10, largetick=100)
	add_to_track(hit_fw_features, datadir+"preliminary_sequencing/hit_fw.fasta", datadir+"preliminary_sequencing/hit_full.fasta", annotation="hit_fw", feature_color=colors.cornflower, forceone=True, sigil="BOX")

	hit_rv_track, hit_rv_features = new_track(gdd, construct_name+" reverse hit", smalltick=10, largetick=100)
	add_to_track(hit_rv_features, datadir+"preliminary_sequencing/hit_rv.fasta", datadir+"preliminary_sequencing/hit_full.fasta", annotation="hit_rv", feature_color=colors.salmon, forceone=True, sigil="BOX")

	primer_track, primer_features = new_track(gdd, construct_name+" primers", smalltick=10, largetick=100)
	add_to_track(primer_features, datadir+"preliminary_sequencing/primer_fw.fasta", datadir+"preliminary_sequencing/hit_full.fasta", annotation="seq_fw_primer", feature_color=colors.green)

	restriction_track, restriction_features = new_track(gdd, construct_name+" restriction sites", smalltick=10, largetick=100)
	for site, name, color in restriction_list:
		index = 0
		while True:
			index  = main_record.seq.find(site, start=index)
			if index == -1 : break
			feature = SeqFeature(FeatureLocation(index, index+len(site)))
			restriction_features.add_feature(feature, color=color, name=name, label=True, label_size=5, label_position="end", label_color=color, label_angle=90)
			index += len(site)

	gdd.draw(format='linear', pagesize="A4", fragments=1, start=1, end=len(main_record))
	gdd.write("output/"+construct_name+"_diagram.pdf", "PDF")

if __name__ == '__main__':
	restriction_list=[
		("GAATTC","EcoRI",colors.chartreuse),
		("CCCGGG","SmaI",colors.orange),
		("GTTAAC","HindI",colors.magenta),
		("GTCGAC","HindI",colors.purple),
		("GAATC","HinfI",colors.fuchsia),
		("GATTC","HinfI",colors.fuchsia),
		("GAGTC","HinfI",colors.fuchsia),
		("GACTC","HinfI",colors.fuchsia),
		("AGCGCT","HaeII",colors.cyan),
		("GGCGCC","HaeII",colors.cyan),
		("GGATCC","BamHI",colors.turquoise)
		]
	cre(datadir="/home/chymera/data/CreBLseq/Cre_aj627603/", construct_name="U80929")
	my_hit(datadir="/home/chymera/data/CreBLseq/Cre_aj627603/", construct_name="my_full_hit")
