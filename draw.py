#!/usr/bin/env python
__author__ = "Horea Christian"

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

def new_track(diagram, track_name, smalltick=100, largetick="", track_no=1, endpoint=0):
	if not largetick:
		largetick=10*smalltick
	if endpoint != 0:
		construct_track = diagram.new_track(
							track_no,
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
							greytrack_labels=1,
							start=0,
							end=endpoint)
	else:
		construct_track = diagram.new_track(
							track_no,
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

def draw_digest(track_features, restriction_dict, restriction_colors=""):
	from itertools import cycle
	from sequence_utils import overhangs

	if not restriction_colors:
		restriction_colors=[
			colors.chartreuse,
			colors.deeppink,
			colors.turquoise,
			colors.fuchsia,
			colors.coral,
			colors.yellow
		]

	restriction_colors = cycle(restriction_colors)
	angle = cycle([120,90,60,30])
	for enzyme, sites in restriction_dict.iteritems():
		current_color = restriction_colors.next()
		current_angle=angle.next()
		overhang_length = len(overhangs(enzyme))
		for site in sites:
			feature = SeqFeature(FeatureLocation(site, site+overhang_length))
			track_features.add_feature(feature, color=current_color, name="      "+enzyme.__name__, label=True, label_size=3, label_position="end", label_color=current_color, label_angle=current_angle)


def check_blast_format(sequence_path):
	from os.path import splitext

	_ , sequence_format = splitext(sequence_path)
	if sequence_format == ".gbk" or sequence_format == ".gb" or sequence_format == ".genbank":
		from write import convert_seq
		sequence_path = convert_seq(sequence_path, "genbank", "fasta")
	return sequence_path

def my_hit(datadir, template_name):
	main_record = SeqIO.read(datadir+"preliminary_sequencing/hit_full.fasta", 'fasta')
	gdd = GenomeDiagram.Diagram(construct_name+' Construct Diagram', x=0.05, track_size=0.2)

	hit_fw_track, hit_fw_features = new_track(gdd, construct_name+" forward hit", smalltick=10, largetick=100)
	add_to_track(hit_fw_features, datadir+"preliminary_sequencing/hit_fw.fasta", datadir+"preliminary_sequencing/hit_full.fasta", annotation="hit_fw", feature_color=colors.cornflower, forceone=True, sigil="BOX")

	hit_rv_track, hit_rv_features = new_track(gdd, construct_name+" reverse hit", smalltick=10, largetick=100)
	add_to_track(hit_rv_features, datadir+"preliminary_sequencing/hit_rv.fasta", datadir+"preliminary_sequencing/hit_full.fasta", annotation="hit_rv", feature_color=colors.salmon, forceone=True, sigil="BOX")

	primer_track, primer_features = new_track(gdd, construct_name+" primers", smalltick=10, largetick=100)
	add_to_track(primer_features, datadir+"preliminary_sequencing/primer_fw.fasta", datadir+"preliminary_sequencing/hit_full.fasta", annotation="seq_fw_primer", feature_color=colors.green)

	gdd.draw(format='linear', pagesize="A4", fragments=1, start=0, end=len(main_record))
	gdd.write("output/"+construct_name+"_diagram.pdf", "PDF")

def construct_on_templates(datadir, construct_names, templates, construct_is_file=False):
	from os.path import isfile, splitext

	construct_files=[]
	truncated_construct_names=[]
	for construct in construct_names:
		actual_name, ext = splitext(construct)
		if ext in [".fasta", ".FASTA", ".gbk", ".GBK", ".genbank"]:
			construct_files += [datadir+construct]
			truncated_construct_names += [actual_name]
		else:
			construct_files += [datadir+construct+".gbk"]
		if type(construct_files) is str:
			construct_files = list(construct_files)
	construct_names = truncated_construct_names

	template_files = [datadir+"{0}.gbk".format(i) for i in templates]
	files_list = construct_files+template_files
	names_list = construct_names+templates

	for entry_name, entry_file in zip(names_list, files_list):
		_, ext = splitext(entry_name)
		if ext in [".fasta", ".FASTA", ".gbk", ".GBK", ".genbank"]:
			pass
		elif not isfile(entry_file):
			from write import write_seq
			write_seq(sequence_write_path=datadir, entrez_id=entry_name, formats=["genbank", "fasta"])

	gdd = GenomeDiagram.Diagram(str(construct_names)+' Hits on Templates', x=0.05, track_size=0.2)

	i = 0
	max_len=0
	for template_file, template in zip(template_files, templates):
		template_record = SeqIO.read(template_file, "gb")
		max_len = max(max_len, len(template_record))
		template_track, template_features = new_track(gdd, ", ".join(map(str, construct_names))+" on "+template, track_no=1, endpoint=len(template_record))
		for feature in template_record.features:
			if feature.type != "gene" and feature.type != "regulatory":
				# Exclude nongene features
				continue
			if len(template_features) % 2 == 0:
				color = colors.thistle
			else:
				color = colors.lightblue
			template_features.add_feature(feature, sigil="ARROW", color=color, label_color=color, label=True, label_size = 14, label_angle=0, arrowshaft_height=1)
		for construct_name, construct_file in zip(construct_names, construct_files):
			add_to_track(template_features, construct_file, template_file, annotation=construct_name, feature_color=colors.deeppink, sigil="BOX", label_angle=60)
		i += 1

	gdd.draw(format='linear', pagesize="A4", fragments=1, start=0, end=max_len)
	gdd.write("output/"+"_".join(map(str, construct_names))+"_template_hits.pdf", "PDF")

if __name__ == '__main__':
	restriction_list=[
		# ("GAATTC","EcoRI",colors.chartreuse),
		# ("CCCGGG","SmaI",colors.orange),
		# ("GTTAAC","HindI",colors.magenta),
		# ("GTCGAC","HindI",colors.purple),
		# ("GAATC","HinfI",colors.fuchsia),
		# ("GATTC","HinfI",colors.fuchsia),
		# ("GAGTC","HinfI",colors.fuchsia),
		# ("GACTC","HinfI",colors.fuchsia),
		# ("AGCGCT","HaeII",colors.cyan),
		# ("GGCGCC","HaeII",colors.cyan),
		("CGGCCG","EagI",colors.chartreuse),
		("GGCGCGCC","AscI",colors.deeppink),
		("GGATCC","BamHI",colors.turquoise)
		]
	cre(datadir="/home/chymera/data/CreBLseq/Cre_aj627603/", construct_name="aj627603")
	# my_hit(datadir="/home/chymera/data/CreBLseq/Cre_aj627603/", construct_name="my_full_hit")
	construct_on_templates(datadir="/home/chymera/data/CreBLseq/Cre_aj627603/", construct_names=["cre_f1.fasta", "cre_r1.fasta", "cre_fw2.fasta"], templates=["AF298785", "AF298789", "aj627603"], construct_is_file=True)
