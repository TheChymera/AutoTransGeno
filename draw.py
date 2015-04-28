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
	from write import check_format

	subject = check_format(subject, "fasta")

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

def draw_digest(track_features, restriction_dict, restriction_colors="", rotate_labels=False):
	from itertools import cycle
	from sequence_utils import overhangs

	#the labels are too close to the track to properly display multiple rotated instances pointing at the same position, this gets prepended:
	name_spacer="      "

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
	if rotate_labels:
		angle = cycle([120,90,60,30])
	for enzyme, sites in restriction_dict.iteritems():
		current_color = restriction_colors.next()
		if rotate_labels:
			current_angle=angle.next()
		else:
			current_angle = 90
		overhang_length = len(overhangs(enzyme))
		for site in sites:
			feature = SeqFeature(FeatureLocation(site, site+overhang_length))
			track_features.add_feature(feature, color=current_color, name=name_spacer+enzyme.__name__, label=True, label_size=3, label_position="end", label_color=current_color, label_angle=current_angle)
