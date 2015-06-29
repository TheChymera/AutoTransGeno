__author__ = "Horea Christian"

from os.path import splitext, basename
from utils import extract_feature, check_fetch_record, write_seq, simple_sequence_merge
from restriction import enzyme_selector
from Bio import SeqIO
from Bio.Graphics import GenomeDiagram
from draw import new_track, draw_digest, add_to_track
from reportlab.lib import colors
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
from itertools import cycle
from Bio.SeqFeature import SeqFeature, FeatureLocation
import numpy as np

def mockup(features_list, write=False, pagesize="A4", scale_fontsize=3, label_size=2, greytrack_fontsize=7, x=0.05, y=0.01, track_size=0.2, track_names="", scale_ticks=False, format="linear", total_len=3000):
	colors_cycle=[colors.orchid, colors.cornflower, colors.lightseagreen, colors.cornflower, colors.salmon]
	colors_cycle=cycle(colors_cycle)

	gdd = GenomeDiagram.Diagram('Construct Diagram', x=x, y=y, track_size=track_size)

	for ix, track_info in enumerate(features_list):
		track_len=0
		for i in track_info:
			track_len += i[1]
		track, features = new_track(gdd, " "+track_names[ix], smalltick=10, scale_fontsize=scale_fontsize, greytrack_fontsize=greytrack_fontsize, scale_ticks=scale_ticks, end=track_len)
		feature_start=0
		for feature_info in track_info:
			if feature_info[0] == "skip":
				feature_start += feature_info[1]
				continue
			feature = SeqFeature(FeatureLocation(feature_start,feature_start+feature_info[1]), strand=feature_info[2])
			if feature_info[0] == "LoxP":
				feature_color = colors.yellow
			elif feature_info[0] == "STOP":
				feature_color = colors.red
			elif feature_info[0] == "Restriction":
				feature_color = colors.chartreuse
			else:
				feature_color = colors_cycle.next()
			features.add_feature(feature,
				name=feature_info[0],
				label=True,
				color=feature_color,
				label_size = label_size,
				label_color=feature_color,
				label_angle=30,
				sigil=feature_info[3],
				arrowshaft_height=1)
			feature_start += feature_info[1]

	gdd.draw(format=format, pagesize=pagesize, fragments=1, start=0, end=total_len)
	if write:
		gdd.write("/home/chymera/src/AutoTransGeno/output/test.pdf", "PDF")

	return gdd

if __name__ == '__main__':
	features_list=[
			[
				["ePet", 2000, +1],
				["cre", 1100, +1]
			],
			[
				["promoter", 100, +1],
				["LoxP", 34, +1],
				["LoxP", 34, +1]
			],
		]
	mockup(features_list, write=True, pagesize="A4", scale_fontsize=3, label_size=2, greytrack_fontsize=7, x=0.05, y=0.01, track_size=0.3, track_names=["On Genme","On Vctor"])
