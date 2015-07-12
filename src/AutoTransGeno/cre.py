__author__ = "Horea Christian"

from os.path import splitext, basename, expanduser
from utils import extract_feature, check_fetch_record, write_seq, simple_sequence_merge, convert_seq, concatenate_overlapping_sequences
from restriction import enzyme_selector
from Bio import SeqIO, SeqFeature
from Bio.Graphics import GenomeDiagram
from draw import new_track, draw_digest, add_to_track
from reportlab.lib import colors
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
from Bio.Alphabet.IUPAC import IUPACAmbiguousDNA
from itertools import cycle


def cre(project_data_dir="", global_data_dir="", feature_name="", construct_name="CRE", extract_from="", main_record_file="", align_with=[], primer_list=[], restriction_interval=[], rb=[], write=False, pagesize="A4", scale_fontsize=3, label_size=2, greytrack_fontsize=7, x=0.05, y=0.01, track_size=0.3):

	#define color palettes:
	primer_colors=[colors.orchid, colors.cornflower, colors.lightseagreen, colors.salmon]

	if extract_from:
		main_record, main_record_file = extract_feature(sequence_id=extract_from, data_dir=global_data_dir, feature_names=feature_name, write_file=True)
	elif main_record_file:
		extract_from = splitext(basename(main_record_file))[0]
		if splitext(basename(main_record_file))[1] == ".fasta":
			main_record = SeqIO.read(main_record_file, "fasta")
		elif splitext(basename(main_record_file))[1] in [".gb", ".gbk", ".genbank"]:
			main_record = SeqIO.read(main_record_file, "gb")
			main_record_file = convert_seq(main_record_file, "genbank", "fasta")
			print main_record_file

	gdd = GenomeDiagram.Diagram(construct_name+' Construct Diagram', x=x, y=y, track_size=track_size)

	genbank_track, genbank_features = new_track(gdd, construct_name+" features", smalltick=10, scale_fontsize=scale_fontsize, greytrack_fontsize=greytrack_fontsize)
	for feature in main_record.features:
		if "Cre" in str(feature.qualifiers):
			color = colors.lavender
		else:
			color = colors.grey
		genbank_features.add_feature(feature, sigil="ARROW", color=color, label_color=color, label=True, label_size=label_size, label_angle=30, arrowshaft_height=1)

	if restriction_interval or rb:
		restriction_dict = enzyme_selector(sequence=main_record, restriction_interval=[0,690], genome_frequency=[700,2000], deterministic_overhangs=True, rb=["XceI", "PsuI"])
		restriction_track, restriction_features = new_track(gdd, construct_name+" restriction sites", smalltick=10, scale_fontsize=scale_fontsize, greytrack_fontsize=greytrack_fontsize)
		draw_digest(restriction_features, restriction_dict)

	# plotting primers
	if primer_list:
		primer_colors = cycle(primer_colors)
		primer_track, primer_features = new_track(gdd, construct_name+" primers", smalltick=10, scale_fontsize=scale_fontsize, greytrack_fontsize=greytrack_fontsize)
		for primer_entry in primer_list:
			primer_color = primer_colors.next()
			add_to_track(primer_features, global_data_dir+"primers/"+primer_entry[0]+".fasta", main_record_file, annotation=primer_entry[0], feature_color=primer_color, label_angle=30, label_size=label_size)
			add_to_track(primer_features, global_data_dir+"primers/"+primer_entry[1]+".fasta", main_record_file, annotation=primer_entry[1], feature_color=primer_color, label_angle=30, label_size=label_size)

	# turn entry names into actual file paths
	align_with=[project_data_dir+entry+".fasta" for entry in align_with]
	for i in align_with:
		hit_track_back, hit_features_back = new_track(gdd, splitext(i)[0][-6:], smalltick=10, end=len(SeqIO.read(i, 'fasta')), scale_fontsize=scale_fontsize, greytrack_fontsize=greytrack_fontsize)
		add_to_track(hit_features_back, main_record_file, i, annotation=" "+construct_name, feature_color=colors.red, label_angle=30, forceone=True, label_size=label_size)
		hit_track, hit_features = new_track(gdd, construct_name+" alignment hits", smalltick=10, scale_fontsize=scale_fontsize, greytrack_fontsize=greytrack_fontsize)
		hsp_list = add_to_track(hit_features, i, main_record_file, annotation=" "+splitext(i)[0][-6:], feature_color=colors.red, label_angle=30, forceone=True, label_size=label_size)

		record = SeqIO.read(i, "fasta")
		#for loop only takes the first alignment, then breaks
		for hsp in hsp_list:
			truncated_record = SeqRecord(record.seq[hsp.query_end:], id = "Region downstream of Cre in ePet-cre mice, i="+i)
			write_seq(sequence_write_path=project_data_dir, record=truncated_record, ID=splitext(basename(i))[0]+"_3-unmatched")
			break

	if align_with:
		record = SeqRecord(main_record.seq+record.seq[hsp.query_end:], id = "Cre and following bases in ePet-cre construct.")
		write_seq(sequence_write_path=project_data_dir, record=record, ID="cre-ff-current")

	gdd.draw(format="linear", pagesize=pagesize, fragments=1, start=0, end=len(main_record))
	if write:
		gdd.write("/home/chymera/src/AutoTransGeno/output/"+construct_name+"_from_"+extract_from+".pdf", "PDF")
		print "/home/chymera/src/AutoTransGeno/output/"+construct_name+"_from_"+extract_from+".pdf"
	return gdd

def write_gb(main_record_file, add="", destination=""):
	destination = expanduser(destination)
	main_record_file = expanduser(main_record_file)
	main_record = SeqIO.read(main_record_file, "fasta")
	add = expanduser(add)
	add = SeqIO.read(add, "fasta")

	main_record.seq.alphabet = IUPACAmbiguousDNA()
	cre_end = len(main_record.seq)
	main_record.seq = Seq(concatenate_overlapping_sequences(main_record.seq.tostring(), add.seq.tostring()), IUPACAmbiguousDNA())
	main_record.name = "ePet-cre"
	main_record.id = "ePet-cre"
	main_record.description = "ePet-cre construct from doi:10.1038/nn.2623"

	my_start_pos = SeqFeature.ExactPosition(0)
	my_end_pos = SeqFeature.ExactPosition(cre_end)

	my_feature_location = SeqFeature.FeatureLocation(my_start_pos,my_end_pos, strand=1)
	my_feature_type = "CDS"
	my_feature = SeqFeature.SeqFeature(my_feature_location,type=my_feature_type)
	my_feature.qualifiers["gene"] = "Cre"
	main_record.features.append(my_feature)


	# Add SV40 PolyA
	my_start_pos = SeqFeature.ExactPosition(cre_end)
	my_end_pos = SeqFeature.ExactPosition(cre_end+118)

	my_feature_location = SeqFeature.FeatureLocation(my_start_pos,my_end_pos, strand=1)
	my_feature_type = "PolyA"
	my_feature = SeqFeature.SeqFeature(my_feature_location,type=my_feature_type)
	my_feature.qualifiers["PolyA"] = "SV40-PolyA"
	main_record.features.append(my_feature)

	SeqIO.write(main_record, destination+".gb", "genbank")
	SeqIO.write(main_record, destination+".fastas", "fasta")

def meld_sequences(ID_list = ["783476", "old_783477", "783477"], base_dir = "/home/chymera/", save_sequences_to=""):
	# creating sequencing meld
	seq_paths = [base_dir + ID + ".fasta" for ID in ID_list]
	sequences = [SeqIO.read(one_path, "fasta") for one_path in seq_paths]
	outp = simple_sequence_merge(sequences=sequences)
	seq_meld = write_seq(sequence=outp,sequence_write_path=save_sequences_to, ID="meld_"+"+".join(ID_list))
	return seq_meld

# def extend_sequence(base_sequence, new_sequence, new_features):


if __name__ == '__main__':
	primer_list=[
			["cre_f1", "cre_r1", colors.orchid],
			["cre_fw3", "cre_rv3", colors.cornflower],
			["cre_fw5", "cre_rv5", colors.lightseagreen],
			["cre_fw7", "cre_rv7", colors.salmon]
		]
	cre(main_record_file="/home/chymera/data/gt.ep/genbank/cre-ff-current.gb", project_data_dir="/home/chymera/data/gt.ep/fasta/", global_data_dir="/home/chymera/data/reference/sequences/", align_with=["meld_783476+old_783477+783477"], primer_list=primer_list, write=True, scale_fontsize=3, label_size=5)

	# write_gb("~/data/reference/sequences/cre_from_aj627603.fasta", add="~/data/gt.ep/fasta/meld_783476+old_783477+783477.fasta", destination="~/data/gt.ep/genbank/cre-ff-current")
	# write_gb("~/data/gt.ep/fasta/cre-ff-current.fasta", add="~/data/gt.ep/fasta/meld_783476+old_783477+783477.fasta", destination="~/data/gt.ep/genbank/cre-ff-current.gb")
