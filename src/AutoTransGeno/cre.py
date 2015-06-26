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

primer_colors=[colors.orchid, colors.cornflower, colors.lightseagreen, colors.salmon]

def cre(data_dir, feature_name, extract_from, save_sequences_to="", align_with=[]):

	if feature_name is not str:
		construct_name = feature_name[0]
	else:
		construct_name = feature_name

	main_record, main_record_file = extract_feature(sequence_id=extract_from, data_dir=data_dir, feature_names=feature_name, write_file=True)

	gdd = GenomeDiagram.Diagram(construct_name+' Construct Diagram', x=0.05, track_size=0.25)

	genbank_track, genbank_features = new_track(gdd, construct_name+" features", smalltick=10)
	for feature in main_record.features:
		if feature.type != "gene" and feature.type != "regulatory":
			# Exclude nongene features
			continue
		if "cre" in str(feature.qualifiers):
			color = colors.lavender
		else:
			color = colors.grey
		genbank_features.add_feature(feature, sigil="ARROW", color=color, label_color=color, label=True, label_size = 14, label_angle=0, arrowshaft_height=1)

	restriction_dict = enzyme_selector(sequence=main_record, restriction_interval=[0,690], genome_frequency=[700,2000], deterministic_overhangs=True, rb=["XceI", "PsuI"])
	restriction_track, restriction_features = new_track(gdd, construct_name+" restriction sites", smalltick=10)
	draw_digest(restriction_features, restriction_dict)

	primer_list=[
			["cre_f1", "cre_r1", colors.orchid],
			["cre_fw3", "cre_rv3", colors.cornflower],
			["cre_fw5", "cre_rv5", colors.lightseagreen],
			["cre_fw7", "cre_rv7", colors.salmon]
		]

	# plotting primers
	primer_track, primer_features = new_track(gdd, construct_name+" primers", smalltick=10)
	for primer_entry in primer_list:
		add_to_track(primer_features, data_dir+"primers/"+primer_entry[0]+".fasta", main_record_file, annotation=primer_entry[0], feature_color=primer_entry[2], label_angle=60)
		add_to_track(primer_features, data_dir+"primers/"+primer_entry[1]+".fasta", main_record_file, annotation=primer_entry[1], feature_color=primer_entry[2], label_angle=60)

	# turn entry names into actual file paths
	align_with=[save_sequences_to+entry+".fasta" for entry in align_with]
	for i in align_with:
		hit_track_back, hit_features_back = new_track(gdd, splitext(i)[0][-6:], smalltick=10, end=len(SeqIO.read(i, 'fasta')))
		add_to_track(hit_features_back, main_record_file, i, annotation=" "+construct_name, feature_color=colors.darkslateblue, label_angle=30, forceone=True)
		hit_track, hit_features = new_track(gdd, construct_name+" alignment hits", smalltick=10)
		hsp_list = add_to_track(hit_features, i, main_record_file, annotation=" "+splitext(i)[0][-6:], feature_color=colors.darksalmon, label_angle=30, forceone=True)

		record = SeqIO.read(i, "fasta")
		#for loop only takes the first alignment, then breaks
		for hsp in hsp_list:
			truncated_record = SeqRecord(record.seq[hsp.query_end:], id = "Region downstream of Cre in ePet-cre mice, i="+i)
			write_seq(sequence_write_path=save_sequences_to, record=truncated_record, ID=splitext(basename(i))[0]+"_3-unmatched")
			break

	record = SeqRecord(main_record.seq+record.seq[hsp.query_end:], id = "Cre and following bases in ePet-cre construct.")
	write_seq(sequence_write_path=save_sequences_to, record=record, ID="cre-ff-current")


	gdd.draw(format="linear", pagesize="A4", fragments=1, start=0, end=len(main_record))
	gdd.write("/home/chymera/src/AutoTransGeno/output/"+construct_name+"_from_"+extract_from+".pdf", "PDF")

	return restriction_dict

def meld_sequences(ID_list = ["783476", "old_783477", "783477"], base_dir = "/home/chymera/", save_sequences_to=""):
	# creating sequencing meld
	seq_paths = [base_dir + ID + ".fasta" for ID in ID_list]
	sequences = [SeqIO.read(one_path, "fasta") for one_path in seq_paths]
	outp = simple_sequence_merge(sequences=sequences)
	seq_meld = write_seq(sequence=outp,sequence_write_path=save_sequences_to, ID="meld_"+"+".join(ID_list))
	return seq_meld

if __name__ == '__main__':
	cre(data_dir="/home/chymera/data/reference/sequences/", feature_name=["cre"], extract_from="aj627603", save_sequences_to="/home/chymera/data/gt.ep/fasta/", align_with=["783476","783477"])
