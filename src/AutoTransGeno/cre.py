__author__ = "Horea Christian"

def cre(datadir, feature_name, extract_from):
	from utils import extract_feature, check_fetch_record
	from restriction import enzyme_selector
	from Bio import SeqIO
	from Bio.Graphics import GenomeDiagram
	from draw import new_track, draw_digest, add_to_track
	from reportlab.lib import colors
	from Bio import Restriction

	if feature_name is not str:
		construct_name = feature_name[0]
	else:
		construct_name = feature_name

	main_record, main_record_file = extract_feature(sequence_id=extract_from, data_dir=datadir, feature_names=feature_name, write_file=True)
	gdd = GenomeDiagram.Diagram(construct_name+' Construct Diagram', x=0.05, track_size=0.2)

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

	primer_track, primer_features = new_track(gdd, construct_name+" primers", smalltick=10)
	add_to_track(primer_features, datadir+"primers/"+"cre_fw2.fasta", main_record_file, annotation="cre_fw2", feature_color=colors.darkorchid, label_angle=60)
	add_to_track(primer_features, datadir+"primers/"+"cre_rv2.fasta", main_record_file, annotation="cre_rv2", feature_color=colors.darkorchid, label_angle=60)
	add_to_track(primer_features, datadir+"primers/"+"cre_fw3.fasta", main_record_file, annotation="cre_fw3", feature_color=colors.salmon, label_angle=60)
	add_to_track(primer_features, datadir+"primers/"+"cre_rv3.fasta", main_record_file, annotation="cre_rv3", feature_color=colors.salmon, label_angle=60)
	add_to_track(primer_features, datadir+"primers/"+"cre_f1.fasta", main_record_file, annotation="cre_f1", feature_color=colors.cornflower, label_angle=60)
	add_to_track(primer_features, datadir+"primers/"+"cre_r1.fasta", main_record_file, annotation="cre_r1", feature_color=colors.cornflower, label_angle=60)
	add_to_track(primer_features, datadir+"primers/"+"cre_fw4.fasta", main_record_file, annotation="cre_fw4", feature_color=colors.green, label_angle=60)
	add_to_track(primer_features, datadir+"primers/"+"cre_rv4.fasta", main_record_file, annotation="cre_rv4", feature_color=colors.green, label_angle=60)
	add_to_track(primer_features, datadir+"primers/"+"cre_fw5.fasta", main_record_file, annotation="cre_fw5", feature_color=colors.green, label_angle=60)
	add_to_track(primer_features, datadir+"primers/"+"cre_rv5.fasta", main_record_file, annotation="cre_rv5", feature_color=colors.green, label_angle=60)
	add_to_track(primer_features, datadir+"primers/"+"cre_fw6.fasta", main_record_file, annotation="cre_fw6", feature_color=colors.green, label_angle=60)
	add_to_track(primer_features, datadir+"primers/"+"cre_rv6.fasta", main_record_file, annotation="cre_rv6", feature_color=colors.green, label_angle=60)
	add_to_track(primer_features, datadir+"primers/"+"cre_fw7.fasta", main_record_file, annotation="cre_fw7", feature_color=colors.green, label_angle=60)
	add_to_track(primer_features, datadir+"primers/"+"cre_rv7.fasta", main_record_file, annotation="cre_rv7", feature_color=colors.green, label_angle=60)

	gdd.draw(format="linear", pagesize="A4", fragments=1, start=0, end=len(main_record))
	gdd.write("/home/chymera/src/AutoTransGeno/output/"+construct_name+"_from_"+extract_from+".pdf", "PDF")

	return restriction_dict
if __name__ == '__main__':
	cre(datadir="/home/chymera/data/reference/sequences/", feature_name=["cre", "Cre", "CRE"], extract_from="aj627603")
