#!/usr/bin/env python
__author__ = "Horea Christian"

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
	# my_hit(datadir="/home/chymera/data/CreBLseq/Cre_aj627603/", construct_name="my_full_hit")
	construct_on_templates(datadir="/home/chymera/data/CreBLseq/Cre_aj627603/", construct_names=["cre_f1.fasta", "cre_r1.fasta", "cre_fw2.fasta"], templates=["AF298785", "AF298789", "aj627603"], construct_is_file=True)
