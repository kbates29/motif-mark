import argparse
import cairo
import math
import re
import numpy as np
import random


def get_arguments():
	parser = argparse.ArgumentParser(description= "Arguments the user wishes to pass")
	parser.add_argument( "-f", "--file",
	help="The sequence in fasta format containing both exons and introns for a given gene. Exons sequences should be capitalized.", 
	required=True, type=str)
	parser.add_argument("-m", "--motif",
	help="The file which contains the motif of interest. Each motif should be lowercase and on a newline.", 
	required=True, type=str)
	return parser.parse_args()
	
args = get_arguments()
file = args.file
motif = args.motif


# searches through sequence for motifs and returns start and stop position of given motif in motif_list
def search_motif(seq, motif_list):
	seq = seq.lower()
	motif_dict = {}
	for motif in motif_list:
		lower_motif = motif
		motif_dict[motif] = []
		motif = re.sub("u", "[t]", motif)
		motif = re.sub("y", "[ct]", motif)
		motif = re.sub("w", "[at]", motif)
		motif = re.sub("s", "[cg]", motif)
		motif = re.sub("m", "[ac]", motif)
		motif = re.sub("k", "[gt]", motif)
		motif = re.sub("r", "[ag]", motif)
		motif = re.sub("b", "[cgt]", motif)
		motif = re.sub("d", "[agt]", motif)
		motif = re.sub("h", "[act]", motif)
		motif = re.sub("v", "[acg]", motif)
		motif = re.sub("n", "[acgt]", motif)

		for match in re.finditer(motif, seq):
			motif_dict[lower_motif].append(match.span())
		
	return(motif_dict)

# searches for exon start and stop position sequence and returns list of spans
def search_exon(seq):
	exon_pos = []
	for match in re.finditer("[A-Z]+", seq):
		exon_pos.append(match.span())
	return(exon_pos)


# populates list with same length as motif list which corresponds to rgb color values 
def colors(num_motif):
	color_list = []
	for i in range(num_motif):
		r = random.random()
		g = random.random()
		b = random.random()
		color_list.append((r,g,b))

	return color_list

#draws gene with introns as lines, exons as black box and motif in colored box all proportaional to size of element
def draw_motifs(gene_dict, exon_dict, position_dict, motif_list, longest_seq, all_colors):
	width = longest_seq + 225
	height = (len(gene_dict) * 110) 
	surface = cairo.SVGSurface("Motifs.svg", width, height)
	context = cairo.Context(surface)
	header_y = 50
	intron_y = 75
	exon_y = 69
	motif_y =69 
	legend_x_text = longest_seq + 80
	legend_y_text = 55
	legend_col_x = longest_seq + 70
	legend_col_y = 47
	color_indexer = 0 
	drawn_legend = False
	context.select_font_face("Arial", cairo.FONT_SLANT_NORMAL, cairo.FONT_WEIGHT_NORMAL)


	for key in gene_dict:
		# writes in header above line
		context.set_font_size(20)
		context.set_source_rgb(0, 0, 0)
		context.set_line_width(5)
		context.move_to(50, header_y)
		context.show_text(key)
		header_y += 100
		
		# add introns
		context.move_to(50, intron_y)
		context.line_to(len(gene_dict[key]) + 50, intron_y)
		context.stroke()
		intron_y += 100


		# add exon
		for pair in exon_dict[key]: # pair is start and stop postion of exon
			context.rectangle(pair[0] + 50, exon_y, (pair[1]-pair[0]), 11)
			context.fill()
			exon_y += 100
	

		# add motifs and puts in legend
		color_indexer = 0
		for pair in position_dict[key]:
			context.set_source_rgba(all_colors[color_indexer][0], all_colors[color_indexer][1], all_colors[color_indexer][2], 1.0)
			color_indexer +=1
			for span in position_dict[key][pair]:
				context.rectangle(span[0] + 50, motif_y, (span[1]-span[0]), 11)
				context.fill()
		motif_y += 100
		
		if drawn_legend == False:
			color_indexer = 0 
			for motif in orig_motif: 
				context.set_source_rgb(0, 0, 0)
				context.set_line_width(1)
				context.rectangle(longest_seq + 60, 40, 125, (len(orig_motif) + 30) * len(orig_motif))
				context.stroke()
				context.set_font_size(12)
				context.move_to(legend_x_text, legend_y_text)
				context.show_text(motif)
				context.set_source_rgba(all_colors[color_indexer][0], all_colors[color_indexer][1], all_colors[color_indexer][2], 1.0)
				context.rectangle(legend_col_x, legend_col_y, 8, 11)
				context.fill()
				legend_y_text += 35
				legend_col_y += 35
				color_indexer +=1
		drawn_legend = True
	
	surface.finish()

# populates a list with all motifs specified in motif input file
with open(motif, "r") as m:
	motif_list = []
	num_motif = 0
	orig_motif = []
	for line in m:
		line = line.strip()
		orig_motif.append(line)
		line = line.lower()
		motif_list.append(line)
		num_motif += 1


#populates dictionary with header line as key and sequence as value
with open(file, "r") as f:
	gene_dict = {}
	longest_seq = 0
	for line in f:
		line = line.strip()
		if line.startswith(">"):
			gene = line[1:]
			gene_dict[gene] = ""
		else:
			gene_dict[gene] += line
			if len(gene_dict[gene]) > longest_seq:
				longest_seq = len(gene_dict[gene])


position_dict = {}
exon_dict = {}
for key in gene_dict:
	position_dict[key] = search_motif(gene_dict[key], motif_list)
	exon_dict[key] = search_exon(gene_dict[key])

all_colors = colors(num_motif)
draw_motifs(gene_dict, exon_dict, position_dict, orig_motif, longest_seq, all_colors)