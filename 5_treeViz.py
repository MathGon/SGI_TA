import os
from ete3 import Tree, NodeStyle, TreeStyle, TextFace, AttrFace
from Bio import SeqIO
from Bio import Entrez

def taxonToColour(gi):
	""" a partir d'un numero de gi, determine le taxon (dans le fichier gi_taxon.csv) et ainsi la couleur a apporter """
	input_handle = open("gi_taxon.csv", "rU")
	for line in input_handle:
		if gi == line.split("\t")[0]:
			tax = line.split(";")[3].lstrip()
			if tax == "Alphaproteobacteria": #proteo
				return [tax,"lightcoral"]
			if tax == "Gammaproteobacteria": #proteo
				return [tax,"indianred"]
			if tax == "Betaproteobacteria": #proteo
				return [tax,"firebrick"]
			if tax == "delta/epsilon subdivisions": #proteo
				return [tax,"salmon"]
				
			if tax == "Bacteroidetes": #Bacteroidetes
				return [tax,"darkkhaki"]
			if tax == "Bacteroidia" : #Bacteroidetes
				return [tax,"olive"]	
			if tax == "Ignavibacteriae": #Bacteroidetes
				return [tax,"darkorange"]
			
			if tax == "Clostridia": #Firmicutes
				return [tax,"royalblue"]
			if tax == "Bacilli": #Firmicutes
				return [tax,"darkblue"]
				
			if tax == "Oscillatoriophycideae": #Cyanobacteria
				return [tax,"limegreen"]
			if tax == "Nostocales": #Cyanobacteria
				return [tax,"green"]
						
			if tax == "Actinobacteria": #Actinobacteria
				return [tax,"darkviolet"]

			if tax == "Fimbriimonadia": #Armatimonadetes
				return [tax,"darkslategrey"]

			if tax == "Deferribacteres": #Deferribacteres
				return [tax,"plum"]

			if tax == "Spirochaetia": # Spirochaetes
				return [tax,"grey"]


def drawTree(treeFile, ShowBool):
	"""
	Draw a tree from a phy file
	"""
	t = Tree(treeFile)
	imgFile = treeFile.replace(".tree", ".tree.pdf")

	# Basic tree style
	ts = TreeStyle()
	ts.show_leaf_name = True
	ts.show_branch_support = True
	ts.scale =  50

	# Draws nodes as small red spheres of diameter equal to 10 pixels
	nstyle = NodeStyle()
	nstyle["shape"] = "sphere"
	nstyle["size"] = 10
	nstyle["fgcolor"] = "darkred"
	#nstyle["faces_bgcolor"] = "pink"

	nstyle2 = NodeStyle()
	nstyle2["shape"] = "sphere"
	nstyle2["size"] = 10
	nstyle2["fgcolor"] = "darkblue"
	
	# Gray dashed branch lines
	nstyle["hz_line_type"] = 1
	nstyle["hz_line_color"] = "#cccccc"

	# Applies the same static style to all nodes in the tree. Note that,
	# if "nstyle" is modified, changes will affect to all nodes
	for n in t.traverse():
		if n.is_leaf():
			if n.name.split("|")[-1] == "GI":
				n.set_style(nstyle)
			if n.name.split("|")[-1] == "plasmid":
				n.set_style(nstyle2)
			gi = n.name.split("|")[1]
			n.name = n.name.split("|")[0] #+ "   " + n.name.split("|")[1]
			n.name = n.name.replace("_tRNA_modification_GTPase_", "")
			n.name = n.name.replace("_DNA", "")
			n.name = " " + n.name + " "
			if n.name[-1] == "_": n.name.rstrip()
			
			taxon, color = taxonToColour(gi)
			n.add_face(TextFace(taxon, fgcolor = color, fsize = 8), column=1, position="branch-right")
			#n.img_style["bgcolor"] = color
			
	if ShowBool == True: #permet de flipper les braches pour avoir des topologies similaires
		t.show(tree_style=ts)
	t.render(imgFile, w=1024, units="mm", tree_style=ts)
	
drawTree("S025S026_fusion.tree", True)
