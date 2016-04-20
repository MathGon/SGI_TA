import os
from ete3 import Tree, NodeStyle, TreeStyle, TextFace, AttrFace
from Bio import SeqIO
from Bio import Entrez

""" Alignement et phylogenie des homologues de S025 et S026 """

def taxonomyGI(gi):
	""" retourne la taxonomy lineage d'un genome par appel de son numero GI """
	Entrez.email = "dfog22@hotmail.com"
	handle = Entrez.efetch(db="nucleotide", rettype="gbwithparts", retmode="text", id=gi)
	seq_record = SeqIO.read(handle, "gb") #using "gb" as an alias for "genbank"
	windowSize = 30000 # regarder a plus ou moins windoSize
	ilotRelatedWords = ["integrase", "recombinase", "tranposase", "metal", "lactamase", "resistance", "antibiotic", "addiction", "trme", "thdf", "mnme", "mercury", "resistance", "resolvase", "phage"]
	
	for feature in seq_record.features:
		#print feature
		if feature.type == "source":
			if "db_xref" in feature.qualifiers:
				try:
					taxonid = feature.qualifiers["db_xref"][0].split(":")[1]
					handle = Entrez.efetch(db="Taxonomy", id=taxonid, retmode="xml")
					records = Entrez.read(handle)
					return records[0]["Lineage"]
				except:
					return "error"
# print taxonomyGI("74419069")

def GIlist(fileName):
	""" A partir d'une liste de sequence fasta contenant le numero de gi dans le header,
		creation d'un fichier .csv contenant les taxons
	"""
	input_handle = open(fileName, "rU")
	records = SeqIO.parse(input_handle, "fasta")
	#GIlist = []
	output_handle = open("gi_taxon.csv", "w")
	for record in records:
		if len(record.id.split("|"))> 1:
			gi = record.id.split("|")[1]
			taxonLineage = taxonomyGI(gi)
			toWrite = gi + "\t" + taxonLineage + "\n"
			output_handle.write(toWrite)
	return GIlist
#GIlist("S025_homo_commun_annotated.fas")	
	
def filterSequences(fileName, giList):
	""" Supprime des deux fichiers d'homologues les sequences dont les gi correspondent a ceux fournit dans la liste """
	input_handle = open(fileName, "rU")
	outfileName = fileName.replace(".aln", "_filtred.aln")
	output_handle = open(outfileName, "w")
	records = SeqIO.parse(input_handle, "fasta")
	for record in records:
		if len(record.id.split("|"))> 1:
			gi = record.id.split("|")[1]
			if gi not in giList:
				 SeqIO.write(record, output_handle, "fasta")
	input_handle.close()
	output_handle.close()

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
print taxonToColour("292673316")
	
def drawTree(treeFile, ShowBool):
	"""
	Draw a tree from a phy file
	"""
	t = Tree(treeFile)
	imgFile = treeFile.replace(".tree", ".tree.png")

	# Basic tree style
	ts = TreeStyle()
	ts.show_leaf_name = True
	ts.show_branch_support = True
	ts.scale =  160

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
	t.render(imgFile, w=393, units="mm", tree_style=ts)

"""
# Alignement
os.system("muscle -in S025_homo_commun_annotated.fas -out S025_homo_commun_annotated.aln")
os.system("muscle -in S026_homo_commun_annotated.fas -out S026_homo_commun_annotated.aln")
"""

# Filtrage des sequences, numero de gi dont les branches sont mal supportees ou seq tronquees
giAfiltrer = ["152026452", "339305362", "311692891", "296012614", "442770505", "255292892","255292743","442770505", "", "255293164", "32446539", "148498119", "675818951", "332337265", "74055513", "569540043", "834771869"]
filterSequences("S025_homo_commun_annotated.aln", giAfiltrer)
filterSequences("S026_homo_commun_annotated.aln", giAfiltrer)

# Phylogenie
os.system("FastTree S025_homo_commun_annotated_filtred.aln > S025_homo_commun_annotated.phy.tree")
os.system("FastTree S026_homo_commun_annotated_filtred.aln > S026_homo_commun_annotated.phy.tree")


drawTree("S025_homo_commun_annotated.phy.tree", True)
drawTree("S026_homo_commun_annotated.phy.tree", True)

