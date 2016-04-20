from Bio import SeqIO
from Bio import Entrez
"""
A partir de la liste des sequences utilises, le but est d'annoter la region dans laquelle ce trouve ce gene: ilot ou pas?
On part de la liste des sequences fasta selectionnees, puis on remonte dans le fichier xxx pour trouver les coordonnees du
gene d'interet et on regarde les annotations a plus ou moins 30kb a la recherche de trace suspectant un ilot
"""

def annotateRegion(gi, start):
	"""
	Regarde sur le genome gi a position start, plus ou moins windowSize bp, 
	si presence de gene lies a ilot
	: gi arg : numero d'accesion ncbi
	: gi type : char
	: ret : "GI" ou "CHR"
	"""
	Entrez.email = "dfog22@hotmail.com"
	handle = Entrez.efetch(db="nucleotide", rettype="gbwithparts", retmode="text", id=gi)
	seq_record = SeqIO.read(handle, "gb") #using "gb" as an alias for "genbank"
	windowSize = 30000 # regarder a plus ou moins windoSize
	ilotRelatedWords = ["integrase", "recombinase", "tranposase", "metal", "lactamase", "resistance", "antibiotic", "addiction", "trme", "thdf", "mnme", "mercury", "resistance", "resolvase", "phage"]
	
	for feature in seq_record.features:
		#print feature
		#if feature.type == "source":
		#	strain = feature.qualifiers["organism"][0].replace(" ", "_")
		if feature.location.start >= start - windowSize and feature.location.start >= start + windowSize:
			if feature.type == "CDS":
				if "product" in feature.qualifiers:
					words = feature.qualifiers["product"][0].split(" ")
					for word in words:
						if word in ilotRelatedWords:
							return "GI"
	return "chr"

def annotationDuGI(gi):
	with open('gi_annotation.csv', 'r') as f:
		lignes = f.readlines()
		for ligne in lignes:
			giligne = ligne.split("\t")[0]
			annotation = ligne.split("\t")[1]
			if giligne == gi:
				print annotation.rstrip()
				return annotation.rstrip()
#annotationDuGI("615567765")

def annotationDuFasta(fastaFileName):
	""" transfert l'annotation depuis annotation.txt dans le header d'un fichier fasta """
	fastaFileIn = open(fastaFileName, "r") 
	outfileName = fastaFileName.replace(".fas", "_annotated.fas")
	fastaFileAnnotated = open(outfileName, "w")
	records = SeqIO.parse(fastaFileIn, "fasta")
	for record in records:
		if len(record.id.split("|"))> 1:
			gi = record.id.split("|")[1]
			annotation = record.id.split("|")[-1]
			if annotation == "unknown":
				annotation = annotationDuGI(gi)
				name = record.id.replace("unknown", annotation)
				toWrite = str(">" + name + "\n" + record.seq +"\n")
				fastaFileAnnotated.write(toWrite)
			else:
				toWrite = str(">" + record.id + "\n" + record.seq +"\n")
				fastaFileAnnotated.write(toWrite)
	fastaFileAnnotated.close()
	
# Genere le fichier annotation.txt contenant numero gi et son assignation
fichierAnnot = open("gi_annotation.csv", "w")
fichierFas = open("S025_homo_commun.fas") #contient les prot selectionnees pour arbre
records = SeqIO.parse(fichierFas, "fasta")
forcage_gi = ["363406024", "126636230", "803440056", "667714421", "268309235", "302567367"] #suite a analyse, ces genomes ne sont pas bien annotes, il faut corriger a la main
for record in records:
	if record.id.split("|")[-1] == "unknown":
		gi = record.id.split("|")[1]
		if gi in forcage_gi:
			toWrite = gi + "\t" + "GI" + "\n"
			fichierAnnot.write(toWrite)
		else:
			pos = record.id.split("|")[2]
			start = int(pos.split("-")[1])
			toWrite = gi + "\t" + annotateRegion(gi, start) + "\n"
			fichierAnnot.write(toWrite)
fichierAnnot.close()
fichierFas.close()

#transfert annotation dans le header fasta des fichier d'homolog
annotationDuFasta("S025_homo_commun.fas")
annotationDuFasta("S026_homo_commun.fas")
