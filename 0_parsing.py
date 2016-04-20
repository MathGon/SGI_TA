from Bio import SeqIO
from Bio.Blast import NCBIXML 
from Bio import Entrez

"""
Parse les resultats de tblastn.xml de S025 et S026
Extrait les matching subjct
"""

def strainNameHomogene(strain):
	""" Renomme le nom des souches pour condenser """
	strain = strain.split(",")[0]
	strain = strain.split("=")[0]
	strain = strain.split("(")[0]
	strain = strain.split("genome")[0]
	words = ["chromosome ", "chromosome", "genomic island ", "variant ", "main ", "complete ", "genomic ", "sequence ", "sequence", "island ", "Salmonella Genomic Island 1 ", "genome ", "strain ", "plasmid ", "str. ", "REU80928 "]
	for word in words:
		strain = strain.replace(word, "")
	strain = strain.replace(":", "-")
	strain = strain.replace(" ", "_")
	if strain[-1] == "_": strain.rstrip()
	if strain[-1] == "_": strain.rstrip()
	return strain
		

blastXMLfile = "S025_tblastn.xml"
result_handle = open(blastXMLfile, "r")
outAlnFile = "S025_homo.fas"
out_handle = open(outAlnFile, "w")
blast_records = NCBIXML.parse(result_handle)
compilResults = []
for record in blast_records:
	for alignment in record.alignments:
		loc = "unknown"
		strain = alignment.title.split("|")[-1].lstrip()
		if "plasmid" in strain:
			loc = "plasmid"
		if "island" in strain:
			loc = "GI"
		strain = strainNameHomogene(strain)
		gi = alignment.title.split("|")[1]
		for hsp in alignment.hsps:
			if hsp.query_start < 100 and hsp.query_end > 500:
				#print hsp
				header = str(">" + strain + "|" + gi + "|" + str(hsp.sbjct_start) + "-" + str(hsp.sbjct_end) +"|" + loc )
				toWrite = header + "\n" + hsp.sbjct + "\n"
				out_handle.write(toWrite)
				break
result_handle.close()
out_handle.close()
				
blastXMLfile = "S026_tblastn.xml"
result_handle = open(blastXMLfile, "r")
outAlnFile = "S026_homo.fas"
out_handle = open(outAlnFile, "w")
blast_records = NCBIXML.parse(result_handle)
compilResults = []
for record in blast_records:
	for alignment in record.alignments:
		loc = "unknown"
		strain = alignment.title.split("|")[-1].lstrip()
		if "plasmid" in strain:
			loc = "plasmid"
		if "island" in strain:
			loc = "GI"
		strain = strainNameHomogene(strain)
		gi = alignment.title.split("|")[1]
		for hsp in alignment.hsps:
			if hsp.query_start < 75 and hsp.query_end > 250:
				#print hsp
				header = str(">" + strain + "|" + gi + "|" + str(hsp.sbjct_start) + "-" + str(hsp.sbjct_end) +"|" + loc )
				toWrite = header + "\n" + hsp.sbjct + "\n"
				out_handle.write(toWrite)
				break
result_handle.close()
out_handle.close()
