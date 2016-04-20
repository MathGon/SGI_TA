from Bio import SeqIO
import os
import matplotlib.pyplot as plt
import numpy as np

"""
Trouver les correspondances entre genomes sources des homologues de 
S025 et S026. Ne garde que les couples dont la distance entre les homologues
est inferieure a 1kb
"""

def generePositionCSV(filein, fileout):
	""" Extrait les positions et genome reference des homologues contenus dans le header du fichier fasta """
	fichierOut = open(fileout, "w")
	
	fichierFas = open(filein)
	records = SeqIO.parse(fichierFas, "fasta")
	for record in records:
		if len(record.id.split("|"))> 1:
			gi = record.id.split("|")[1]
			pos = record.id.split("|")[2]
			start, end  = pos.split("-")
			toWrite =  "\t".join([gi, str(start), str(end), "\n"])
			fichierOut.write(toWrite)
		else:
			print record.id
	fichierFas.close()
	fichierOut.close()

def filtreFastaRecords(fileIn, listeGIagarder):
	""" Ecrit un nouveau fichier fasta contenant uniquement ceux dont le gi est present dans une liste """
	outFileName = fileIn.replace(".fas", "_commun.fas")
	input_handle = open(fileIn, "rU")
	output_handle = open(outFileName, "w")
	records = SeqIO.parse(input_handle, "fasta")
	for record in records:
		gi = record.id.split("|")[1]
		if gi in listeGIagarder:
			SeqIO.write(record, output_handle, "fasta")
	input_handle.close()
	output_handle.close()


# Genere un fichier des gi + positions pour chacun des couples d'homologues
generePositionCSV("S025_homo.fas", "S025_homo_position.csv")
generePositionCSV("S026_homo.fas", "S026_homo_position.csv")

# Compare la liste des gi et trouve les communs
os.system("cut S025_homo_position.csv -f1 | sort > 1.txt")
os.system("cut S026_homo_position.csv -f1 | sort > 2.txt")
os.system("comm -12 1.txt 2.txt > commun.txt")
os.system("rm 1.txt 2.txt")

# Ne garde que les homologues communs dans les fichiers fasta
aGarder = []
with open("commun.txt", "r") as f:
	for line in f: 
		aGarder.append(line.rstrip())
	
filtreFastaRecords("S025_homo.fas", aGarder)
filtreFastaRecords("S026_homo.fas", aGarder)

# Repete l'operation extraction gi et positions mais desormais sur les couples
generePositionCSV("S025_homo_commun.fas", "S025_homo_position.csv")
generePositionCSV("S026_homo_commun.fas", "S026_homo_position.csv")
os.system("sort S025_homo_position.csv > S025_homo_position_sorted.csv")
os.system("sort S026_homo_position.csv > S026_homo_position_sorted.csv")

# Calcul des distances
input_handle1 = open("S025_homo_position_sorted.csv", "rU")
input_handle2 = open("S026_homo_position_sorted.csv", "rU")
lines1 = input_handle1.readlines()
lines2 = input_handle2.readlines()
listDist = []
for i in range(len(lines1)-2):
	start1 = int(lines1[i].split("\t")[1])
	start2 = int(lines2[i].split("\t")[1])
	end1 = int(lines1[i].split("\t")[2])
	end2 = int(lines2[i].split("\t")[2])
	# les deux homologues ne sont pas toujours annotes sur le meme brin
	if end1 < end2:
		dist = start2 - end1
	else:
		dist = start1 - end2
	listDist.append(dist)

#plot distances
listDist = sorted(listDist)
plt.hist(listDist, 200)
plt.savefig('distance_between_homologs.png')




