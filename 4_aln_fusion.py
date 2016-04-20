from Bio import SeqIO

fastaFileS025 = open("S025_homo_commun_annotated_filtred.aln", "r") 
recordsS025 = SeqIO.parse(fastaFileS025, "fasta")

fichierFusion = open("S025S026_fusion.aln", "w")

for recordS025 in recordsS025:
	S025id = recordS025.id.split("|")[1]
	fastaFileS026 = open("S026_homo_commun_annotated_filtred.aln", "r") 
	recordsS026 = SeqIO.parse(fastaFileS026, "fasta")
	for recordS026 in recordsS026:
		S026id = recordS026.id.split("|")[1]
		if S025id == S026id:
			startS025, endS025 = recordS025.id.split("|")[2].split("-")
			startS026, endS026 = recordS026.id.split("|")[2].split("-")
			locs = [startS025, endS025, startS026, endS026]
			header = ">" + recordS026.id.split("|")[0] + "|" + S026id + "|" + min(locs) + "-" + max(locs) + "|" + recordS026.id.split("|")[-1]
			toWrite = str(header + "\n" + recordS025.seq + recordS026.seq  + "\n")
			fichierFusion.write(toWrite)
fastaFileS025.close()
fastaFileS026.close()
fichierFusion.close()


			
		
