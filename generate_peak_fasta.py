from argparse import ArgumentParser

input_ = ArgumentParser()
input_.add_argument("-g", dest = "genome", required = True)
input_.add_argument("-p", dest = "peaks", required = True)
args = input_.parse_args()

print("--- reading chromosomes ---")

#read genome and safe chromosomes in dictionary
datei = open(args.genome , "r")
genome = datei.readlines()
datei.close()

#generate dictionary with all chromosomes and their corresponding sequence
dic_genome = {}
currentSequence = ""

for line in genome:
	if line.startswith(">"):
		if currentSequence: #check if current sequence is not empty
			dic_genome[currentChromosome] = currentSequence
		currentChromosome = line.split()[0][1:].lower()
		currentSequence = ""
	else:
		currentSequence += line.strip()
dic_genome[currentChromosome] = currentSequence #write the last sequence in dictionary

print("Chromosomes in genome data: " + str(dic_genome.keys()))

#read file and safe as list
datei = open(args.peaks , "r")
lines = datei.readlines()
datei.close()

#generate list for new file
ls_final = []
for line in lines[1:]:
	line = line.split("\t") #line[0] = Chromosome, line[1] = start, line[2] = stop
	ls_final.append(">" + line[0]+":"+line[1]+"-"+line[2])
	chrX = line[0].lower()
	ls_final.append(dic_genome[chrX][(int(line[1])-1):int(line[2])])

#write in new output file (fasta)
datei = open("peak_sequences.fasta" ,"w")
datei.write("\n".join(ls_final))
datei.close()
