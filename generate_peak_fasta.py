from argparse import ArgumentParser

input_ = ArgumentParser()
input_.add_argument("-g", dest = "genome", required = True)
input_.add_argument("-p", dest = "peaks", required = True)
args = input_.parse_args()

print("--- reading chromosomes ---")

#read genome and safe chromosomes in dictionary
temp_file = open(args.genome , "r")
genome = temp_file.readlines()
temp_file.close()

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

#give correct key names
dic_genome["chr1"] = dic_genome.pop("1")
dic_genome["chr2"] = dic_genome.pop("2")
dic_genome["chr3"] = dic_genome.pop("3")
dic_genome["chr4"] = dic_genome.pop("4")
dic_genome["chr5"] = dic_genome.pop("5")

print("Chromosomes in genome data: " + str(list(dic_genome.keys())))

#read file and safe as list
temp_file = open(args.peaks , "r")
peaks = temp_file.readlines()
temp_file.close()

#generate list for new file
ls_final = []
for line in peaks:
	line = line.split("\t") #line[0] = Chromosome, line[1] = start, line[2] = stop
	ls_final.append(">" + line[0]+":"+line[1]+"-"+line[2])
	chrX = line[0].lower()
	ls_final.append(dic_genome[chrX][(int(line[1])-1):int(line[2])])

#write in new output file (fasta)
temp_file = open("peak_sequences.fasta" ,"w")
temp_file.write("\n".join(ls_final))
temp_file.close()
