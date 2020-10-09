import pickle
import pandas as pd
import numpy as np
import translate
from Bio.Seq import reverse_complement
from sklearn import preprocessing
from argparse import ArgumentParser

input_ = ArgumentParser()
input_.add_argument("-f", dest = "fimo", required = True)
input_.add_argument("-g", dest = "genome", required = True)
args = input_.parse_args()

#read genome and safe chromosomes in dictionary
datei = open(args.genome , "r")
genome = datei.readlines()
datei.close()

#generate dictionary with all chromosomes and their corresponding sequence
dic_genome = {}
currentSequence = ""

for line in genome:
	if line.startswith(">"):
		if currentSequence: #check if currentsequence is not empty
			dic_genome[currentChromosome] = currentSequence
		currentChromosome = line.split()[0][1:].lower()
		currentSequence = ""
	else:
		currentSequence += line.strip()
dic_genome[currentChromosome] = currentSequence #write the last sequence in dictionary

# generate list of shapes
all_shapes = ['Stagger', 'Rise', 'Opening', 'Buckle', 'MGW', 'Tilt', 'HelT', 'Roll', 'Shear', 'Slide', 'Stretch', 'ProT', 'Shift']

print("-- reading motif occurrences --")

fimo_tsv_file = args.fimo
whole_fimo = pd.read_csv(fimo_tsv_file, sep="\t", comment="#")
whole_fimo["sequence_name"] = whole_fimo["sequence_name"].str.lower()

# remove chloroplast and mitochondra data
whole_fimo = whole_fimo.drop(whole_fimo[(whole_fimo["sequence_name"]=="chloroplast") | (whole_fimo["sequence_name"]=="mitochondria")].index)

# binding column for classifier
whole_fimo["binding"] = 0
whole_fimo.loc[whole_fimo["signal_value"]>0, "binding"] = 1

min_max_scaler = preprocessing.MinMaxScaler(feature_range=(0,1000))
whole_fimo["signal_value"] = min_max_scaler.fit_transform(whole_fimo["signal_value"].to_numpy().reshape(-1,1))

border_length = 30

# get length of motif + borders
motif_len = len(whole_fimo["matched_sequence"].iloc[0]) + (30*2)

# generate list with headers for the shape train dataframe
header = whole_fimo.columns.tolist()
for head in all_shapes:
  for i in range(0,(motif_len)):
    header.append(head + "_" + str(i))

print("---- generate shape dataframe ----")

#---- generate shape dataframe ----
temp_dict = {}
for index, row in whole_fimo.iterrows():
  # get shape values of sequence
  if row["strand"] == "-":
    extracted_sequence = reverse_complement(dic_genome[row["sequence_name"]][((row['start']-1)-(border_length+2)):(row['stop']+(border_length+2))])
  else:
    extracted_sequence = dic_genome[row["sequence_name"]][((row['start']-1)-(border_length+2)):(row['stop']+(border_length+2))]
  if len(extracted_sequence) < 5:
    continue
  sequence_shape_df = translate.seq_to_shape(extracted_sequence)[2:-2].drop(columns = "EP").to_numpy().round(decimals=3)
  values = row.to_list()
  for i in range(13):
    for list_ in sequence_shape_df:
      values.append(list_[i])
  
  # check if lengths of values and df columns is same
  if (len(values) != len(header)):
    continue
  #append the shape of the sequence to the dap dataframe
  temp_dict[index] = values

  if (index+1)%100 == 0:
    print(str(index+1) + "/" + str(len(whole_fimo)) + " shapes have been calculated")

# generate shape df
shapes_df = pd.DataFrame.from_dict(temp_dict, orient = "index")
shapes_df.columns = header

# drop NAs (shapes that could not be calculated)
shapes_df = shapes_df.dropna(axis=0)

shapes_df.to_csv("shape_fimo_forward.csv", index=False)
