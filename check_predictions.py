import pandas as pd
import numpy as np
import translate
import sys
from joblib import load
from argparse import ArgumentParser
from sklearn import preprocessing

input_ = ArgumentParser()
input_.add_argument("-b", dest = "border", type = int, default = 4)
input_.add_argument("-s", dest = "sequence", required = True)
args = input_.parse_args()

def get_values_of_shape(sequence_):
  sequence_shape_df = translate.seq_to_shape(sequence_)[2:-2].drop(columns = "EP").to_numpy().round(decimals=3)
  values = []
  for i in range(13):
    for list_ in sequence_shape_df:
      values.append(list_[i])
  return values

sequence_border = args.border
all_shapes = ['Stagger', 'Rise', 'Opening', 'Buckle', 'MGW', 'Tilt', 'HelT', 'Roll', 'Shear', 'Slide', 'Stretch', 'ProT', 'Shift']

fimo_df = pd.read_csv("fimo_with_signal_values.tsv", sep="\t", comment="#")
all_core_motifs = list(set(fimo_df["matched_sequence"].to_list()))

core_motif_len = len(all_core_motifs[0])

#feature sets
features = []
for shape in all_shapes:
  for i in range((30 - sequence_border),(30 + core_motif_len + sequence_border)):
    features.append(shape + "_" + str(i))

#load model
regressor = load("rf_regressor.joblib")

sequence_ = args.sequence
sequence_ = sequence_.upper()

#check if sequence has core motif
if sequence_[(sequence_border+2):-(sequence_border+2)] not in all_core_motifs:
  sys.exit("sequence does not contain core motif")

shape_of_sequence = get_values_of_shape(sequence_)

predicted_affinity = regressor.predict([shape_of_sequence])[0].round(decimals=3)

#normalize signal values the same way they were used for learning and calculate mean
min_max_scaler = preprocessing.MinMaxScaler(feature_range=(0,1000))
fimo_df["signal_value"] = min_max_scaler.fit_transform(fimo_df["signal_value"].to_numpy().reshape(-1,1))
mean_singal_value = np.mean(fimo_df["signal_value"].to_numpy())

print("predicted affinity: " + str(predicted_affinity) + "\naverage motif affinity in genome: " + str(mean_singal_value))
