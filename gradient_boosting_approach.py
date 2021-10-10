import matplotlib
# set non interactive backend so that qsub works
matplotlib.use('Agg')

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

from sklearn.ensemble import GradientBoostingRegressor
from sklearn.model_selection import train_test_split
from sklearn.metrics import average_precision_score
from sklearn.metrics import precision_recall_curve
from sklearn.metrics import plot_precision_recall_curve
from sklearn.metrics import mean_squared_error
from joblib import dump
from argparse import ArgumentParser


input_ = ArgumentParser()
input_.add_argument("-b", dest = "border", type = int, default = 4)
args = input_.parse_args()

# generate list of shapes
all_shapes = ['Stagger', 'Rise', 'Opening', 'Buckle', 'MGW', 'Tilt', 'HelT', 'Roll', 'Shear', 'Slide', 'Stretch', 'ProT', 'Shift']

border_length = args.border

#----- read fimo and peak data, check whether dimer or not, generate learning set ---------
fimo_tsv_file_without_dimer = "fimo_without_dimer_hits.tsv"
fimo_tsv_file_without_single = "fimo_without_single_hits.tsv"

peak_file_single = "peak_with_exactly_one_fimo.tsv"
peak_file_dimer = "peak_with_at_least_two_fimo.tsv"

peaks_single = pd.read_csv(peak_file_single, sep="\t")
peaks_dimer = pd.read_csv(peak_file_dimer, sep="\t")

if len(peaks_single) >= len(peaks_dimer):
  peaks = peaks_single.copy()
  whole_fimo = pd.read_csv(fimo_tsv_file_without_dimer, sep="\t", comment="#")
else:
  peaks = peaks_dimer.copy()
  whole_fimo = pd.read_csv(fimo_tsv_file_without_single, sep="\t", comment="#")

# remove chloroplast and mitochondrai data
whole_fimo = whole_fimo.drop(whole_fimo[(whole_fimo["sequence_name"]=="chloroplast") | (whole_fimo["sequence_name"]=="mitochondria")].index)

# generate whole set
shape_fimo = pd.read_csv("shape_fimo_forward.csv")
whole_shape_fimo = pd.merge(whole_fimo, shape_fimo, on=['motif_id', 'motif_alt_id', 'sequence_name', 'start', 'stop', 'strand','score', 'p-value', 'q-value', 'matched_sequence']).copy()

# generate subset for large fimos
# if dataframe too unbalance, generate new one with ratio 1:5 or 1:3 if df too large
if (len(whole_shape_fimo.query("binding==0"))/len(whole_shape_fimo.query("binding==1")) > 5):
  processed_shape_fimo = whole_shape_fimo.query("binding==1").copy()
  if len(whole_shape_fimo.query("binding==1")) > 20000:
    processed_shape_fimo = processed_shape_fimo.append(whole_shape_fimo.query("binding==0").sample(n = (len(whole_shape_fimo.query("binding==1"))*3),random_state = 123).copy())
  else:
    processed_shape_fimo = processed_shape_fimo.append(whole_shape_fimo.query("binding==0").sample(n = (len(whole_shape_fimo.query("binding==1"))*5),random_state = 123).copy())
  processed_shape_fimo = processed_shape_fimo.sort_index()
else:
  processed_shape_fimo = whole_shape_fimo.copy()


core_motif_len = len(processed_shape_fimo["matched_sequence"].iloc[0])

header = []
for shape in all_shapes:
  for i in range((30 - border_length),(30 + core_motif_len + border_length)):
    header.append(shape + "_" + str(i))

X_shapes = processed_shape_fimo[header]

y_shapes = processed_shape_fimo["signal_value"]
y_shapes_binding = processed_shape_fimo["binding"]

# generate train/test and validation set (80% train and 20% validate) (stratify = y, so that the splitting is precisely balanced) (random_state = seed for reproduction)
X_train, X_test, y_train, y_test = train_test_split(X_shapes, y_shapes, test_size=0.2, stratify = y_shapes_binding, random_state=123)
#X_train.shape, y_train.shape, X_test.shape, y_test.shape

# generate sample weights, so that true positives have the weight of the ratio between true positive and false positive fimos
sample_weight_bound = len(processed_shape_fimo.query("binding == 0")) / len(processed_shape_fimo.query("binding > 0"))
print("the bound sample weight is", sample_weight_bound)
sample_weight = np.array([1 if val == 0 else sample_weight_bound for val in y_train])


boosting_regressor = GradientBoostingRegressor(random_state=42)

boosting_regressor.fit(X_train, y_train, sample_weight = sample_weight)

dump(boosting_regressor, "boosting_regressor.joblib")

# --------- PRC whole subset ------------
y_scores_regressor = boosting_regressor.predict(X_shapes)
y_true_regressor = processed_shape_fimo["binding"]

precision, recall, _ = precision_recall_curve(y_true_regressor, y_scores_regressor)
avg_precision_value_regressor = average_precision_score(y_true_regressor, y_scores_regressor)

plt.step(recall, precision, where='post')
plt.xlabel('Recall')
plt.ylabel('Precision')
plt.ylim([0.0, 1.05])
plt.xlim([0.0, 1.0])
plt.title('Precision-Recall curve\navg precision: {:.4f}'.format(avg_precision_value_regressor))
#plt.show()
plt.savefig("PRC_boosting_regressor_whole_subset_border" + str(border_length) + ".png", dpi=300, format="png")
plt.clf()
plt.close()

# --------- PRC validation set ------------
y_scores_regressor = boosting_regressor.predict(X_test)
y_true_regressor = [1 if val > 0 else 0 for val in y_test]

precision, recall, _ = precision_recall_curve(y_true_regressor, y_scores_regressor)
avg_precision_value_regressor_validation = average_precision_score(y_true_regressor, y_scores_regressor)

plt.step(recall, precision, where='post')
plt.xlabel('Recall')
plt.ylabel('Precision')
plt.ylim([0.0, 1.05])
plt.xlim([0.0, 1.0])
plt.title('Precision-Recall curve\navg precision: {:.4f}'.format(avg_precision_value_regressor_validation))
#plt.show()
plt.savefig("PRC_boosting_regressor_validation_border" + str(border_length) + ".png", dpi=300, format="png")
plt.clf()
plt.close()

final_header = ["AUPRC_boosting_regressor","AUPRC_boosting_validation", "border"]
final_values = [avg_precision_value_regressor, avg_precision_value_regressor_validation, border_length]

final_df = pd.DataFrame([final_values], columns=final_header)
final_df.to_csv("results_boosting_forward.tsv", sep="\t", index=False)
