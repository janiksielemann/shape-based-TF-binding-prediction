import pandas as pd
from argparse import ArgumentParser

input_ = ArgumentParser()
input_.add_argument("-f", dest = "fimo", required = True)
input_.add_argument("-p", dest = "peaks", required = True)
args = input_.parse_args()

peak_file = args.peaks
fimo_file = args.fimo

peak_df = pd.read_csv(peak_file, sep="\t", header = None)
fimo_df = pd.read_csv(fimo_file, sep="\t", comment="#")

peak_df_without_dimer = peak_df.copy()
peak_df_with_exactly_one_fimo = peak_df.copy()
peak_df_with_at_least_two_fimo = peak_df.copy()
fimo_without_dimer = fimo_df.copy()
fimo_without_single = fimo_df.copy()

fimo_df["signal_value"] = 0

for index, row in peak_df.iterrows():
  chromosome_ = row[0].lower()
  start_ = row[1]
  stop_ = row[2]
  all_fimos_within_peak = fimo_df.loc[(fimo_df["sequence_name"].str.lower()==chromosome_) & (fimo_df["start"]>=start_) & (fimo_df["stop"]<=stop_), : ]
  fimo_df.loc[(fimo_df["sequence_name"].str.lower()==chromosome_) & (fimo_df["start"]>=start_) & (fimo_df["stop"]<=stop_), "signal_value"] = row[6]
  
  if len(all_fimos_within_peak) > 1:
    peak_df_without_dimer = peak_df_without_dimer.drop(index)
    peak_df_with_exactly_one_fimo = peak_df_with_exactly_one_fimo.drop(index)
    try:
      fimo_without_dimer = fimo_without_dimer.drop(all_fimos_within_peak.index)
    except KeyError:
      pass
  elif len(all_fimos_within_peak) == 0:
    peak_df_with_exactly_one_fimo = peak_df_with_exactly_one_fimo.drop(index)
    peak_df_with_at_least_two_fimo = peak_df_with_at_least_two_fimo.drop(index)
  elif len(all_fimos_within_peak) == 1:
    peak_df_with_at_least_two_fimo = peak_df_with_at_least_two_fimo.drop(index)
    try:
      fimo_without_single = fimo_without_single.drop(all_fimos_within_peak.index)
    except KeyError:
      pass
  if (index+1)%100 == 0:
    print(str(index+1) + "/" + str(len(peak_df)) + " peaks have been checked")


fimo_df.to_csv("fimo_with_signal_values.tsv", sep = "\t", index= False)
fimo_without_dimer.to_csv("fimo_without_dimer_hits.tsv", sep = "\t", index= False)
fimo_without_single.to_csv("fimo_without_single_hits.tsv", sep = "\t", index= False)

peak_df_without_dimer.to_csv("peaks_without_dimer_hits.tsv", sep="\t",  index= False)
peak_df_with_exactly_one_fimo.to_csv("peak_with_exactly_one_fimo.tsv", sep="\t",  index= False)
peak_df_with_at_least_two_fimo.to_csv("peak_with_at_least_two_fimo.tsv", sep="\t",  index= False)
