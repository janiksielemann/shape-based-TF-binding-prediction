import pandas as pd
from argparse import ArgumentParser

input_ = ArgumentParser()
input_.add_argument("-f", dest = "fimo_file", required = True)
input_.add_argument("-p", dest = "peak_file", required = True)
args = input_.parse_args()

final_df = pd.DataFrame(columns=["total_fimo_hits","total_peaks","peaks_hit","peaks_not_hit"])

chromosomes = ["chr1","chr2","chr3","chr4","chr5"]
total_fimo_hits = 0 
total_peaks = 0
peaks_hit = 0
peaks_not_hit = 0

for chromosome_ in chromosomes:
  
  df_fimo = pd.read_csv(args.fimo_file, sep="\t", comment="#")
  df_fimo["sequence_name"] = df_fimo["sequence_name"].replace({"1":"chr1", "2":"chr2", "3":"chr3", "4":"chr4", "5":"chr5"})
  df_fimo = df_fimo.query("sequence_name == @chromosome_")
  
  df_peaks = pd.read_csv(args.peak_file, sep="\t", header=None)
  df_peaks.columns =["chromosome", "start", "end", "name", "score", "strand", "signalValue", "pValue", "qValue", "peak"]
  df_peaks = df_peaks.query("chromosome == @chromosome_")
  
  half_length_of_motif = int(len(df_fimo.iloc[0]["matched_sequence"])/2)
  all_fimo_indices = df_fimo["start"].values + half_length_of_motif
  
  def check_hit(row, indices = pd.Series(all_fimo_indices)):
    return indices.between(row["start"],row["end"]).any()
  
  df_peaks["hit"] = df_peaks.apply(check_hit, axis = 1)
  
  
  total_fimo_hits += len(df_fimo) 
  total_peaks += len(df_peaks)
  peaks_hit += len(df_peaks.query("hit == True"))
  peaks_not_hit += len(df_peaks.query("hit == False"))
  

final_df = pd.DataFrame([[total_fimo_hits,total_peaks,peaks_hit,peaks_not_hit]], columns=final_df.columns)

final_df.to_csv("results_motif_precision_recall.csv")

