# shape-based-TF-binding-prediction
A workflow to learn protein-DNA binding specificity from experimental peak data based on DNA shape.  
<br><br>
Required software:  
- MEME-Suite (v5.0 +)  
- Python (v3.7 +)  

Required Python modules:  
- Pandas  
- Numpy  
- Scikit-learn  

Required data:  
- Peak file in ENCODE narrowPeak format  
- Genome of corresponding organism in fasta format  
<br><br>

## Usage

**Generate peak fasta**  
>generate_peak_fasta.py -g TAIR10_genome.fas -p chr1-5_GEM_events.narrowPeak

**Run MEME-chip on peak sequences**  
>meme-chip -oc meme_out -dna peak_sequences.fasta

**Run FIMO on a thaliana genome**  
>fimo --thresh 5e-4 --motif 1 --max-strand --max-stored-scores 1000000 --oc fimo_out meme_out/combined.meme TAIR10_genome.fas

**Dimer analysis**  
>get_dimer_peaks.py -f fimo_out/fimo.tsv -p chr1-5_GEM_events.narrowPeak

**Generate shape fimo**  
>generate_shape_fimo_forward.py

**Train regressor**  
>train_regressor.py  

**Check prediction (example for the transcription factor HY5)**  

>check_predictions.py -s ACGCTGACGTGGCTGCAT  
>predicted affinity: **87.491**  
>average motif affinity in genome: **11.055650525713467**  

>check_predictions.py -s TTTTTTACGTGGTTTTTT  
>predicted affinity: **7.201**  
>average motif affinity in genome: **11.055650525713467**  
