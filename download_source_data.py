import subprocess

# download DAP-seq data
subprocess.run(["wget","http://neomorph.salk.edu/dap_web/pages/dap_data_v4/fullset/dap_download_may2016_peaks.zip"])

subprocess.run(["unzip", "dap_download_may2016_peaks.zip"])


# download A. thaliana genome
subprocess.run(["wget","https://www.arabidopsis.org/download_files/Genes/TAIR10_genome_release/TAIR10_chromosome_files/TAIR10_chr_all.fas"])

