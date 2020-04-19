import os
import numpy as np
import pandas as pd

all_peaks_file = "data/allPeaks_light.hg19.10.bed"

website_TFs_file = "data/website.TFs.srx"
website_TFs = []
with open(website_TFs_file) as f:
    for line in f:
        website_TFs.append(line.strip())

filtered_peaks_dir = "data/TFs_peaks_per_srx"
os.makedirs(filtered_peaks_dir, exist_ok=True)

with open(all_peaks_file, 'r') as fin:
    for line in fin:
        res = line.strip().split('\t')
        if res[3] in website_TFs:
            with open(os.path.join(filtered_peaks_dir, res[3]+'.bed'), 'a') as fout:
                fout.write(line)
