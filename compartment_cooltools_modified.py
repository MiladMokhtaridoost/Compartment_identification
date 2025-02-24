###############################################
## cooltools modified my Milad Mokhtaridoost (milad.mokhtaridoost@sickkids.ca)
##############################################

#import numpy as np
#import matplotlib.pyplot as plt
#import pandas as pd
import os, subprocess

import cooler
import cooltools.lib.plotting

from packaging import version
if version.parse(cooltools.__version__) < version.parse('0.5.1'):
    raise AssertionError("tutorials rely on cooltools version 0.5.1 or higher,"+
                         "please check your cooltools version and update to the latest")

import cooltools

import sys
input = sys.argv[1]
prefix = sys.argv[2]
resolution = sys.argv[3]

#input = 1000000
#prefix = "Pancreas_PatientF_Schmitt"
#resolution = "Pancreas_PatientF_Schmitt"

clr = cooler.Cooler(f"/hpf/largeprojects/pmaass/3D-flow/normalized_data_4DNuc_pipeline/human/{input}/{prefix}.pairs.res{resolution}.cool")
#clr = cooler.Cooler(f"Z:/3D-flow/normalized_data_4DNuc_pipeline/human/{input}/{prefix}.pairs.res{resolution}.cool")


## fasta sequence is required for calculating binned profile of GC conent
if not os.path.isfile('./hg38.fa'):
    ## note downloading a ~1Gb file can take a minute
    subprocess.call('wget https://hgdownload.cse.ucsc.edu/goldenpath/hg38/bigZips/hg38.fa.gz', shell=True)
    subprocess.call('gunzip hg38.fa.gz', shell=True)

import bioframe
bins = clr.bins()[:]
hg38_genome = bioframe.load_fasta('./hg38.fa');
## note the next command may require installing pysam
gc_cov = bioframe.frac_gc(bins[['chrom', 'start', 'end']], hg38_genome)
gc_cov.to_csv('hg38_gc_cov_100kb.tsv',index=False,sep='\t')

# obtain first 3 eigenvectors
cis_eigs = cooltools.eigs_cis(
                        clr,
                        gc_cov,
                        n_eigs=3,
                        )

# cis_eigs[0] returns eigenvalues, here we focus on eigenvectors
eigenvector_track = cis_eigs[1][['chrom','start','end','E1']]

#eigenvector_track.to_csv(f"{input}_compartments{resolution}.csv")
eigenvector_track.to_csv(f"/hpf/largeprojects/pmaass/3D-flow/normalized_data_4DNuc_pipeline/human/{input}/{input}_compartments{resolution}.csv", encoding='utf-8', index=False)
#eigenvector_track.to_csv(f"Z:/3D-flow/normalized_data_4DNuc_pipeline/human/{input}/{input}_compartments{resolution}.csv", encoding='utf-8', index=False)
