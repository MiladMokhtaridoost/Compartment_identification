###############################################
## cooltools modified my Milad Mokhtaridoost 
##############################################

import os, subprocess
import cooler
import cooltools.lib.plotting
import cooltools
import sys
import bioframe

from packaging import version
if version.parse(cooltools.__version__) < version.parse('0.5.1'):
    raise AssertionError("tutorials rely on cooltools version 0.5.1 or higher,"+
                         "please check your cooltools version and update to the latest")

input = sys.argv[1] ##### pathway to .cool data
prefix = sys.argv[2] ##### dataset name 
resolution = sys.argv[3] ##### resolution of HiC data

clr = cooler.Cooler(f"{input}/{prefix}.pairs.res{resolution}.cool")

## fasta sequence is required for calculating binned profile of GC conent
if not os.path.isfile('./hg38.fa'):
    ## note downloading a ~1Gb file can take a minute
    subprocess.call('wget https://hgdownload.cse.ucsc.edu/goldenpath/hg38/bigZips/hg38.fa.gz', shell=True)
    subprocess.call('gunzip hg38.fa.gz', shell=True)

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

eigenvector_track = cis_eigs[1][['chrom','start','end','E1']]

### save as csv file
eigenvector_track.to_csv(f"{input}/{prefix}_compartments{resolution}.csv", encoding='utf-8', index=False)
