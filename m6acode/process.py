# BASEDIR=$(git rev-parse --show-toplevel)
# ANALYSIS="${BASEDIR}/profiles/METTL3_polyA/analysis"
# RESULTS="${BASEDIR}/profiles/METTL3_polyA/results"
# BED="${BASEDIR}/nanocompore_pipelines/METTL3_KD_polyA/results/nanocompore/out_sig_sites_GMM_logit_pvalue_thr_0.05.bed"
# IMAGE="${BASEDIR}/nanocompore_pipelines/singularity_images/nanocompore_4443e07.img"

import sys 
basedir = "/hps/nobackup/enright/tom/nanocompore_paper_analyses/"
results = basedir + "/nanocompore_pipelines/METTL3_KD_polyA/results/"
dbpath = basedir + "/nanocompore_pipelines/METTL3_KD_polyA/results/nanocompore/out_SampComp.db"
fasta = basedir + "/nanocompore_pipelines/METTL3_KD_polyA/results/references/reference_transcriptome.fa"
outdir = basedir + "/m6acode"
from nanocompore.SampCompDB import SampCompDB
from nanocompore.SampComp import SampComp
from nanocompore.Whitelist import Whitelist
from sklearn.preprocessing import StandardScaler
import numpy as np
from collections import *
from nanocompore.TxComp import count_reads_in_cluster
from nanocompore.TxComp import normalised_ordered_counter
from nanocompore.common import *

fn_dict = {"WT": {"WT_1":results+"WT_1/nanopolishcomp/out_eventalign_collapse.tsv", "WT_2":results+"WT_2/nanopolishcomp/out_eventalign_collapse.tsv"}, "KD":{"KD_1":results+"KD_1/nanopolishcomp/out_eventalign_collapse.tsv", "KD_2":results+"KD_2/nanopolishcomp/out_eventalign_collapse.tsv"}}

wl = Whitelist(eventalign_fn_dict=fn_dict, fasta_fn=fasta, select_ref_id = ["ENST00000646664"])
s=SampComp(eventalign_fn_dict=fn_dict, fasta_fn=fasta, outpath=outdir+"/out", whitelist=wl, logit=True)

db=s()

db = SampCompDB(outdir+"/out/out_SampComp.db", fasta_fn=fasta)

poi=1424
tx='ENST00000646664'
model=db[tx][poi]['txComp']['GMM_model']['model']

data = db[tx][poi]['data']
condition_labels = tuple(["WT", "KD"])
sample_labels = ['WT_1', 'WT_2', 'KD_1', 'KD_2']

global_intensity = np.concatenate(([v['intensity'] for v in data[condition_labels[0]].values()]+[v['intensity'] for v in data[condition_labels[1]].values()]), axis=None)
global_dwell = np.concatenate(([v['dwell'] for v in data[condition_labels[0]].values()]+[v['dwell'] for v in data[condition_labels[1]].values()]), axis=None)
global_reads = np.concatenate(([v['reads'] for v in data[condition_labels[0]].values()]+[v['reads'] for v in data[condition_labels[1]].values()]), axis=None)
global_dwell = np.log10(global_dwell)
# Scale the intensity and dwell time arrays
X = StandardScaler().fit_transform([(i, d) for i,d in zip(global_intensity, global_dwell)])
Y = [ k for k,v in data[condition_labels[0]].items() for _ in v['intensity'] ] + [ k for k,v in data[condition_labels[1]].items() for _ in v['intensity'] ]

y_pred = model.predict(X)
logr=dict()
labels=[]
counters = dict()
for lab in sample_labels:
    counters[lab] = Counter(y_pred[[i==lab for i in Y]])
    cluster_counts = count_reads_in_cluster(counters)

for sample,counter in counters.items():
       labels.append(sample)
       ordered_counter = [ counter[i] for i in range(2)]
       total = sum(ordered_counter)
       normalised_ordered_counter = [ i/total for i in ordered_counter ]
       # Loop through ordered_counter and divide each value by the first
       logr[sample] = np.log(normalised_ordered_counter[0]/(1-normalised_ordered_counter[0]))

if(np.mean([v for k,v in logr.items() if "KD" in k])<0):
    mod_cluster=0
else:
    mod_cluster=1








