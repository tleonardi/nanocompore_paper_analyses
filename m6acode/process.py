import numpy as np
from nanocompore.SampCompDB import SampCompDB
from sklearn.preprocessing import StandardScaler
from collections import *
from nanocompore.common import *

db = SampCompDB("out/out_SampComp.db", fasta_fn="reference_transcriptome.fa")
modified_reads=dict()
poi_list = [1533, 650, 1322]
outfile = open("out.tsv", "w")
for poi in poi_list:
    tx='ENST00000331789'
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
    y_pred_prob = model.predict_proba(X)
    y_pred = model.predict(X)
    logr=dict()
    labels=[]
    counters = dict()
    for lab in sample_labels:
        counters[lab] = Counter(y_pred[[i==lab for i in Y]])
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
    for read, p, lab in zip(global_reads, y_pred_prob, Y):
        outfile.write("\t".join([lab, str(poi), read, str(p[mod_cluster])])+"\n")
outfile.close()
