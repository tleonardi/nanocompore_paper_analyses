import csv
from collections import *
from itertools import chain
import numpy as np
import sys

dataset = sys.argv[1]
datasets = [dataset]
# Basedire where nanocompore output is found
basedir="/hps/nobackup/enright/tom/in_silico_dataset/analysis/nanocompore/"
# Nucleotide tolerance to consider a peak a positive hit
tol=3

# Modified positions annotation
mod_pos_file = "/hps/nobackup/enright/nanopore/analyses/nanocompore_paper_analyses/in_silico_dataset/data/simulated_datasets2/dataset_0001/reads_1_pos.tsv"
mod_pos=dict()
with open(mod_pos_file) as tsvfile:
    reader = csv.reader(tsvfile, delimiter='\t')
    next(reader)
    for row in reader:
        mod_pos[row[0]] = [ int(i) for i in row[1].split(";")]

# Parse length of references
ref_file="/hps/nobackup/enright/nanopore/analyses/nanocompore_paper_analyses/in_silico_dataset/data/references/random_guided_weight.fa.fai"
ref_lens=dict()
with open(ref_file) as tsvfile:
    reader = csv.reader(tsvfile, delimiter='\t')
    for row in reader:
        ref_lens[row[0]] = int(row[1])


thresholds=[10**-i for i in np.arange(0,10,0.1)] + [0]
for ds in datasets:
    positives_list=defaultdict(lambda: defaultdict(lambda: defaultdict(set)))
    with open(basedir+ds+"/out_nanocompore_results.tsv") as tsvfile:
        reader = csv.reader(tsvfile, delimiter='\t')
        headers = next(reader)
        nc=namedtuple("nc", headers)
        for row in reader:
            row = nc(*row)
            for thr in thresholds:
                    if np.isnan(float(row.GMM_anova_pvalue)): 
                        GMM_anova_pvalue = 1
                    else:
                        GMM_anova_pvalue = float(row.GMM_anova_pvalue)
                    if np.isnan(float(row.KS_dwell_pvalue)): 
                        KS_dwell_pvalue = 1
                    else:
                        KS_dwell_pvalue = float(row.KS_dwell_pvalue)
                    if np.isnan(float(row.KS_intensity_pvalue)):
                        KS_intensity_pvalue = 1
                    else:
                        KS_intensity_pvalue = float(row.KS_intensity_pvalue)
                    if(GMM_anova_pvalue<=thr):
                        positives_list[row.ref_id]["GMM"][thr].add(int(row.pos))
                    if(KS_dwell_pvalue<=thr):
                        positives_list[row.ref_id]["KS_dwell"][thr].add(int(row.pos))
                    if(KS_intensity_pvalue<=thr):
                        positives_list[row.ref_id]["KS_int"][thr].add(int(row.pos))
    
    
    results = defaultdict(lambda: defaultdict(lambda: defaultdict(int)))
    for ref,real_positives in mod_pos.items():
        ref_len = ref_lens[ref]
        real_negatives = [i for i in range(0,ref_len) if i not in real_positives]
        for met in positives_list[ref].keys():
            for thr in thresholds:
                for real_positive in real_positives:
                    if real_positive in positives_list[ref][met][thr]:
                        # Increment TP counter
                        results[met][thr]["TP"]+=1
                        # Remove real_positive and neighbours from positives_list
                    for i in range(real_positive-tol, real_positive+tol+1):
                        try:
                            positives_list[ref][met][thr].remove(i)
                        except KeyError:
                            pass
                for real_negative in real_negatives:
                    if real_negative in positives_list[ref][met][thr]:
                        results[met][thr]["FP"]+=1
        
    total_real_pos = 0
    for l in mod_pos.values():
        total_real_pos += len(l)
    
    total_sites = 0
    for l in ref_lens.values():
        total_sites+= l
    total_real_neg=total_sites-(total_real_pos*(2*tol+1))
    
    for met in results.keys():
        for thr in thresholds:
            if thr in results[met].keys():
                print(f"{ds}\t{met}\t{thr}\t{results[met][thr]['TP']/total_real_pos}\t{results[met][thr]['FP']/total_real_neg}")
            else:
                print(f"{ds}\t{met}\t{thr}\t0\t0")

