#!/bin/python
import numpy as np
from collections import *
from matplotlib import cm                                                                                                           
from matplotlib import colors
from string import ascii_letters

nanocomp_res="/home/tleonardi/programming/bioinformatics/nanocompore_paper_analyses/7sk/data/KD/nanocompore/out_nanocompore_results.tsv"
nanocomp_bed="/home/tleonardi/programming/bioinformatics/nanocompore_paper_analyses/7sk/data/KD/references/reference_transcriptome.bed"
sto="/home/tleonardi/programming/bioinformatics/nanocompore_paper_analyses/ncRNAs_structures/7SK/out.sto"
rfam_id="X04236.1/1-328"

nanocomp_start=0
real_start=52995620
symbol="M"
pval_col=8
thr=0.0001

with open(nanocomp_bed) as f:
    for line in f:
        line=line.split('\t')
        if(line[3] == "ENST00000636484"): nanocomp_start=int(line[1])

offset=nanocomp_start-real_start+2
print(f"offset: {offset}")

pos=dict()
with open(nanocomp_res) as f:
    f.readline()
    for line in f:
        line=line.split('\t')
        if(float(line[pval_col])<thr): 
            pos[int(line[0])+offset] = float(line[pval_col])

print(f"positions to annotate: {pos}")


all=defaultdict(list)
for kmer in pos.keys():
    for p in range(kmer,kmer+5):
        all[p].append(pos[kmer])

min_p={k:-np.log10(min(v)) for k,v in all.items()}

# Make color map
cmap = cm.get_cmap('YlOrRd', len(ascii_letters)) 
norm=colors.Normalize(vmin=min(min_p.values()), vmax=max(min_p.values())+3, clip=True) 
cols=dict()
for i in range(len(ascii_letters)):
    cols[ascii_letters[i]] = cmap(i)

print(cols)
# Center kmers
pos = [k for i in pos.keys() for k in range(i,i+5)]

used_key=list()
with open(sto) as f:
    for line in f:
        line=line.split()
        if(len(line)>1):
            if(line[0]==rfam_id):
                seq=line[1]
                counter=1
                annot=""
                annot_p=""
                for i in seq:
                    if(i!="-"):
                        if counter in pos:
                            annot+=symbol
                            c = cmap(norm(min_p[counter]))
                            for k,v in cols.items():
                                if v==c:
                                    annot_p+=k
                                    used_key.append(k)
                        else:
                            annot+="."
                            annot_p+="."
                        counter+=1
                    else:
                        annot+="."
                        annot_p+="."
print(seq)
print(annot)
print(annot_p)

for k,v in cols.items():
    if k in used_key:
        r=int(v[0]*255)
        g=int(v[1]*255)
        b=int(v[2]*255)
        print(f"{k} rgb:{r},{g},{b}")
#REAL_START="52995620"
#NANOCOMP_start=$(awk '$4=="ENST00000636484"{print $2}' $NANOCOMP_BED)
#
#offset=$(echo $NANOCOMP_start-$REAL_START | bc)

