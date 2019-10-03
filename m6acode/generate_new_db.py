#!/bin/bash

basedir = "/hps/nobackup/enright/tom/nanocompore_paper_analyses/"
results = basedir + "/nanocompore_pipelines/METTL3_KD_polyA/results/"
dbpath = basedir + "/nanocompore_pipelines/METTL3_KD_polyA/results/nanocompore/out_SampComp.db"
fasta = basedir + "/nanocompore_pipelines/METTL3_KD_polyA/results/references/reference_transcriptome.fa"
outdir = basedir + "/m6acode"
from nanocompore.SampComp import SampComp
from nanocompore.Whitelist import Whitelist

fn_dict = {"WT": {"WT_1":results+"WT_1/nanopolishcomp/out_eventalign_collapse.tsv", "WT_2":results+"WT_2/nanopolishcomp/out_eventalign_collapse.tsv"}, "KD":{"KD_1":results+"KD_1/nanopolishcomp/out_eventalign_collapse.tsv", "KD_2":results+"KD_2/nanopolishcomp/out_eventalign_collapse.tsv"}}

wl = Whitelist(eventalign_fn_dict=fn_dict, fasta_fn=fasta, select_ref_id = ["ENST00000331789"])
s=SampComp(eventalign_fn_dict=fn_dict, fasta_fn=fasta, outpath=outdir+"/out", whitelist=wl, logit=True)
db=s()

