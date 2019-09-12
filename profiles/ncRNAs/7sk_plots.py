#!/bin/python
# usage: plot_results.py dbpath
import sys
basedir = "/hps/nobackup/enright/tom/nanocompore_paper_analyses/"
dbpath = basedir + "/nanocompore_pipelines/METTL3_KD_ncRNAs/results/nanocompore/out_SampComp.db"
fasta = basedir + "/nanocompore_pipelines/METTL3_KD_ncRNAs/results/references/reference_transcriptome.fa"
outdir = basedir + "/profiles/ncRNAs"
from nanocompore.SampCompDB import SampCompDB
s=SampCompDB(db_fn=dbpath, fasta_fn=fasta)

id = "ENST00000636484"
offset= 69

p1=s.plot_signal(ref_id=id, plot_style='seaborn-whitegrid', figsize=[20,10], start=35+offset, end=68+offset)
p1[0].savefig(outdir+"/7sk_mettl3_signal_plot.svg")

p2=s.plot_pvalue(ref_id=id, plot_style='seaborn-whitegrid', figsize=[20,8], tests="GMM_logit_pvalue")
p2[0].savefig(outdir+"/7sk_mettl3_pvalue_plot.svg")

p3=s.plot_position(ref_id=id, pos=41+offset, plot_style='seaborn-whitegrid', figsize=[14,10], pointSize=10)
p3[0].savefig(outdir+"/7sk_mettl3_pos41_plot.svg")

p4=s.plot_position(ref_id=id, pos=63+offset, plot_style='seaborn-whitegrid', figsize=[14,10], pointSize=10)
p4[0].savefig(outdir+"/7sk_mettl3_pos63_plot.svg")

dbpath = basedir + "/nanocompore_pipelines/IVT_7SK/results/nanocompore/out_SampComp.db"
fasta = basedir + "/nanocompore_pipelines/IVT_7SK/results/references/reference_transcriptome.fa"
outdir = basedir + "/profiles/ncRNAs"
from nanocompore.SampCompDB import SampCompDB
s=SampCompDB(db_fn=dbpath, fasta_fn=fasta)

id = "ENST00000636484"
offset= 69

p1=s.plot_signal(ref_id=id, plot_style='seaborn-whitegrid', figsize=[20,10], start=35+offset, end=68+offset)
p1[0].savefig(outdir+"/7sk_ivt_signal_plot.svg")

p2=s.plot_pvalue(ref_id=id, plot_style='seaborn-whitegrid', figsize=[20,8], tests="GMM_logit_pvalue")
p2[0].savefig(outdir+"/7sk_ivt_pvalue_plot.svg")

p3=s.plot_position(ref_id=id, pos=41+offset, plot_style='seaborn-whitegrid', figsize=[14,10], pointSize=10)
p3[0].savefig(outdir+"/7sk_ivt_pos41_plot.svg")

p4=s.plot_position(ref_id=id, pos=63+offset, plot_style='seaborn-whitegrid', figsize=[14,10], pointSize=10)
p4[0].savefig(outdir+"/7sk_ivt_pos63_plot.svg")

p4=s.plot_position(ref_id=id, pos=248+offset, plot_style='seaborn-whitegrid', figsize=[14,10], pointSize=10)
p4[0].savefig(outdir+"/7sk_ivt_pos248_plot.svg")

p4=s.plot_position(ref_id=id, pos=247+offset, plot_style='seaborn-whitegrid', figsize=[14,10], pointSize=10)
p4[0].savefig(outdir+"/7sk_ivt_pos247_plot.svg")

p1=s.plot_signal(ref_id=id, plot_style='seaborn-whitegrid', figsize=[20,10], start=228+offset, end=250+offset)
p1[0].savefig(outdir+"/7sk_ivt_signal_plot_228-250.svg")
