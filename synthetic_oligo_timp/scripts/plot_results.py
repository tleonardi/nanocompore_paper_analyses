#!/bin/python
# usage: plot_results.py dbpath
import sys

dbpath = sys.argv[1]

from nanocompore.SampCompDB import SampCompDB

s=SampCompDB(dbpath+"/out_", fasta_fn=dbpath+"/ref.fa")

p1=s.plot_signal(ref_id='FLuc_Control_Plasmid:88-1805', plot_style='seaborn-whitegrid', figsize=[18,10])
p1[0].savefig(dbpath+"/signal_plot.svg")

p2=s.plot_pvalue(ref_id='FLuc_Control_Plasmid:88-1805', plot_style='seaborn-whitegrid', figsize=[12,8])
p2[0].savefig(dbpath+"/pvalue_plot.svg")

p3=s.plot_position(ref_id='FLuc_Control_Plasmid:88-1805', pos=4, plot_style='seaborn-whitegrid', figsize=[14,10])
p3[0].savefig(dbpath+"/pos4_plot.svg")


