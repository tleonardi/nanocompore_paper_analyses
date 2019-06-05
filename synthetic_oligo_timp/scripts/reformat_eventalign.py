#!/bin/python
# usage: reformat_eventalign.py eventalign.txt outpath
import sys



#  TAAGGTTAACCTGCAGG GGACTGTAGTCCAACTTGAGGA
#  _________________ 11111--------------222
#  CTGTAAAAAAATGAGGACTGTAGTCCAACTTGAGGACTGT
#  22--       ***33333**************44444**
#
# 1: 1718
# 2: 1737
# 3: 1754
# 4: 1773


unmod=list()
mod=list()
mod_ref_pos=1754
unmod_ref_pos=1737
eventalign_file = sys.argv[1]
outpath = sys.argv[2]

with open(eventalign_file) as infile:
	header = infile.readline()
	unmod.append(header.strip().split('\t'))
	mod.append(header.strip().split('\t'))
	for line in infile:
		line=line.strip().split('\t')
		if(line[0] == 'FLuc_Control_Plasmid:88-1805'):
			if(int(line[1]) in range(unmod_ref_pos-3, unmod_ref_pos+4)):
				line[1] = str(int(line[1]) - unmod_ref_pos+3)
				unmod.append(line)
			elif(int(line[1]) in range(mod_ref_pos-3, mod_ref_pos+4)):
				line[1] = str(int(line[1]) - mod_ref_pos+3)
				mod.append(line)

with open(outpath+'/unmod_eventalign.txt', 'w') as unmodfile:
	for l in unmod:
		unmodfile.write('\t'.join(l)+'\n')

with open(outpath+'/mod_eventalign.txt', 'w') as modfile:
	for l in mod:
		modfile.write('\t'.join(l)+'\n')
