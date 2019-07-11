
nanocomp_res_kd="/home/tleonardi/programming/bioinformatics/nanocompore_paper_analyses/7sk/data/KD/nanocompore/out_nanocompore_results.tsv"
nanocomp_res_ivt="/home/tleonardi/programming/bioinformatics/nanocompore_paper_analyses/7sk/data/IVT/nanocompore/out_nanocompore_results.tsv"
nanocomp_res_dkc="/home/tleonardi/programming/bioinformatics/nanocompore_paper_analyses/7sk/data/DKC1/nanocompore/out_nanocompore_results.tsv"
nanocomp_bed="/home/tleonardi/programming/bioinformatics/nanocompore_paper_analyses/7sk/data/IVT/references/reference_transcriptome.bed"                                                                                                                                              
template="7sk.sto" 
tx="ENST00000636484"
thr="0.01"
pval_col=8
real_start=52995620

r2r  --GSC-weighted-consensus RF00100.stockholm.txt $template 3 0.97 0.9 0.75 4 0.97 0.9 0.75 0.5 0.1
cat <(head -n -1 $template) 7sk_custom_annots.txt <(echo //) > 7sk_consensus_structure.sto
cat <(head -n -1 $template) 7sk_custom_annots.txt <(python create_annotations.py $nanocomp_res_kd $nanocomp_bed $template $tx 0.01 $pval_col $real_start) <(echo //) > 7sk_mettl3.sto
cat <(head -n -1 $template) 7sk_custom_annots.txt <(python create_annotations.py $nanocomp_res_ivt $nanocomp_bed $template $tx 0.01 $pval_col $real_start) <(echo //) > 7sk_ivt.sto
cat <(head -n -1 $template) 7sk_custom_annots.txt <(python create_annotations.py $nanocomp_res_dkc $nanocomp_bed $template $tx 0.1 $pval_col $real_start) <(echo //) > 7sk_dkc1.sto

r2r 7sk.meta 7sk.pdf
