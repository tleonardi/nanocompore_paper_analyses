
library("tidyverse")
library("rtracklayer")
library("biomaRt")

ROOTDIR=system2("git", args=c("rev-parse", "--show-toplevel"), stdout=T)
NANOCOMPORE=paste0(ROOTDIR, "/nanocompore_pipelines/METTL3_KD_polyA/results/nanocompore/")
BASEDIR=paste0(ROOTDIR, "/profiles/METTL3_polyA/")
RESULTS=paste0(ROOTDIR, "/profiles/METTL3_polyA/results/")

# Load nanocompore results
nanocompore <- read_tsv(paste0(NANOCOMPORE, "/out_nanocompore_results.tsv"), col_types="icicccddddddddcicdc")


tot_transcripts <- nanocompore %>% pull(ref_id) %>% unique %>% length

n_sig <- data.frame(Thr=seq(0,1,0.01))
for(i in n_sig$Thr){
	n_sig[n_sig$Thr==i, "Sites"] <- filter(nanocompore, GMM_logit_pvalue<i) %>% nrow
	n_sig[n_sig$Thr==i, "Tx"] <- filter(nanocompore, GMM_logit_pvalue<i) %>% pull(ref_id) %>% unique() %>% length
}

pdf(paste0(RESULTS, "/basic_stats.pdf"))
mutate(n_sig, Ave=Sites/Tx) %>% reshape2::melt(id.var="Thr") %>% filter(variable!="Sites") %>% ggplot(aes(x=Thr, y=value, colour=variable)) + geom_line() + theme_bw() + scale_colour_discrete(name="", labels=c("Number of sig Tx","Number of sig\nsites per tx")) + xlab("p-value threshold (GMM method)") + ylab("Number of sites/transcripts")
dev.off()


# Total number of transcripts: 752
# nanocompore %>% pull(ref_id) %>% unique %>% length

# Applying a p-value threshold of 0.01 (GMM+logit method, see Materials and Methods) we found 6,021 significant sites in 437 distinct transcripts, averaging at 13.8 m6A sites per transcript (Figure 4A,B).  
# mutate(n_sig, Ave=Sites/Tx) %>% filter(Thr==0.01)

#filter(nanocompore, ref_id=="ENST00000646664") %>% dplyr::filter(`GMM_logit_pvalue`<0.01) %>% nrow
#filter(nanocompore, ref_id %in% c("ENST00000425660","ENST00000462494","ENST00000646664")) %>% dplyr::filter(`GMM_logit_pvalue`<0.01) %>% arrange(GMM_logit_pvalue) %>% dplyr::select(pos, genomicPos, ref_kmer, ref_id, GMM_logit_pvalue)

