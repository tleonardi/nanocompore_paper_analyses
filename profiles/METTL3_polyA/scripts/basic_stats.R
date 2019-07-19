
library("tidyverse")
library("rtracklayer")
library("biomaRt")

ROOTDIR=system2("git", args=c("rev-parse", "--show-toplevel"), stdout=T)
BASEDIR=paste0(ROOTDIR, "/profiles/METTL3_polyA/")
RESULTS=paste0(ROOTDIR, "/profiles/METTL3_polyA/results/")

# Load nanocompore results
nanocompore <- read_tsv(paste0(BASEDIR, "data/nanocompore/nanocompore_results.txt"), col_types="cicccdddddddc")


tot_transcripts <- nanocompore %>% pull(ref_id) %>% unique %>% length

n_sig <- data.frame(Thr=seq(0,1,0.01))
for(i in n_sig$Thr){
	n_sig[n_sig$Thr==i, "Sites"] <- filter(nanocompore, GMM_pvalue<i) %>% nrow
	n_sig[n_sig$Thr==i, "Tx"] <- filter(nanocompore, GMM_pvalue<i) %>% pull(ref_id) %>% unique() %>% length
}

pdf(paste0(RESULTS, "/basic_stats.pdf"))
mutate(n_sig, Ave=Sites/Tx) %>% reshape2::melt(id.var="Thr") %>% filter(variable!="Sites") %>% ggplot(aes(x=Thr, y=value, colour=variable)) + geom_line() + theme_bw() + scale_colour_discrete(name="", labels=c("Number of sig Tx","Number of sig\nsites per tx")) + xlab("p-value threshold (GMM method)") + ylab("Number of sites/transcripts")
dev.off()
