library(Guitar)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)

library(tidyverse)
library(cowplot)
library(lemon)
library(ggrepel)
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene

f="/home/tleonardi/programming/bioinformatics/nanocompore_paper_analyses/nanocompore_pipelines/METTL3_KD_polyA/results/nanocompore/out_sig_sites_GMM_logit_pvalue_thr_0.05.bed"
ncmp_results_sel <- pipe(paste0("tail -n +2 ", f, "| bedparse convertChr --assembly hg38 --target ucsc")) %>% read_tsv(, col_names=c("chr", "start", "end", "name", "score", "strand"))
# Use all differential sites

sites05 <- ncmp_results_sel %>% with(., GRanges(seqnames=chr, ranges=IRanges(start=start, end=end), strand=strand))
sites001 <- filter(ncmp_results_sel, score>2) %>% with(., GRanges(seqnames=chr, ranges=IRanges(start=start, end=end), strand=strand))

diffSites <- list(sites001)

names(diffSites) <- c("p-value<0.01")

# Generate guitar plot
pdf("/home/tleonardi/programming/bioinformatics/nanocompore_paper_analyses/profiles/METTL3_polyA/results/mRNA_profile.pdf")
GuitarPlot(stGRangeLists = diffSites, txTxdb=txdb, pltTxType=c("mrna"), headOrtail=FALSE, enableCI=FALSE)
dev.off()

#library(RNAModR)
#BuildTx("hg38")
#
#sites001 <- filter(ncmp_results_sel, score>2) %>% with(., GRanges(seqnames=chr, ranges=IRanges(start=start+3, end=start+3), strand=strand))
#mcols(sites001)$id <- "site"
#mcols(sites001)$score <- filter(ncmp_results_sel, score>2) %>% pull(score)
#PlotSectionDistribution(FilterTxLoc(SmartMap(sites001, id="nanocompore"),filter = c("5'UTR", "CDS", "3'UTR")))
