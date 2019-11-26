library("tidyverse")
library("rtracklayer")
library("biomaRt")

ROOTDIR=system2("git", args=c("rev-parse", "--show-toplevel"), stdout=T)
NANOCOMPORE=paste0(ROOTDIR, "/nanocompore_pipelines/METTL3_KD_polyA/results/nanocompore/")
BASEDIR=paste0(ROOTDIR, "/profiles/METTL3_polyA/")
RESULTS=paste0(ROOTDIR, "/profiles/METTL3_polyA/results/")

# Load nanocompore results
nanocompore <- read_tsv(paste0(NANOCOMPORE, "/out_nanocompore_results.tsv"), col_types="icicccddddddddcicdc")

# Annotate gene names
ensembl = useMart("ensembl",dataset="hsapiens_gene_ensembl")
id2name <- getBM(attributes=c('ensembl_transcript_id', 'hgnc_symbol'), 
      	   filters = 'ensembl_transcript_id', 
      	   values = unique(nanocompore$ref_id), 
      	   mart = ensembl)
nanocompore <- left_join(nanocompore, id2name, by=c("ref_id"="ensembl_transcript_id"))

dplyr::select(nanocompore, pos, genomicPos, ref_id, ref_kmer, GMM_logit_pvalue, hgnc_symbol) %>%
	write_tsv(paste0(RESULTS,"annotated_sites.txt"))

volcano <- nanocompore %>% mutate(LOR=case_when(Logit_LOR=="NC"~0, TRUE~as.numeric(Logit_LOR))) %>% { ggplot(., aes(x=LOR, y=-log10(GMM_logit_pvalue))) + 
	geom_point(size=0.9, alpha=0.8) + 
	ggrepel::geom_text_repel(data=top_n(., 20, -log10(GMM_logit_pvalue)), size=4, aes(label=paste0(hgnc_symbol, "\n", ref_kmer, " (", pos, ")"))) +
	xlab("Logistic regression\nodds ratio") +
	ylab("Nanocompore p-value (-log10)") +
	theme_bw(22)
}


volcano_abs <- nanocompore %>% mutate(LOR=case_when(Logit_LOR=="NC"~0, TRUE~as.numeric(Logit_LOR))) %>% filter(!is.na(GMM_logit_pvalue), GMM_logit_pvalue<0.1) %>% { ggplot(., aes(x=abs(LOR), y=-log10(GMM_logit_pvalue))) + 
	geom_point(alpha=0.8) + 
	ggrepel::geom_text_repel(data=top_n(., 15, -log10(GMM_logit_pvalue)), size=6, aes(label=paste0(hgnc_symbol, "\n", ref_kmer, " (", pos, ")"))) +
	xlab("Logistic regression\nodds ratio") +
	ylab("Nanocompore p-value (-log10)") +
	theme_bw(22)
}


pdf(paste0(RESULTS,"/volcano_plot.pdf"), height=10, width=12)
print(volcano)
dev.off()

pdf(paste0(RESULTS,"/volcano_abs_lor_plot.pdf"), height=10, width=12)
print(volcano_abs)
dev.off()
