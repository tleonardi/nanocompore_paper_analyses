C
library("tidyverse")
library("rtracklayer")
library("biomaRt")

ROOTDIR=system2("git", args=c("rev-parse", "--show-toplevel"), stdout=T)
BASEDIR=paste0(ROOTDIR, "/METTL3_polyA/")
RESULTS=paste0(ROOTDIR, "/METTL3_polyA/results/")

# Load nanocompore results
nanocompore <- read_tsv(paste0(BASEDIR, "data/nanocompore/nanocompore_results.txt"), col_types="cicccdddddddc")

# Annotate gene names
ensembl = useMart("ensembl",dataset="hsapiens_gene_ensembl")
id2name <- getBM(attributes=c('ensembl_transcript_id', 'hgnc_symbol'), 
      	   filters = 'ensembl_transcript_id', 
      	   values = unique(nanocompore$ref_id), 
      	   mart = ensembl)
nanocompore <- left_join(nanocompore, id2name, by=c("ref_id"="ensembl_transcript_id"))


volcano <- nanocompore %>% { ggplot(., aes(x=LOR, y=-log10(GMM_pvalue))) + 
	geom_point(size=0.9, alpha=0.8) + 
	ggrepel::geom_text_repel(data=top_n(., 20, -log10(GMM_pvalue)), size=4, aes(label=paste0(hgnc_symbol, "\n", ref_kmer, " (", pos, ")"))) +
	xlab("Logistic regression\nodds ratio") +
	ylab("Nanocompore p-value (-log10)") +
	theme_bw(22)
}


volcano_abs <- nanocompore %>% { ggplot(., aes(x=abs(LOR), y=-log10(GMM_pvalue))) + 
	geom_point(size=0.9, alpha=0.8) + 
	ggrepel::geom_text_repel(data=top_n(., 20, -log10(GMM_pvalue)), size=4, aes(label=paste0(hgnc_symbol, "\n", ref_kmer, " (", pos, ")"))) +
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
