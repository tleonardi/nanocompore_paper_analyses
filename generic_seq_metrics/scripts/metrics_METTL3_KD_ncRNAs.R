library(tidyverse)
library(UpSetR)

BASEDIR="/home/tleonardi/programming/bioinformatics/nanocompore_paper_analyses"

counts <- read_tsv(paste0(BASEDIR, "/generic_seq_metrics/results/METTL3_KD_ncRNAs/tx_counts.txt"), col_names=c("Sample", "n_reads", "Tx")) %>%
	mutate(Sample=factor(Sample, levels=c("WT_1", "WT_2", "KD_1", "KD_2")))

# Average number of reads
counts %>% group_by(Tx) %>% summarise(mean=mean(n_reads))
