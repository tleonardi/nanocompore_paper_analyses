library(tidyverse)

BASEDIR="/home/tleonardi/programming/bioinformatics/nanocompore_paper_analyses"

counts <- read_tsv(paste0(BASEDIR, "/generic_seq_metrics/results/METTL3_KD_polyA/tx_counts.txt"), col_names=c("Sample", "n_reads", "Tx")) %>%
	mutate(Sample=factor(Sample, levels=c("WT_1", "WT_2", "KD_1", "KD_2")))

pdf(paste0(BASEDIR, "/generic_seq_metrics/results/METTL3_KD_polyA/cum_frac.pdf"))
counts %>% group_by(Sample, n_reads) %>% summarise(n_molecules=n()) %>% arrange(Sample, n_reads) %>% mutate(cumsum=cumsum(n_molecules)/sum(n_molecules)) %>%
	ggplot(aes(x=n_reads, y=cumsum, colour=Sample)) + geom_line() +  scale_x_log10() + xlab("Coverage per molecule\n(number of reads)") + ylab("Cumulative fraction of reads") + theme_bw()
dev.off()

pdf(paste0(BASEDIR, "/generic_seq_metrics/results/METTL3_KD_polyA/tot_reads.pdf"))
counts %>% group_by(Sample) %>% summarise(n_reads=sum(n_reads)) %>% ggplot(aes(x=Sample, y=n_reads, label=n_reads)) + geom_col() + geom_text(aes(y=n_reads+50000)) + ylab("Total number of mapped reads") + theme_bw()
dev.off()

pdf(paste0(BASEDIR, "/generic_seq_metrics/results/METTL3_KD_polyA/tot_tx.pdf"))
counts %>% group_by(Sample) %>% summarise(n_tx=n()) %>% ggplot(aes(x=Sample, y=n_tx, label=n_tx)) + geom_col() + geom_text(aes(y=n_tx+500)) + ylab("Total number of transcripts sequences") + theme_bw()
dev.off()

library(UpSetR)

intersections <- list(WT1=filter(counts, Sample=="WT_1") %>% pull(Tx), WT2=filter(counts, Sample=="WT_2") %>% pull(Tx), KD_1=filter(counts, Sample=="KD_1") %>% pull(Tx), KD_2=filter(counts, Sample=="KD_2") %>% pull(Tx))

pdf(paste0(BASEDIR, "/generic_seq_metrics/results/METTL3_KD_polyA/intersection_sets.pdf"))
upset(fromList(intersections))
dev.off()

pdf(paste0(BASEDIR, "/generic_seq_metrics/results/METTL3_KD_polyA/between_reps_scatterplot.pdf"), width=14)
counts %>% group_by(Sample) %>% mutate(n_reads=1e6*n_reads/sum(n_reads)) %>% separate(Sample, c("Condition", "Rep")) %>% reshape2::dcast(Tx+Condition~Rep, value.var="n_reads", fill=0) %>% ggplot(aes(x=`1`, y=`2`)) + geom_point(size=0.5, alpha=0.2) + facet_wrap(~Condition) + geom_abline(slope=1, intercept=0) + theme_bw() + scale_x_log10() + scale_y_log10() + xlab("Replicate 1") + ylab("Replicate 2")
dev.off()

pdf(paste0(BASEDIR, "/generic_seq_metrics/results/METTL3_KD_polyA/aggregated_scatterplot.pdf"))
counts %>% group_by(Sample) %>% mutate(n_reads=1e6*n_reads/sum(n_reads)) %>% separate(Sample, c("Condition", "Rep")) %>% reshape2::dcast(Tx+Condition~Rep, value.var="n_reads", fill=0) %>% mutate(mean_counts=(`1`+`2`)/2) %>% reshape2::dcast(Tx~Condition, value.var="mean_counts") %>% ggplot(aes(x=WT, y=KD)) + geom_point(size=0.5, alpha=0.2) + geom_abline(slope=1, intercept=0) + theme_bw() + scale_x_log10() + scale_y_log10() 
dev.off()
