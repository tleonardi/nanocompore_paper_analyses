library(tidyverse)
library(UpSetR)

BASEDIR="/home/tleonardi/programming/bioinformatics/nanocompore_paper_analyses"

counts <- read_tsv(paste0(BASEDIR, "/generic_seq_metrics/results/TRM5_KO_polyA/tx_counts.txt"), col_names=c("Sample", "n_reads", "Tx")) %>%
	mutate(Sample=factor(Sample, levels=c("WT_1", "WT_2", "KO_1", "KO_2")))

pdf(paste0(BASEDIR, "/generic_seq_metrics/results/TRM5_KO_polyA/cum_frac.pdf"))
counts %>% group_by(Sample, n_reads) %>% summarise(n_molecules=n()) %>% arrange(Sample, n_reads) %>% mutate(cumsum=cumsum(n_molecules)/sum(n_molecules)) %>%
	ggplot(aes(x=n_reads, y=cumsum, colour=Sample)) + geom_line() +  scale_x_log10() + xlab("Coverage per molecule\n(number of reads)") + ylab("Cumulative fraction of reads") + theme_bw(18)
dev.off()

# approximately 50% of the total number of reads mapping to transcripts with coverage less than or equal to 3x.
#counts %>% group_by(Sample, n_reads) %>% summarise(n_molecules=n()) %>% arrange(Sample, n_reads) %>% mutate(cumsum=cumsum(n_molecules)/sum(n_molecules)) %>% filter(n_reads<11) %>% group_by(n_reads) %>% summarise(cumcum=mean(cumsum)) 

pdf(paste0(BASEDIR, "/generic_seq_metrics/results/TRM5_KO_polyA/tot_reads.pdf"))
counts %>% group_by(Sample) %>% summarise(n_reads=sum(n_reads)) %>% ggplot(aes(x=Sample, y=n_reads, label=n_reads)) + geom_col() + geom_text(aes(y=n_reads+50000)) + ylab("Total number of mapped reads") + theme_bw(18)
dev.off()

pdf(paste0(BASEDIR, "/generic_seq_metrics/results/TRM5_KO_polyA/tot_tx.pdf"))
counts %>% group_by(Sample) %>% summarise(n_tx=n()) %>% ggplot(aes(x=Sample, y=n_tx, label=n_tx)) + geom_col() + geom_text(aes(y=n_tx+100)) + ylab("Total number of transcripts sequences") + theme_bw(18)
dev.off()


intersections <- list(WT1=filter(counts, Sample=="WT_1") %>% pull(Tx), WT2=filter(counts, Sample=="WT_2") %>% pull(Tx), KO_1=filter(counts, Sample=="KO_1") %>% pull(Tx), KO_2=filter(counts, Sample=="KO_2") %>% pull(Tx))

pdf(paste0(BASEDIR, "/generic_seq_metrics/results/TRM5_KO_polyA/intersection_sets.pdf"), width=8)
upset(fromList(intersections), text.scale=1.3)
dev.off()

pdf(paste0(BASEDIR, "/generic_seq_metrics/results/TRM5_KO_polyA/between_reps_scatterplot.pdf"), width=9)
counts %>% group_by(Sample) %>% mutate(n_reads=1e6*n_reads/sum(n_reads)) %>% separate(Sample, c("Condition", "Rep")) %>% reshape2::dcast(Tx+Condition~Rep, value.var="n_reads", fill=0) %>% ggplot(aes(x=`1`, y=`2`, colour=Condition)) + geom_point(size=0.5, alpha=0.3) + geom_abline(slope=1, intercept=0) + theme_bw(20) + scale_x_log10() + scale_y_log10() + xlab("Replicate 1") + ylab("Replicate 2") + geom_smooth(method=lm)
dev.off()


# good correlation between replicates (R2  of 0.757 and  0.937 for WT and KD respectively)
# counts %>% group_by(Sample) %>% mutate(n_reads=1e6*n_reads/sum(n_reads)) %>% separate(Sample, c("Condition", "Rep")) %>% reshape2::dcast(Tx+Condition~Rep, value.var="n_reads", fill=0) %>% group_by(Condition) %>% do(a = lm(`1`~`2`,data=.)) %>% glance(., a)


pdf(paste0(BASEDIR, "/generic_seq_metrics/results/TRM5_KO_polyA/aggregated_scatterplot.pdf"))
counts %>% group_by(Sample) %>% mutate(n_reads=1e6*n_reads/sum(n_reads)) %>% separate(Sample, c("Condition", "Rep")) %>% reshape2::dcast(Tx+Condition~Rep, value.var="n_reads", fill=0) %>% mutate(mean_counts=(`1`+`2`)/2) %>% reshape2::dcast(Tx~Condition, value.var="mean_counts") %>% ggplot(aes(x=WT, y=KO)) + geom_point(size=0.5, alpha=0.2) + geom_abline(slope=1, intercept=0) + theme_bw(18) + scale_x_log10() + scale_y_log10() + geom_smooth(method=lm)
dev.off()

#  we observed very good correlation between WT and KD after averaging the replicates (R2  of 0.969)
# counts %>% group_by(Sample) %>% mutate(n_reads=1e6*n_reads/sum(n_reads)) %>% separate(Sample, c("Condition", "Rep")) %>% reshape2::dcast(Tx+Condition~Rep, value.var="n_reads", fill=0) %>% mutate(mean_counts=(`1`+`2`)/2) %>% reshape2::dcast(Tx~Condition, value.var="mean_counts") %>% do(a = lm(WT~KO,data=.)) %>% glance(., a)

#counts %>% group_by(Sample) %>%  reshape2::dcast(Tx~Sample, value.var="n_reads", fill=0) %>% filter(WT_1>30, WT_2>30, KD_1>30, KD_2>30) %>% nrow
