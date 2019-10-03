library("tidyverse")

ROOTDIR=system2("git", args=c("rev-parse", "--show-toplevel"), stdout=T)

BASEDIR=paste0(ROOTDIR, "/profiles/METTL3_polyA/analysis/motifs/")
RESULTS=paste0(ROOTDIR, "/profiles/METTL3_polyA/results/")

syl <- read_tsv(paste0(BASEDIR, "/sylamer_full.out")) %>% reshape2::melt() %>% mutate(variable=as.numeric(as.character(variable)))
scores <- read_delim(paste0(BASEDIR, "/full_sites_ext.fa.ids.txt"), delim="_", col_names=F)
score_thr <- (filter(scores, X3>-log10(0.01)) %>% nrow )

top_motifs <- filter(syl, value>10) %>% pull(upper) %>% unique
top_auc <- group_by(syl, upper) %>% summarise(AUC=sum(value)) %>% top_n(1, AUC) %>% pull(upper)
top_100_auc <- group_by(syl, upper) %>% summarise(AUC=sum(value)) %>% top_n(100, AUC) %>% pull(upper)

syl <- syl %>% mutate(label=case_when(upper%in%top_auc~upper, T~"Other")) %>% mutate(label=gsub("T", "U", label))
syl[syl$label=="Other", "label"] <- NA

#ggplot(syl, aes(x=variable, y=value, group=upper, colour=label)) + geom_line()

pdf(paste0(RESULTS, "/sylamer.pdf"), width=12)
filter(syl, upper%in%top_100_auc) %>% ggplot(aes(x=variable, y=value, group=upper, colour=label)) + geom_line(size=0.3) +geom_vline(xintercept=score_thr, size=0.5, linetype=2) + theme_classic(24) + xlab("Ranked sequences") + ylab("Hypergeometric p-value\n(-log10)") + labs(colour="Motif") + theme(axis.line = element_line(colour = 'black', size = 0.5))
dev.off() 
