library("tidyverse")
library("Repitools")
library("GenomicRanges")
library("GenomicAlignments")
library("cowplot")
library("RColorBrewer")
library("reshape2")

ROOTDIR=system2("git", args=c("rev-parse", "--show-toplevel"), stdout=T)
BASEDIR=paste0(ROOTDIR, "/METTL3_polyA/")
RESULTS=paste0(ROOTDIR, "/METTL3_polyA/results/")


file_annot <- tribble(~Name, ~File,
		"Input_Ctrl_rep1", "MOLM13_Input_Ctrl_rep1.bam",
		"Input_sh68_rep1", "MOLM13_Input_sh68_rep1.bam",
		"Input_sh70_rep1", "MOLM13_Input_sh70_rep1.bam",
		"IP_Ctrl_rep1", "MOLM13_IP_Ctrl_rep1.bam",
		"IP_sh68_rep1", "MOLM13_IP_sh68_rep1.bam",
		"IP_sh70_rep1", "MOLM13_IP_sh70_rep1.bam",
		"Input_Ctrl_rep2", "MOLM13_Input_Ctrl_rep2.bam",
		"Input_sh68_rep2", "MOLM13_Input_sh68_rep2.bam",
		"Input_sh70_rep2", "MOLM13_Input_sh70_rep2.bam",
		"IP_Ctrl_rep2", "MOLM13_IP_Ctrl_rep2.bam",
		"IP_sh68_rep2", "MOLM13_IP_sh68_rep2.bam",
		"IP_sh70_rep2", "MOLM13_IP_sh70_rep2.bam") %>%
	      mutate(File=paste0(BASEDIR, "/data/MOLM13_m6A_RIP/", File))
bam <- file_annot$File
names(bam) <- file_annot$Name

# Load GMM peaks
peaks <- read_tsv(paste0(BASEDIR, "data/nanocompore/GMM_pvalue_thr1.bed"), col_names=F, skip=1) %>% 
	 rename(X1="Chr", X2="Start", X3="End", X4="Name", X5="Score", X6="Strand") %>%
	 mutate(Name=paste0(Name, "_", Start))

# Extract Significant and not significant peaks
sig_peaks <- filter(peaks, Score>-log10(0.01))
non_sig_peaks <- top_n(peaks, nrow(sig_peaks), rev(Score))
peaks_subset <- rbind(sig_peaks, non_sig_peaks)

# Identify the midpoint
peaks_mid <- GRanges(seqnames=peaks_subset$Chr, ranges=IRanges(start=peaks_subset$Start+2, end=peaks_subset$Start+2), strand=peaks_subset$Strand, name=peaks_subset$Name)

# Calculate coverage of Sites
rangeW=100
all_peak_counts=list()
for (sname in names(bam)) {
    message(paste0("Processing: ", sname))
    peak_count <- featureScores(bam[[sname]], peaks_mid, up = rangeW, down = rangeW, dist = "base", s.width=100, freq = 10, use.strand=T)
    for (i in 1:length(peak_count)) peak_count@scores[[i]] <- 1e6*peak_count@scores[[i]]
    all_peak_counts[[sname]] <- peak_count
    gc()
}

point_tibble<-function(df,name){
        df1 <-  df@scores[[1]]
        df1 <-  as_tibble(df1) %>%
                mutate(Name= map2(rownames(df1), paste0("##", 1:nrow(df1)), paste0) %>% unlist) %>%
                gather(-Name, key="position", value="Signal") %>%
                mutate(position = as.numeric(position))
        df1 <-  mutate(df1, sample = rep(name, nrow(df1)))
        return(df1)
}
# Summarise the coverage
merged_points <- 
	map2(all_peak_counts, as.list(names(all_peak_counts)), point_tibble) %>% 
	bind_rows() %>% 
	mutate(sample=factor(sample, levels=names(bam))) %>%     
	mutate(Condition=gsub("_rep[12]", "", sample)) %>%
	mutate(IP=gsub("(Input|IP).+", "\\1", sample)) %>%
	mutate(Treatment=gsub(".+(Ctrl|sh[0-9]+).+", "\\1", sample)) %>%
	mutate(Rep=gsub(".+(rep[12])", "\\1", sample)) %>%
	dcast(Name+position+Rep+Treatment~IP, value.var="Signal") %>%
	mutate(Enrichment=log2((IP+1)/(Input+1))) %>%
	group_by(Name, position, Treatment)

# Annotate scores
peaks_scores=mutate(peaks, Sig=Score>17) %>% select(Name, Score, Sig) 
pdf(paste0(RESULTS,"/RIP_coverage_of_nanocompore_sites.pdf"), width=17)
pdf("tmp.pdf", width=17)
ungroup(merged_points) %>% 
	filter(IP+Input>0) %>%
	mutate(Name=gsub("##[0-9]+$", "", Name)) %>% 
	left_join(peaks_scores) %>%
	mutate(Sig=case_when(Sig==T~"Nanocompore p-value<0.01", T~"Nanocompore p-value>=0.01")) %>%
	group_by(position, Treatment, Rep, Sig) %>% 
	summarise(meanEnrichment=mean(Enrichment), sd=sd(Enrichment)) %>% 
	ggplot(aes(x=position, y=meanEnrichment, colour=Treatment)) + geom_line() + facet_grid(Rep~Sig) + geom_ribbon(aes(ymin=meanEnrichment-sd, ymax=meanEnrichment+sd, fill=Treatment), alpha=0.1, size=0.1)
dev.off()

pdf(paste0(RESULTS,"/RIP_vs_nanocompore_scores.pdf"), width=17)
ungroup(merged_points) %>% 
	 mutate(Name=gsub("##[0-9]+$", "", Name)) %>% 
	 left_join(peaks_scores) %>%
	 filter(position>-11, position<11) %>%
	 group_by(Name, Rep, Treatment) %>%
	 summarise(Enrichment=mean(Enrichment), Score=unique(Score), Sig=unique(Sig)) %>%
	 ggplot(aes(y=Score, x=Enrichment, colour=Sig)) +geom_point(alpha=0.5) + facet_grid(Rep~Treatment) + xlab("m6A RIP enrichment ratio\nlog2(IP/Input)") + ylab("Nanocompore p-value (-log10)")
dev.off()


pdf(paste0(RESULTS,"/RIP_vs_nanocompore_boxplot.pdf"), width=17)
ungroup(merged_points) %>% 
	mutate(Name=gsub("##[0-9]+$", "", Name)) %>% 
	left_join(peaks_scores) %>%
	mutate(Sig=case_when(Sig==T~"Nanocompore p-value<0.01", T~"Nanocompore p-value>=0.01")) %>%
	filter(position>-11, position<11) %>%
	group_by(Name, Rep, Treatment) %>%
	summarise(Enrichment=mean(Enrichment), Score=unique(Score), Sig=unique(Sig)) %>%
	ggplot(aes(fill=Sig, y=Enrichment, x=Treatment)) +geom_boxplot() + facet_wrap(~Rep) + xlab("m6A RIP enrichment ratio\nlog2(IP/Input)") + ylab("Nanocompore p-value (-log10)")
dev.off()

pdf("tmp.pdf")
rip_de_peaks_ranges <- GRanges(seqnames=rip_de_peaks$chr, ranges=IRanges(start=rip_de_peaks$chromStart, end=rip_de_peaks$chromEnd), strand=rip_de_peaks$strand)
olaps <- findOverlaps(nanocompore_sig_sites_ranges, rip_de_peaks_ranges)
olap_rows <- subjectHits(olaps) %>% unique
rip_de_peaks[olap_rows, "Nanocompore"] <- "Yes"
rip_de_peaks <- mutate(rip_de_peaks, Nanocompore=case_when(is.na(Nanocompore)~"No", T~Nanocompore))
ggplot(rip_de_peaks, aes(x=diff.log2.fc, y=-diff.lg.fdr, colour=Nanocompore)) + geom_point(alpha=0.2) + ylim(0,20)
dev.off()


nanocomp_counts <- read_tsv(paste0(BASEDIR, "analysis/rip_clip_olaps/nanocompore_counts.txt"))
rip_clip_counts <- read_tsv(paste0(BASEDIR, "analysis/rip_clip_olaps/rip_clip_counts.txt"))

pdf(paste0(RESULTS,"/RIP_vs_nanocompore_boxplot.pdf"), width=12)
rip_clip_counts %>% mutate(wo_nanocompore=Total-w_nanocompore) %>% select(-Total) %>% melt %>%
	mutate(Type=fct_recode(Type,  `m6A miCLIP (Vu et al.)`="CLIP", `m6A meRIP (Barbieri et al.)`="RIP", `m6A meRIP METTL3 dependent peaks (Barbieri et al.)`="RIP_M3")) %>%
	mutate(Validation=fct_recode(Validation, 
				     `All meRIP\nsites`="RIP_only", 
				     `meRIP sites\nwith miCLIP`="With_CLIP",
				     `All miCLIP\nsites`="CLIP_only",
				     `miCLIP sites\nwith meRIP`="With_RIP"
				     )
        ) %>%
	mutate(variable=fct_recode(variable,  `p-value>=0.05`="wo_nanocompore", `p-value<0.05`="w_nanocompore")) %>%
	ggplot(aes(x=Validation, y=value, fill=fct_reorder(variable, rev(value)))) + 
		geom_col(colour="black") + 
		facet_wrap(~Type, ncol=1, scales="free") + 
		scale_fill_manual(name="Nanocompore", values=c("grey", "#de2d26")) +
		xlab("") + ylab("Count") +
		theme_bw()
dev.off()


pdf("tmp.pdf")
nanocomp_counts %>% 
	ggplot(aes(x=Type, y=N)) + geom_bar()
