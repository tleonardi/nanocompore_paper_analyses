```{r}
library("tidyverse")
```

```{r}
file_annot <- tribble(~Name, ~File,
	      "KD68_Rep1", "data/m6A_KD_68_Rep1.bw",
	      "KD68_Rep2", "data/m6A_KD_68_Rep2.bw",
	      "KD70_Rep1", "data/m6A_KD_70_Rep1.bw",
	      "KD70_Rep2", "data/m6A_KD_70_Rep2.bw",
	      "CT_Rep1", "data/m6A_CT_Rep1.bw",
	      "CT_Rep2", "data/m6A_CT_Rep2.bw")

bam <- list()

for( n in file_annot$Name) {
	bam[[n]] <- rtracklayer::import(as.character(file_annot[file_annot$Name==n, "File"]), format="BigWig")
}


peaks <- import("data/nanocompore_results_gmm_ucsc.txt", format="bed")
peaks_df <- data.frame(Name=peaks$name, pvalue=peaks$score)
start(peaks) <- start(peaks)+2
end(peaks) <- end(peaks)-2
all_counts <- list()
for (sname in names(bam)) {
	message(paste0("Processing: ", sname))
	hits <- findOverlaps(peaks, bam[[sname]])
	peaks$covScore[queryHits(hits)] <- bam[[sname]]$score[subjectHits(hits)]
	names <- peaks$name[queryHits(hits)]
	CovScores <- bam[[sname]]$score[subjectHits(hits)]
	#pvalues <- peaks$score[queryHits(hits)]
	all_counts[[sname]] <- data.frame(Name=names, Score=CovScores)
}

all_counts_df <- map_df(all_counts, I, .id="Sample")

mutate(all_counts_df, Sig=Pvalue<0.01) %>%
ggplot(aes(x=as.factor(Sig), colour=Sample, y=Score)) + geom_boxplot()


mutate(all_counts_df, Sig=Pvalue<0.01) %>%
ggplot(aes(colour=as.factor(Sig), x=Score)) + geom_density() + facet_wrap(~Sample)


mutate(all_counts_df, Sig=ntile(Pvalue,5)) 


enrichment_fc <- mutate(all_counts_df, Condition=gsub("_Rep[12]", "", Sample)) %>% 
	group_by(Name, Condition) %>% 
	summarise(Score=mean(Score)) %>% 
	reshape2::dcast(Name~Condition) %>% 
	mutate(KD68fc=KD68/CT, KD70fc=KD70/CT) %>% 
	select(Name, KD_sh68=KD68fc, KD_sh70=KD70fc) %>% 
	reshape2::melt() %>% 
	left_join(., peaks_df)


mutate(enrichment_fc, Sig=case_when(pvalue<0.01~"p-value<0.01", T~"p-value>=0.01")) %>% 
	mutate(variable=fct_recode(variable, "KD\nshRNA 1"="KD_sh68", "KD\nshRNA 2"="KD_sh70")) %>% 
	ggplot(aes(x=variable, fill=Sig, y=log2(value))) + 
		geom_boxplot(outlier.size=0) +
		#geom_jitter(alpha=0.1, size=0.1, position=position_jitterdodge(jitter.width=0.2)) +
		theme_bw() + 
		ylab("m6A RIP KD vs Ctrl\nfold enrichment") + 
		xlab("Sample") + 
		scale_fill_discrete(name="")

mutate(enrichment_fc, ntiles=ntile(pvalue,30)) %>% ggplot(aes(x=variable, fill=factor(ntiles), y=log2(value))) + geom_boxplot()

ggplot(enrichment_fc, aes(y=-log10(pvalue), x=log2(value))) + geom_point()  + facet_wrap(~variable)


rangeW=500
all_peak_counts=list()
peaks <- import("data/bed_files/sig_sites_GMM_pvalue_thr0.05.ucsc.bed", format="bed")
peaks_mid <- GRanges(seqnames=seqnames(peaks), ranges=IRanges(start=start(peaks)+2, end=start(peaks)+2))
for (sname in names(bam)) {
    message(paste0("Processing: ", sname))
    peak_count <- featureScores(bam[[sname]], peaks_mid, up = rangeW, down = rangeW, dist = "base", freq = 2, use.strand=F)
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

merged_points <-  map2(all_peak_counts, as.list(names(all_peak_counts)), point_tibble) %>% 
	bind_rows() %>% 
	mutate(sample=factor(sample, levels=names(bam))) %>%     
	mutate(Condition=gsub("_Rep[12]", "", sample)) %>%
	group_by(Name, position, Condition) %>%
	summarise(Signal=mean(Signal)) %>% 
	reshape2::dcast(Name+position~Condition, value.var="Signal")

merged_points_fc <- mutate(merged_points, KD68fc=log2((KD68+1)/(CT+1)), KD70fc=log2((KD70+1)/(CT+1))) %>%
	select(Name, position, KD_sh68=KD68fc, KD_sh70=KD70fc) %>%
	reshape2::melt(id.vars=c("Name", "position"))




library("Repitools")
library("GenomicRanges")
library("GenomicAlignments")
library("cowplot")
library("RColorBrewer")

ARSs <-import(paste0(BASEDIR,"/../scripts/ARS_Namshik.bed"), format = "BED")
# Remove ARS1216.5 as it's adjacent to a repetitive rDNA
ARSs <- ARSs[mcols(ARSs)$name!="ARS1216.5",]

ARSs_midpoint <- ARSs
start(ARSs_midpoint) <- floor((end(ARSs)-start(ARSs))/2)+start(ARSs)
end(ARSs_midpoint) <- floor((end(ARSs)-start(ARSs))/2)+start(ARSs)


filt_annot <- mutate(annot, Condition=paste(Condition, Time, IP, Rep, sep="_"))

name_list <- as.character(filt_annot$Sample)
bam=list()
for(i in name_list){
        bam[[i]] <- BAM2GRanges(as.character(filt_annot[filt_annot$Sample==i, "bam"]))
}
names(bam) <- name_list

bam<-GRangesList(bam)
bam<-mergeReplicates(bam, filt_annot$Condition)
gc()


```
