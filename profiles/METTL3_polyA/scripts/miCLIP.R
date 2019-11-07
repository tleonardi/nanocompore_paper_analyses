library(tidyverse)
library(csaw)
library(Repitools)

ROOTDIR=system2("git", args=c("rev-parse", "--show-toplevel"), stdout=T)
NANOCOMPORE=paste0(ROOTDIR, "/nanocompore_pipelines/METTL3_KD_polyA/results/nanocompore/")
BASEDIR=paste0(ROOTDIR, "/profiles/METTL3_polyA/")
RESULTS=paste0(ROOTDIR, "/profiles/METTL3_polyA/results/")
BAMDIR=paste0(BASEDIR, "/analysis/miCLIP/transcriptome_bam/")
annot <- tribble(~IP, ~Genotype, ~Rep, ~File,
"Input", "KO", "Rep1", "input-molm13-mt3-ko-total-rna-20180921.bam",
"Input", "WT", "Rep1", "input-molm13-wt-total-rna-20180921.bam",
"m6A", "KO", "Rep1", "m6a-molm13-mt3-ko-1-20180522.bam",
"m6A", "KO", "Rep2", "m6a-molm13-mt3-ko-2-20180522.bam",
"m6A", "WT", "Rep1", "m6a-molm13-wt1-20180515.bam",
"m6A", "WT", "Rep11", "m6a-molm13-wt1-20180522.bam",
"m6A", "WT", "Rep2", "m6a-molm13-wt2-20180515.bam",
"m6A", "WT", "Rep22", "m6a-molm13-wt2-20180522.bam")


point_tibble_peaks <-function(df,name){
        df1 <-  as_tibble(df@scores[[1]])
        df1$Name <- names(df@anno)
        df1 <- df1 %>% gather(-Name, key="position", value="Signal") %>%
              mutate(position = as.numeric(position))
        df1 <-  mutate(df1, sample = name)
        return(df1)
}

# Load total depth and nanocompore transcriptome depth
tot_depth <- read_tsv(paste0(BASEDIR, "/analysis/miCLIP/tot_depth.txt"), col_names=c("File", "TotDepth")) %>% mutate(File=paste0(File, ".bam"))
tx_depth <- read_tsv(paste0(BASEDIR, "/analysis/miCLIP/transcriptome_depth.txt"), col_names=c("File", "TxDepth")) %>% mutate(File=paste0(File, ".bam"))


genome_bams <- data.frame(genomebam=list.files(paste0(BASEDIR, "/data/miCLIP_BAMs/"), pattern="bam$")) %>%
	       mutate(File=gsub("-jufastqgz.+", "", genomebam), File=paste0(File, ".bam"))

annot <- left_join(annot, genome_bams)
annot <- mutate(annot, Sample=paste(IP, Genotype, Rep, sep="_"), bam=paste0(BAMDIR,File)) %>% left_join(tot_depth) %>% left_join(tx_depth)

nanocompore <- read_tsv(paste0(NANOCOMPORE, "/out_nanocompore_results.tsv"), col_types="icicccddddddddcicdc")
nanocompore_ss <- filter(nanocompore, GMM_logit_pvalue<0.01)
# Granges is 1-based: pos+1 is first nt. pos+1+5 is last nt of kmer
nanocompore_ss_gr <- reduce(with(nanocompore_ss, GRanges(ref_id, strand="*", IRanges(pos+1, pos+1+5), id=ref_kmer, score=GMM_logit_pvalue)), min.gapwidth=10)

# Annotate the midpoint of nanocompore sites
nanocompore_ss_gr_mid <- nanocompore_ss_gr
start(nanocompore_ss_gr_mid) <- floor((end(nanocompore_ss_gr_mid)-start(nanocompore_ss_gr_mid))/2)+start(nanocompore_ss_gr_mid)
end(nanocompore_ss_gr_mid) <- start(nanocompore_ss_gr_mid)
names(nanocompore_ss_gr_mid) <- paste(seqnames(nanocompore_ss_gr_mid), start(nanocompore_ss_gr_mid), sep="_")
strand(nanocompore_ss_gr_mid) <- "+"

# Load BAMs
name_list <- as.character(annot$Sample)
bam=list()
for(i in name_list){
        bam[[i]] <- BAM2GRanges(as.character(annot[annot$Sample==i, "bam"]))
	start_sites <- start(bam[[i]])-10
	start_sites[start_sites<1] <- 1
	start(bam[[i]]) <- start_sites
	end(bam[[i]]) <-  start_sites+10
	levs <- seqlevels(bam[[i]])
	levs <- gsub('\\(.+', "", levs)
	bam[[i]] <- renameSeqlevels(bam[[i]], levs)
	#end(bam[[i]]) <- start(bam[[i]])+1
}
names(bam) <- name_list
bam<-GRangesList(bam)
gc()



### ESTIMATE NORMALISATION FACTORS
genome_bams <- paste0(BASEDIR, "/data/miCLIP_BAMs/", annot$genomebam)
fw_only <- readParam(forward=TRUE)
binned <- windowCounts(genome_bams, ext=1, bin=TRUE, width=500, param=fw_only)
binned <- normFactors(binned, se.out=binned)
annot <- mutate(annot, scaling_factor=colData(binned)$norm.factors*colData(binned)$totals)
annot <- mutate(annot, tt=colData(binned)$totals)


# Compute coverage
peaks_coverage <- vector("list", length(bam))
names(peaks_coverage) <- names(bam)
rangeW=1500
for (sname in names(bam)) {
        message(paste0("Processing: ", sname))
        current <- featureScores(bam[[sname]], nanocompore_ss_gr_mid, up = rangeW, down = rangeW, dist = "base", freq = 50)
	tx_depth <- as.numeric(annot[annot$Sample==sname, "TxDepth"])
	scaling_factor <- as.numeric(annot[annot$Sample==sname, "scaling_factor"])
	# We multiply by tx_depth (to undo the normalisation done by featureScores) and divide by the TMM scaling factor
        for(i in 1:length(current)) current@scores[[i]] <- 1e6*current@scores[[i]]*tx_depth/scaling_factor
        peaks_coverage[[sname]] <- current
        gc()
}

merged <- map2(peaks_coverage, as.list(names(peaks_coverage)), point_tibble_peaks) %>% bind_rows()

pdf(paste0(RESULTS, "/miCLIP_profile_splitReps.pdf"))
merged %>% rename(sample="samp") %>% 
	separate(samp, into=c("Ab", "Genotype", "Rep"), remove=F) %>%
	group_by(position, samp, Ab, Genotype, Rep) %>% 
	summarise(Signal=mean(Signal)) %>% 
	ggplot(aes(x=position, y=Signal, colour=Genotype, group=samp)) + geom_line() + facet_wrap(~Ab, ncol=1) + theme_bw(22)
dev.off()

pdf(paste0(RESULTS, "/miCLIP_profile.pdf"))
merged %>% rename(sample="samp") %>% 
	separate(samp, into=c("Ab", "Genotype", "Rep"), remove=F) %>%
	group_by(position, Ab, Genotype, Name) %>% 
	summarise(Signal=mean(Signal)) %>% 
	summarise(sd=sd(Signal)/sqrt(length(Signal)), Signal=mean(Signal)) %>% 
	ggplot(aes(x=position, y=Signal, colour=Genotype)) + geom_line()  + theme_bw(22) + geom_ribbon(aes(fill=Genotype, ymin=Signal-sd, ymax=Signal+sd), colour=NA, alpha=0.1) + facet_wrap(~Ab, ncol=1) + ylab("Signal") + xlab("Position relative to Nanocompore peaks")
dev.off()


pdf(paste0(RESULTS, "/miCLIP_profile_input_norm.pdf"), width=10)
merged %>% rename(sample="samp") %>% 
	separate(samp, into=c("Ab", "Genotype", "Rep"), remove=F) %>%
	group_by(position, Ab, Genotype, Name) %>% 
	summarise(Signal=mean(Signal)) %>% 
	reshape2::dcast(position+Genotype+Name~Ab, value.var="Signal") %>%
	mutate(Signal=(m6A+1)/(Input+1)) %>%
	group_by(position, Genotype) %>% 
	summarise(sd=sd(Signal, na.rm=T)/sqrt(length(Signal)), Signal=mean(Signal, na.rm=T)) %>% 
	ggplot(aes(x=position, y=Signal, colour=Genotype)) + 
		geom_line()  + 
		theme_bw(24) + 
		geom_ribbon(aes(fill=Genotype, ymin=Signal-sd, ymax=Signal+sd), colour=NA, alpha=0.1) + 
		ylab("Signal over Input") + 
		xlab("Position relative to Nanocompore peaks")
dev.off()

# P-value (Mann–Whitney U test): 7.733821e-62
pvalue <- merged %>% rename(sample="samp") %>% 
	separate(samp, into=c("Ab", "Genotype", "Rep"), remove=F) %>%
	filter(abs(position)<100) %>%
	filter(Ab=="m6A") %>%
	group_by(Genotype, Name, Rep) %>% 
	summarise(Signal=mean(Signal)) %>% 
	summarise(Signal=mean(Signal)) %>% 
	reshape2::dcast(Name~Genotype) %>%
	ungroup() %>%
	summarise(w=wilcox.test(KO,WT, data=., paired=F)$p.value)
	

myscale <- function(x){
  (x - mean(x, na.rm=TRUE)) / sd(x, na.rm=TRUE)
}

scaled_agg_repl <- 
	merged %>%
	separate(sample, into=c("Ab", "Genotype", "Rep"), remove=T) %>%
	group_by(position, Name, Ab, Genotype) %>% 
	summarise(Signal=mean(Signal)) %>%
	group_by(Name) %>% 
	mutate(Signal=myscale(Signal))

pdf(paste0(RESULTS, "/miCLIP_heatmap.pdf"))
filter(scaled_agg_repl) %>%
	filter(Ab=="m6A") %>%
	ggplot(aes(x=position, y=Name, fill=Signal)) + 
	geom_tile() + 
	facet_wrap(~Genotype) +
	scale_fill_distiller(palette="GnBu", oob=scales::squish, limits=c(min(scaled_agg_repl$Signal), quantile(scaled_agg_repl$Signal, 0.95, na.rm=T)), direction=-1, name="miCLIP signal\n(z-score)") +
	ylab("") +
	xlab("Position from\nNanocompore peak (nt)") +
	theme_bw(16) + 
	theme(axis.line = ggplot2::element_line(size = 0)) +
	scale_y_discrete(breaks=NULL)
dev.off()



tss_gr_mid <- nanocompore_ss_gr
start(tss_gr_mid) <- 1
end(tss_gr_mid) <- 1
strand(tss_gr_mid) <- "+"
names(tss_gr_mid) <- paste(seqnames(tss_gr_mid), start(tss_gr_mid), sep="_")
lengths <- group_by(nanocompore_ss, ref_id) %>% summarise(len=max(pos))
seqlengths(tss_gr_mid)[lengths$ref_id] <- lengths$len

tss_coverage <- vector("list", length(bam))
names(tss_coverage) <- names(bam)
rangeW=500
for (sname in names(bam)) {
        message(paste0("Processing: ", sname))
        current <- featureScores(bam[[sname]], tss_gr_mid, up = 0, down = rangeW, dist = "base", freq = 50)
	tx_depth <- as.numeric(annot[annot$Sample==sname, "TxDepth"])
	scaling_factor <- as.numeric(annot[annot$Sample==sname, "scaling_factor"])
        for(i in 1:length(current)) current@scores[[i]] <- 1e6*current@scores[[i]]*tx_depth/scaling_factor
        tss_coverage[[sname]] <- current
        gc()
}

point_tibble_percent <-function(df,name){
        df1 <-  as_tibble(df@scores[[1]])
        df1$Name <- names(df@anno)
        df1 <- df1 %>% gather(-Name, key="position", value="Signal")
        df1 <-  mutate(df1, sample = name)
        return(df1)
}


merged_tss <- map2(tss_coverage, as.list(names(tss_coverage)), point_tibble_percent) %>% bind_rows()

tss_coverage_percent <- vector("list", length(bam))
names(tss_coverage_percent) <- names(bam)
rangeW=100
for (sname in names(bam)) {
        message(paste0("Processing: ", sname))
        current <- featureScores(bam[[sname]], tss_gr_mid, up = 0, down = rangeW, dist = "percent", freq = 5)
	tx_depth <- as.numeric(annot[annot$Sample==sname, "TxDepth"])
	scaling_factor <- as.numeric(annot[annot$Sample==sname, "scaling_factor"])
        for(i in 1:length(current)) current@scores[[i]] <- 1e6*current@scores[[i]]*tx_depth/scaling_factor
        tss_coverage_percent[[sname]] <- current
        gc()
}

merged_tss_percent <- map2(tss_coverage_percent, as.list(names(tss_coverage)), point_tibble_percent) %>% bind_rows()


pdf(paste0(RESULTS, "/miCLIP_profile_TSS.pdf"))
merged_tss_percent %>% rename(sample="samp") %>% 
	separate(samp, into=c("Ab", "Genotype", "Rep"), remove=F) %>%
	group_by(position, Ab, Genotype, Name) %>% 
	summarise(Signal=mean(Signal)) %>% 
	reshape2::dcast(position+Genotype+Name~Ab, value.var="Signal") %>%
	mutate(Signal=(m6A+1)/(Input+1)) %>%
	mutate(position=as.numeric(gsub(" %", "", position))) %>%
	group_by(position, Genotype) %>% 
	summarise(sd=sd(Signal)/sqrt(length(Signal)), Signal=mean(Signal)) %>% 
	ggplot(aes(x=position, y=Signal, colour=Genotype)) + 
		geom_line() + 
		geom_point() +
		geom_ribbon(aes(fill=Genotype, ymin=Signal-sd, ymax=Signal+sd), colour=NA, alpha=0.1) + 
		theme_bw(22) + 
		ylab("miCLIP signal") + 
		xlab("Position from TSS (%)")
merged_tss %>% rename(sample="samp") %>% 
	separate(samp, into=c("Ab", "Genotype", "Rep"), remove=F) %>%
	group_by(position, Ab, Genotype, Name) %>% 
	summarise(Signal=mean(Signal)) %>% 
	reshape2::dcast(position+Genotype+Name~Ab, value.var="Signal") %>%
	mutate(Signal=(m6A+1)/(Input+1)) %>%
	mutate(position=as.numeric(gsub(" %", "", position))) %>%
	group_by(position, Genotype) %>% 
	summarise(sd=sd(Signal)/sqrt(length(Signal)), Signal=mean(Signal)) %>% 
	ggplot(aes(x=position, y=Signal, colour=Genotype)) + 
		geom_line() + 
		geom_point() +
		geom_ribbon(aes(fill=Genotype, ymin=Signal-sd, ymax=Signal+sd), colour=NA, alpha=0.1) + 
		theme_bw(22) + 
		ylab("miCLIP signal") + 
		xlab("Position from TSS (nt)")
dev.off()
