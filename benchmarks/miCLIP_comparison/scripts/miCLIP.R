library(tidyverse)
library(csaw)
library(Repitools)

ROOTDIR=system2("git", args=c("rev-parse", "--show-toplevel"), stdout=T)
BASEDIR=paste0(ROOTDIR, "/benchmarks/miCLIP_comparison/analysis/")
RESULTS=paste0(ROOTDIR, "/benchmarks/miCLIP_comparison/results/")
BAMDIR=paste0(ROOTDIR, "/profiles/METTL3_polyA/analysis/miCLIP/transcriptome_bam/")
PEAKS=paste0(ROOTDIR, "/benchmarks/peak_intersections/")

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
tot_depth <- read_tsv(paste0(ROOTDIR, "/profiles/METTL3_polyA/analysis/miCLIP/tot_depth.txt"), col_names=c("File", "TotDepth")) %>% mutate(File=paste0(File, ".bam"))
tx_depth <- read_tsv(paste0(ROOTDIR, "/profiles/METTL3_polyA/analysis/miCLIP/transcriptome_depth.txt"), col_names=c("File", "TxDepth")) %>% mutate(File=paste0(File, ".bam"))


genome_bams <- data.frame(genomebam=list.files(paste0(ROOTDIR, "/profiles/METTL3_polyA/data/miCLIP_BAMs/"), pattern="bam$")) %>%
	       mutate(File=gsub("-jufastqgz.+", "", genomebam), File=paste0(File, ".bam"))

annot <- left_join(annot, genome_bams)
annot <- mutate(annot, Sample=paste(IP, Genotype, Rep, sep="_"), bam=paste0(BAMDIR,File)) %>% left_join(tot_depth) %>% left_join(tx_depth)

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
genome_bams <- paste0(ROOTDIR, "/profiles/METTL3_polyA/data/miCLIP_BAMs/", annot$genomebam)
fw_only <- readParam(forward=TRUE)
binned <- windowCounts(genome_bams, ext=1, bin=TRUE, width=500, param=fw_only)
binned <- normFactors(binned, se.out=binned)
annot <- mutate(annot, scaling_factor=colData(binned)$norm.factors*colData(binned)$totals)
annot <- mutate(annot, tt=colData(binned)$totals)




peak_files <- tibble(bedfile=paste0(PEAKS, list.files(PEAKS, pattern="bed$"))) %>% mutate(name=gsub(".bed", "", basename(bedfile)))


all_peaks <- dplyr::select(peak_files, bedfile, name) %>% split(.$name) %>% 
        map(~read_tsv(.x[[1]], col_names=c("Tx", "Start", "End")) %>% mutate(Name=.x[[2]])) %>%
        purrr::reduce(rbind) %>%
	mutate(uniquePeakID=paste0("Peak",1:n()))


# Granges is 1-based: pos+1 is first nt. pos+1+5 is last nt of kmer
all_peaks_gr <- with(all_peaks, GRanges(Tx, strand="*", IRanges(Start+1, Start+1+5), id=uniquePeakID, score=0))

# Annotate the midpoint of nanocompore sites
all_peaks_gr_mid <- all_peaks_gr
start(all_peaks_gr_mid) <- floor((end(all_peaks_gr)-start(all_peaks_gr))/2)+start(all_peaks_gr)
end(all_peaks_gr_mid) <- start(all_peaks_gr_mid)
names(all_peaks_gr_mid) <- all_peaks_gr$id
strand(all_peaks_gr_mid) <- "+"



# Compute coverage
peaks_coverage <- vector("list", length(bam))
names(peaks_coverage) <- names(bam)
rangeW=1500
for (sname in names(bam)) {
        message(paste0("Processing: ", sname))
        current <- featureScores(bam[[sname]], all_peaks_gr_mid, up = rangeW, down = rangeW, dist = "base", freq = 50)
	tx_depth <- as.numeric(annot[annot$Sample==sname, "TxDepth"])
	scaling_factor <- as.numeric(annot[annot$Sample==sname, "scaling_factor"])
	# We multiply by tx_depth (to undo the normalisation done by featureScores) and divide by the TMM scaling factor
        for(i in 1:length(current)) current@scores[[i]] <- 1e6*current@scores[[i]]*tx_depth/scaling_factor
        peaks_coverage[[sname]] <- current
        gc()
}

merged <- map2(peaks_coverage, as.list(names(peaks_coverage)), point_tibble_peaks) %>% bind_rows()
rm(peaks_coverage)
gc()
merged <- left_join(all_peaks, merged, by=c("uniquePeakID"="Name"))


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
	summarise(sd=sd(Signal)/sqrt(length(Signal)), Signal=mean(Signal)) %>% 
	ggplot(aes(x=position, y=Signal, colour=Genotype)) + geom_line()  + theme_bw(22) + geom_ribbon(aes(fill=Genotype, ymin=Signal-sd, ymax=Signal+sd), colour=NA, alpha=0.1) + facet_grid(~Ab, ncol=1) + ylab("Signal") + xlab("Position relative to Nanocompore peaks")
dev.off()


pdf(paste0(RESULTS, "/miCLIP_profile_input_norm.pdf"), width=10)
merged %>% rename(sample="samp") %>% 
	separate(samp, into=c("Ab", "Genotype", "Rep"), remove=F) %>%
	group_by(position, Ab, Genotype, Name, uniquePeakID) %>% 
	summarise(Signal=mean(Signal)) %>% 
	reshape2::dcast(position+Genotype+Name+uniquePeakID~Ab, value.var="Signal") %>%
	mutate(Signal=(m6A+1)/(Input+1)) %>%
	group_by(position, Genotype, Name) %>% 
	summarise(sd=sd(Signal, na.rm=T)/sqrt(length(Signal)), Signal=mean(Signal, na.rm=T)) %>% 
	ggplot(aes(x=position, y=Signal, colour=Genotype)) + 
		geom_line()  + 
		theme_bw(24) + 
		geom_ribbon(aes(fill=Genotype, ymin=Signal-sd, ymax=Signal+sd), colour=NA, alpha=0.1) + 
		facet_wrap(~Name) +
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


