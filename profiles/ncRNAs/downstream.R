library(tidyverse)
ROOTDIR=system2("git", args=c("rev-parse", "--show-toplevel"), stdout=T)
BASEDIR=paste0(ROOTDIR, "/profiles/ncRNAs/")
DATA=paste0(ROOTDIR, "nanocompore_pipelines/METTL3_KD_ncRNAs/results/nanocompore")

kd <- read_tsv(paste0(ROOTDIR, "/nanocompore_pipelines/METTL3_KD_ncRNAs/results/nanocompore/out_nanocompore_results.tsv"))
ivt <- read_tsv(paste0(ROOTDIR, "/nanocompore_pipelines/IVT_7SK/results/nanocompore/out_nanocompore_results.tsv"))
dkc1 <- read_tsv(paste0(ROOTDIR, "/nanocompore_pipelines/DKC1_ncRNAs/results/nanocompore/out_nanocompore_results.tsv"))

offsets <- tribble(~ref_id, ~offset, ~name,
		   "ENST00000636484", 69, "7SK",
		   "ENST00000602361", 51, "RMRP",
		   "ENST00000618664", 86, "U2",
	   	   "ENST00000516869", 0, "RPPH1")


kd_logit <- select(kd, pos, ref_id, METTL3_KD=GMM_logit_pvalue) %>% mutate(METTL3_KD=case_when(is.na(METTL3_KD)~1, T~METTL3_KD))
ivt_logit <- select(ivt, pos, ref_id, IVT=GMM_logit_pvalue) %>% mutate(IVT=case_when(is.na(IVT)~1, T~IVT))
dkc1_logit <- select(dkc1, pos, ref_id, DKC1_KD=KS_intensity_pvalue) %>% mutate(DKC1_KD=case_when(is.na(DKC1_KD)~1, T~DKC1_KD))

all_datasets <- full_join(kd_logit, ivt_logit, by=c("ref_id", "pos")) %>% 
	full_join(dkc1_logit, by=c("ref_id", "pos")) %>%
	left_join(offsets) %>% mutate(pos=pos-offset, ref_id=paste0(name, " (", ref_id, ")")) %>% select(-offset, -name)

## Logit method
pdf("combined_samples.pdf",width=14)
all_datasets %>%
	reshape2::melt(id.vars=c("pos", "ref_id")) %>%
	ggplot(aes(x=pos, y=-log10(value), colour=variable)) + 
		geom_line(size=0.4) + 
		facet_wrap(~ref_id, scales="free", ncol=1) +
		geom_hline(yintercept=-log10(0.01), linetype=2) +
		xlab("Position along RNA") +
		ylab("-log10(p-value)") +
		scale_color_discrete(name="") +
		scale_x_continuous(breaks=seq(0,max(kd_logit$pos), 20)) +
		theme_bw()
dev.off()


pdf("combined_samples_capped20.pdf",width=18)
all_datasets %>%
	reshape2::melt(id.vars=c("pos", "ref_id")) %>%
	mutate(value=case_when(value<10e-20~10e-20, TRUE~value)) %>%
	ggplot(aes(x=pos, y=-log10(value), colour=variable)) + 
		geom_line(size=0.4) + 
		facet_wrap(~ref_id, ncol=1, scales="free") +
		geom_hline(yintercept=-log10(0.01), linetype=2) +
		xlab("Position along RNA") +
		ylab("-log10(p-value)") +
		scale_color_discrete(name="") +
		scale_x_continuous(breaks=seq(0,max(kd_logit$pos), 20)) +
		theme_bw()
dev.off()

pdf("split_samples.pdf",width=18)
all_datasets %>%
	reshape2::melt(id.vars=c("pos", "ref_id")) %>%
	ggplot(aes(x=pos, y=-log10(value), colour=variable)) + 
		geom_line(size=0.4) + 
		facet_grid(variable~ref_id, scales="free") +
		geom_hline(yintercept=-log10(0.01), linetype=2) +
		xlab("Position along RNA") +
		ylab("-log10(p-value)") +
		scale_color_discrete(name="") +
		scale_x_continuous(breaks=seq(0,max(kd_logit$pos), 20)) +
		theme_bw()
dev.off()

pdf("dkc1_all_tests.pdf")
select(dkc1, pos, ref_id, GMM_logit_pvalue, KS_dwell_pvalue, KS_intensity_pvalue) %>% 
left_join(offsets) %>% mutate(pos=pos-offset, ref_id=paste0(name, " (", ref_id, ")")) %>% select(-offset, -name) %>%
reshape2::melt(id.vars=c("pos", "ref_id")) %>% mutate(value=case_when(is.na(value)~1, T~value)) %>%
	ggplot(aes(x=pos, y=-log10(value), colour=variable)) + 
		geom_line(size=0.4) + 
		facet_wrap(~ref_id, scales="free", ncol=1) +
		geom_hline(yintercept=-log10(0.01), linetype=2) +
		xlab("Position along RNA") +
		ylab("-log10(p-value)") +
		scale_color_discrete(name="") +
		scale_x_continuous(breaks=seq(0,max(dkc1$pos), 20)) +
		theme_bw()

select(dkc1, pos, ref_id, GMM_logit_pvalue, KS_dwell_pvalue, KS_intensity_pvalue) %>% 
left_join(offsets) %>% mutate(pos=pos-offset, ref_id=paste0(name, " (", ref_id, ")")) %>% select(-offset, -name) %>%
reshape2::melt(id.vars=c("pos", "ref_id")) %>% mutate(value=case_when(is.na(value)~1, T~value)) %>%
	ggplot(aes(x=pos, y=-log10(value), colour=variable)) + 
		geom_line(size=0.4) + 
		facet_wrap(~ref_id+variable, scales="free", ncol=3) +
		geom_hline(yintercept=-log10(0.01), linetype=2) +
		xlab("Position along RNA") +
		ylab("-log10(p-value)") +
		scale_color_discrete(name="") +
		scale_x_continuous(breaks=seq(0,max(dkc1$pos), 20)) +
		theme_bw()

select(dkc1, pos, ref_id, ref_kmer, GMM_logit_pvalue, KS_dwell_pvalue, KS_intensity_pvalue) %>% 
left_join(offsets) %>% mutate(pos=pos-offset, ref_id=paste0(name, " (", ref_id, ")")) %>% select(-offset, -name) %>%
filter(pos %in% seq(240,255), grepl("7SK", ref_id)) %>%
mutate(ref_kmer=gsub("T", "U", ref_kmer), pos=paste0(pos, "\n", ref_kmer)) %>% select(-ref_kmer) %>%
reshape2::melt(id.vars=c("pos", "ref_id")) %>% mutate(value=case_when(is.na(value)~1, T~value)) %>%
{
	ggplot(., aes(x=pos, y=-log10(value), colour=variable, group=ref_id)) + 
		geom_line(size=0.4) + 
		facet_wrap(~ref_id+variable, scales="free", ncol=1) +
		geom_hline(yintercept=-log10(0.01), linetype=2) +
		xlab("Position along RNA") +
		ylab("-log10(p-value)") +
		scale_color_discrete(name="") +
		theme_bw()
}

dev.off()


# 24 sig sites at the 0.01
# filter(all_datasets, METTL3_KD<0.01)  %>% data.frame %>% nrow

select(dkc1, pos, ref_id, ref_kmer, cluster_counts, GMM_logit_pvalue, KS_dwell_pvalue, KS_intensity_pvalue) %>%      
left_join(offsets) %>% mutate(pos=pos-offset, ref_id=paste0(name, " (", ref_id, ")")) %>% select(-offset, -name) %>%
filter(pos %in% seq(310,330), grepl("7SK", ref_id)) %>% as.data.frame

select(dkc1, pos, ref_id, ref_kmer, cluster_counts, GMM_logit_pvalue, KS_dwell_pvalue, KS_intensity_pvalue) %>%      
left_join(offsets) %>% mutate(pos=pos-offset, ref_id=paste0(name, " (", ref_id, ")")) %>% select(-offset, -name) %>%
filter(pos %in% seq(240,255), grepl("7SK", ref_id)) %>% as.data.frame

select(dkc1, pos, ref_id, ref_kmer, cluster_counts, GMM_logit_pvalue, KS_dwell_pvalue, KS_intensity_pvalue) %>%      
left_join(offsets) %>% mutate(pos=pos-offset, ref_id=paste0(name, " (", ref_id, ")")) %>% select(-offset, -name) %>%
filter(grepl("RMRP", ref_id), GMM_logit_pvalue<0.01)
