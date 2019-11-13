library(tidyverse)
library(ggpubr)
library(UpSetR)
library(SuperExactTest)
ROOTDIR=system2("git", args=c("rev-parse", "--show-toplevel"), stdout=T)
BASEDIR=paste0(ROOTDIR, "/m6acode/")

# [1533, 650, 1322] 
p <- read_tsv(paste0(BASEDIR,"out_actin/out.tsv"), col_names=c("lab", "pos", "read", "p")) %>% mutate(lab=gsub("_[12]", "", lab), pos=paste0("A",pos)) 

wt <- filter(p, lab=="WT") %>%
	reshape2::dcast(read~pos) %>% 
	select(-read) %>% 
	filter(!is.na(A650), !is.na(A1322), !is.na(A1533)) %>%
	ggheatmap::ggheatmap(colorPalette="GnBu", orderCol=F, scaleName="p modified") +
	ggtitle("WT")

kd <- 	filter(p, lab=="KD") %>%
	reshape2::dcast(read~pos) %>% 
	select(-read) %>% 
	filter(!is.na(A650), !is.na(A1322), !is.na(A1533)) %>%
	ggheatmap::ggheatmap(colorPalette="GnBu", orderCol=F, scaleName="p modified") +
	ggtitle("KD")

thr=0.75
int_tibble <- filter(p, lab=="WT") %>%
	reshape2::dcast(read~pos) %>% 
	filter(!is.na(A650), !is.na(A1322), !is.na(A1533)) %>% as_tibble()

intersections <- list(Pos650=filter(int_tibble, A650>thr) %>% pull(read), 
		      Pos1322=filter(int_tibble, A1322>thr) %>% pull(read), 
		      Pos1533=filter(int_tibble, A1533>thr) %>% pull(read))
n=filter(p, lab=="WT", p>thr) %>% pull(read) %>% unique %>% length

#mat <- filter(p, lab=="WT") %>%
#	reshape2::dcast(read~pos) %>% 
#	filter(!is.na(A650), !is.na(A1322), !is.na(A1533)) %>%
#	mutate(A650=1*(A650>thr), A1322=1*(A1322>thr), A1533=1*(A1533>thr))



mset2df <- function(mset){
	set_names=mset$set.names
	olap_sizes <- mset$overlap.sizes
	olap_names <- names(olap_sizes) %>% 
			strsplit("") %>% 
			lapply(., function(x) paste(set_names[x==1], collapse=" & ")) %>% 
			unlist
	olap_exp <- mset$overlap.expected
	p <- mset$P.value
	msetdf <- tibble(Names=olap_names, Size=olap_sizes, Expected=olap_exp, pvalue=p) %>%
		  mutate(Names=fct_reorder(Names, desc(Size)))
	return(msetdf)
}

pdf(paste0(BASEDIR, "/beta-actin.pdf"), width=14, height=10)
ggarrange(cowplot::plot_grid(wt, kd, ncol=2, nrow=1, align="v", axis="lrt", labels=c("WT", "KD")))
ggplot(p, aes(x=p, fill=lab)) + geom_density(colour="black", alpha=0.4) + facet_wrap(~pos, ncol=1) + theme_bw(22)
upset(fromList(intersections), text.scale=1.3)

Result=supertest(intersections,n=n)
mset2df(Result) %>% ggplot(aes(x=Names, y=Size, fill=-log10(pvalue))) + 
	geom_col(colour="black") + 
	geom_col(aes(y=Expected), colour="red", linetype="dotted") + 
	geom_text(aes(y=Size+60, label=paste0(Size,"\n(",round(Size/n,3)*100,"%)"))) + 
	scale_fill_distiller(palette = "OrRd", direction=1, na.value="white") + 
	theme_bw(22) + 
	theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) + xlab("")
dev.off()


dat <- read_tsv(paste0(BASEDIR,"out_actin/out_data.tsv"), col_names=c("pos", "lab", "cluster", "prob", "intensity", "dwell")) %>% mutate(lab=gsub("_[12]", "", lab), pos=paste0("A",pos+2)) %>%
	mutate(pos=factor(pos, levels=c("A652", "A1324", "A1535"))) %>%
	mutate(lab=factor(lab, levels=c("WT", "KD"))) %>%
	group_by(pos) %>% 
	mutate(norm_int=(intensity-mean(intensity))/sd(intensity), norm_dwell=(dwell-mean(dwell))/sd(dwell))

pdf(paste0(BASEDIR, "/beta-actin-scatterplots.pdf"), width=15, height=7)
ggplot(dat, aes(x=norm_int, y=norm_dwell)) +  geom_point(aes(colour=prob), size=0.5, alpha=0.4) + 
	stat_density2d(geom="density2d", aes(alpha=..nlevel..), size=0.3, colour="black", contour=TRUE, show.legend = F) +
	scale_colour_distiller(palette = "Spectral", name="m6A probability") +
	facet_grid(lab~pos, scales="free") + theme_bw(22) + xlim(-3,3)  + ylim(-3,3) +
	xlab("Current intensity (z-score normalised)") +
	ylab("Dwell time (log10, z-score normalised)")
dev.off() 


mat <- filter(p, lab=="WT") %>% reshape2::dcast(read~pos) %>%
        select(-read) %>% filter(!is.na(A650), !is.na(A1322), !is.na(A1533)) %>% as_tibble()

matkd <- filter(p, lab=="KD") %>% reshape2::dcast(read~pos) %>%
        select(-read) %>% filter(!is.na(A650), !is.na(A1322), !is.na(A1533)) %>% as_tibble()

getprob <- function(x, mod_freqs){
                if(x[1]) {p1 <- mod_freqs[["A650"]]} else  {p1 <-(1 - mod_freqs[["A650"]])}
                if(x[2]) {p2 <- mod_freqs[["A1322"]]} else {p2 <-(1 - mod_freqs[["A1322"]])}
                if(x[3]) {p3 <- mod_freqs[["A1533"]]} else {p3 <-(1 - mod_freqs[["A1533"]])}
                return(p1*p2*p3)
}

calc_p <-function(thr, full=F){
        binpos <- (mat> thr) %>% as_tibble()
        obs_n <- binpos %>% group_by(A650, A1322, A1533) %>% summarise(n=n())
        n_tot <- sum(obs_n$n)
        freq <- apply(binpos, 2, mean)
        obs_n$exp <- round(n_tot * apply(obs_n, 1, function(x) getprob(x, freq)))
        obs_n <- mutate(obs_n, chi = ((n-exp)^2)/exp)
	if(!full) return(pchisq(sum(obs_n$chi), df=7, lower.tail=F))
	else return(obs_n)
}
calc_p_v=Vectorize(calc_p)

data.frame(thr=seq(0.1,1,0.05)) %>% mutate(p=calc_p_v(thr)) %>% mutate(adj=p.adjust(p, method="BH"))
calc_p(0.75, full=T)

pdf(paste0(BASEDIR, "/beta-actin-expected_freq.pdf"))
calc_p(0.75, full=T) %>% 
	ggplot(aes(x=paste(A650, A1322, A1533), y=n)) + 
		geom_col() + 
		geom_col(aes(y=exp), colour="blue", fill="blue", alpha=0.1 ,linetype="dotted") + 
		theme_bw(22) + ylab("Number of molecules") + 
		theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
dev.off()
