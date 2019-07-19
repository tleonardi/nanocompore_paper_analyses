#!/bin/bash
BASEDIR=$(git rev-parse --show-toplevel)
CONDA_ENVS="$BASEDIR/conda_environments/"
ANALYSIS="${BASEDIR}/profiles/METTL3_polyA/analysis"
RESULTS="${BASEDIR}/profiles/METTL3_polyA/results"
TRANSCRIPTOME="${BASEDIR}/profiles/METTL3_polyA/data/references/reference_transcriptome.bed"
RIP_PEAKS="${BASEDIR}/profiles/METTL3_polyA/data/rip/MOLM13_m6A_RIP/diff_peak.tsv"

mkdir -p $BASEDIR
mkdir -p $ANALYSIS/Ule_miCLIP
mkdir -p $ANALYSIS/rip_coverage

# Create the conda environment
if [[ ! -d "$CONDA_ENVS/homer" ]]; then
	conda env create -p "$CONDA_ENVS/homer" --file "$CONDA_ENVS/homer.yml"
	source activate "$CONDA_ENVS/homer"
else
	source activate "$CONDA_ENVS/homer"
fi

tail -n +2 ${BASEDIR}/profiles/METTL3_polyA/data/nanocompore/bed_files/sig_sites_GMM_pvalue_thr0.05.bed | bedparse convertChr --assembly hg38 --target ucsc > ${BASEDIR}/profiles/METTL3_polyA/data/nanocompore/bed_files/sig_sites_GMM_pvalue_thr0.05.ucsc.bed
tail -n +2 ${BASEDIR}/profiles/METTL3_polyA/data/nanocompore/bed_files/sig_sites_GMM_pvalue_context_2_thr0.05.bed | bedparse convertChr --assembly hg38 --target ucsc > ${BASEDIR}/profiles/METTL3_polyA/data/nanocompore/bed_files/sig_sites_GMM_pvalue_context_2_thr0.05.ucsc.bed
tail -n +2 ${BASEDIR}/profiles/METTL3_polyA/data/nanocompore/bed_files/sig_sites_KS_intensity_pvalue_thr0.05.bed | bedparse convertChr --assembly hg38 --target ucsc > ${BASEDIR}/profiles/METTL3_polyA/data/nanocompore/bed_files/sig_sites_KS_intensity_pvalue_thr0.05.ucsc.bed


         
# Filter RIP sites that overlap Nanocompore transcripts                                                                                                                                                                                                                                     
bedparse filter --annotation <(awk 'BEGIN{OFS=FS="\t"}NR>1{print $3}' ${BASEDIR}/profiles/METTL3_polyA/data/nanocompore/nanocompore_results.txt) -c 1 $TRANSCRIPTOME | bedparse convertChr --assembly hg38 --target ucsc > $ANALYSIS/rip_coverage/nanocompore_targeted_tx.bed      
bedparse bed12tobed6  $ANALYSIS/rip_coverage/nanocompore_targeted_tx.bed > $ANALYSIS/rip_coverage/nanocompore_targeted_tx.bed6                                                                                                                                                              
                                                                                                                                                                                                                                                                                            
cat <(head -n 1 $RIP_PEAKS) <(bedtools intersect -u -s -a $RIP_PEAKS -b $ANALYSIS/rip_coverage/nanocompore_targeted_tx.bed6) > $ANALYSIS/rip_coverage/rip_diff_peaks_nanocompore_tx.txt                                                                                                     
                                                                                                                                                                                                                                                                                            
# Filter METTL3 dep RIP sites that overlap Nanocompore transcripts                                                                                                                                                                                                                          
awk 'BEGIN{OFS=FS="\t"}$14<-1 && $16<-1 && $18<0' $RIP_PEAKS > $ANALYSIS/rip_coverage/rip_diff_peaks_METTL3dep.bed                                                                                                                                                                          
cat <(head -n 1 $RIP_PEAKS) <(bedtools intersect -u -s -a $ANALYSIS/rip_coverage/rip_diff_peaks_METTL3dep.bed -b $ANALYSIS/rip_coverage/nanocompore_targeted_tx.bed6) > $ANALYSIS/rip_coverage/rip_diff_peaks_METTL3dep_nanocompore_tx.txt     


# Filter miCLIP sites that overlap Nanocompore transcripts
THR=0
for f in  ${BASEDIR}/METTL3_polyA/data/Ule_miCLIP/*bedgraph; do
	NAME=$(basename $f .bedgraph)
	awk 'BEGIN{OFS=FS="\t"}/^chr/{if($4<0){strand="-"; $4=-$4}else{strand="+"}; print $1, $2-2, $3+2, "peak_"NR, $4, strand}'  $f | sort -k1,1 -k2,2n | bedtools merge -s -c 5,6 -o sum,distinct | awk -v thr=$THR 'BEGIN{OFS=FS="\t"}$4>thr{print $1,$2,$3,"Peak"NR, $4, $5}' > ${ANALYSIS}/Ule_miCLIP/${NAME}.bedgraph
	maxscore=$(cut -f5 ${ANALYSIS}/Ule_miCLIP/${NAME}.bedgraph | sort -n | tail -1)
	(echo "track name=$NAME useScore=1"
	awk -vmaxscore=$maxscore 'BEGIN{OFS=FS="\t"}{print $1,$2,$3,$4,int(($5/maxscore)*1000), $6}' ${ANALYSIS}/Ule_miCLIP/${NAME}.bedgraph )> ${ANALYSIS}/Ule_miCLIP/${NAME}.bedgraph.track
	mkdir -p  ${ANALYSIS}/Ule_miCLIP/strand_split/
	awk 'BEGIN{OFS=FS="\t"}$6=="+"{print $1,$2,$3,$5}' ${ANALYSIS}/Ule_miCLIP/${NAME}.bedgraph | sort -k1,1 -k2,2n > ${ANALYSIS}/Ule_miCLIP/strand_split/${NAME}.bedgraph.plus
	awk 'BEGIN{OFS=FS="\t"}$6=="-"{print $1,$2,$3,$5}' ${ANALYSIS}/Ule_miCLIP/${NAME}.bedgraph | sort -k1,1 -k2,2n > ${ANALYSIS}/Ule_miCLIP/strand_split/${NAME}.bedgraph.minus
done

PREF="${ANALYSIS}/Ule_miCLIP/strand_split/"
bedtools unionbedg -i ${PREF}/m6A_MOLM13_1.genome.bedgraph.plus ${PREF}/m6A_MOLM13_2.genome.bedgraph.plus ${PREF}/m6A_MOLM13_3.genome.bedgraph.plus ${PREF}/m6A_MOLM13_4.genome.bedgraph.plus ${PREF}/m6A_MOLM13_M3_KO_1.genome.bedgraph.plus ${PREF}/m6A_MOLM13_M3_KO_2.genome.bedgraph.plus > ${PREF}/plus.unionbedg
bedtools unionbedg -i ${PREF}/m6A_MOLM13_1.genome.bedgraph.minus ${PREF}/m6A_MOLM13_2.genome.bedgraph.minus ${PREF}/m6A_MOLM13_3.genome.bedgraph.minus ${PREF}/m6A_MOLM13_4.genome.bedgraph.minus ${PREF}/m6A_MOLM13_M3_KO_1.genome.bedgraph.minus ${PREF}/m6A_MOLM13_M3_KO_2.genome.bedgraph.minus > ${PREF}/minus.unionbedg

(
awk 'BEGIN{OFS=FS="\t"}{print $1,$2,$3,".",0,"+",$4,$5,$6,$7,$8,$9}' ${PREF}/plus.unionbedg
awk 'BEGIN{OFS=FS="\t"}{print $1,$2,$3,".",0,"-",$4,$5,$6,$7,$8,$9}' ${PREF}/minus.unionbedg
) | sort -k1,1 -k2,2n > ${PREF}/combined.unionbedg


awk 'BEGIN{OFS=FS="\t"}{
	meanKD=($11+$12)/2
	meanWT=($7+$8+$9+$10)/4
	logFC=log((meanKD+1)/(meanWT+1))/log(2)

	if((($7>0)+($8>0)+($9>0)+($10>0))>2 && ($7+$8+$9+$10)>10 && logFC<-1){
		print $1,$2,$3,$4,logFC,$6
	}
}' ${PREF}/combined.unionbedg > ${ANALYSIS}/Ule_miCLIP/sig_clusters.bed

bedtools intersect -u -s -a ${ANALYSIS}/Ule_miCLIP/sig_clusters.bed -b $ANALYSIS/rip_coverage/nanocompore_targeted_tx.bed6 > ${ANALYSIS}/Ule_miCLIP/sig_clusters_target_olap.bed


mkdir ${ANALYSIS}/rip_clip_olaps
# CLIP 
bedtools intersect -u -s -split -a ${ANALYSIS}/Ule_miCLIP/sig_clusters_target_olap.bed  -b $ANALYSIS/rip_coverage/rip_diff_peaks_nanocompore_tx.txt  > ${ANALYSIS}/rip_clip_olaps/clip_w_rip.bed

bedtools intersect -u -s -split -a ${ANALYSIS}/Ule_miCLIP/sig_clusters_target_olap.bed  -b ${BASEDIR}/profiles/METTL3_polyA/data/nanocompore/bed_files/sig_sites_GMM_pvalue_thr0.05.ucsc.bed  > ${ANALYSIS}/rip_clip_olaps/clip_w_nanocompore.bed

bedtools intersect -u -s -split -a ${ANALYSIS}/rip_clip_olaps/clip_w_rip.bed  -b ${BASEDIR}/profiles/METTL3_polyA/data/nanocompore/bed_files/sig_sites_GMM_pvalue_thr0.05.ucsc.bed   > ${ANALYSIS}/rip_clip_olaps/clip_w_rip_w_nanocompore.bed

# RIP
bedtools intersect -u -s -split -a $ANALYSIS/rip_coverage/rip_diff_peaks_nanocompore_tx.txt  -b ${ANALYSIS}/Ule_miCLIP/sig_clusters_target_olap.bed  > ${ANALYSIS}/rip_clip_olaps/rip_w_clip.bed

bedtools intersect -u -s -split -a $ANALYSIS/rip_coverage/rip_diff_peaks_nanocompore_tx.txt  -b ${BASEDIR}/profiles/METTL3_polyA/data/nanocompore/bed_files/sig_sites_GMM_pvalue_thr0.05.ucsc.bed  > ${ANALYSIS}/rip_clip_olaps/rip_w_nanocompore.bed

bedtools intersect -u -s -split -a ${ANALYSIS}/rip_clip_olaps/rip_w_clip.bed  -b ${BASEDIR}/profiles/METTL3_polyA/data/nanocompore/bed_files/sig_sites_GMM_pvalue_thr0.05.ucsc.bed   > ${ANALYSIS}/rip_clip_olaps/rip_w_clip_w_nanocompore.bed


# RIP METTL3 dep
bedtools intersect -u -s -split -a $ANALYSIS/rip_coverage/rip_diff_peaks_METTL3dep_nanocompore_tx.txt  -b ${ANALYSIS}/Ule_miCLIP/sig_clusters_target_olap.bed > ${ANALYSIS}/rip_clip_olaps/rip_METTL3dep_w_clip.bed

bedtools intersect -u -s -split -a $ANALYSIS/rip_coverage/rip_diff_peaks_METTL3dep_nanocompore_tx.txt  -b ${BASEDIR}/profiles/METTL3_polyA/data/nanocompore/bed_files/sig_sites_GMM_pvalue_thr0.05.ucsc.bed  > ${ANALYSIS}/rip_clip_olaps/rip_METTL3dep_w_nanocompore.bed

bedtools intersect -u -s -split -a ${ANALYSIS}/rip_clip_olaps/rip_METTL3dep_w_clip.bed  -b ${BASEDIR}/profiles/METTL3_polyA/data/nanocompore/bed_files/sig_sites_GMM_pvalue_thr0.05.ucsc.bed   > ${ANALYSIS}/rip_clip_olaps/rip_METTL3dep_w_clip_w_nanocompore.bed



(cat <<EOF
Type	Validation	Total	w_nanocompore
RIP	RIP_only	$(wc -l $ANALYSIS/rip_coverage/rip_diff_peaks_nanocompore_tx.txt|cut -d ' ' -f1)	$(wc -l ${ANALYSIS}/rip_clip_olaps/rip_w_nanocompore.bed|cut -d ' ' -f1)
RIP	With_CLIP	$(wc -l ${ANALYSIS}/rip_clip_olaps/rip_w_clip.bed|cut -d ' ' -f1)	$(wc -l ${ANALYSIS}/rip_clip_olaps/rip_w_clip_w_nanocompore.bed|cut -d ' ' -f1)
RIP_M3	RIP_only	$(wc -l $ANALYSIS/rip_coverage/rip_diff_peaks_METTL3dep_nanocompore_tx.txt|cut -d ' ' -f1)	$(wc -l ${ANALYSIS}/rip_clip_olaps/rip_METTL3dep_w_nanocompore.bed|cut -d ' ' -f1)
RIP_M3	With_CLIP	$(wc -l ${ANALYSIS}/rip_clip_olaps/rip_METTL3dep_w_clip.bed|cut -d ' ' -f1)	$(wc -l ${ANALYSIS}/rip_clip_olaps/rip_METTL3dep_w_clip_w_nanocompore.bed|cut -d ' ' -f1)
CLIP	CLIP_only	$(wc -l ${ANALYSIS}/Ule_miCLIP/sig_clusters_target_olap.bed|cut -d ' ' -f1)	$(wc -l ${ANALYSIS}/rip_clip_olaps/clip_w_nanocompore.bed|cut -d ' ' -f1)
CLIP	With_RIP	$(wc -l ${ANALYSIS}/rip_clip_olaps/clip_w_rip.bed|cut -d ' ' -f1)	$(wc -l ${ANALYSIS}/rip_clip_olaps/clip_w_rip_w_nanocompore.bed|cut -d ' ' -f1)
EOF
)| perl -pe 's/ +/\t/g' >${ANALYSIS}/rip_clip_olaps/rip_clip_counts.txt


bedtools merge -s -c 6 -o distinct -i <(sort -k1,1 -k2,2n ${BASEDIR}/profiles/METTL3_polyA/data/nanocompore/bed_files/sig_sites_GMM_pvalue_thr0.05.ucsc.bed) | awk 'BEGIN{OFS=FS="\t"}{print $1,$2,$3,".",0,$4}' | sort -k1,1 -k2,2n > ${BASEDIR}/profiles/METTL3_polyA/data/nanocompore/bed_files/sig_sites_GMM_pvalue_thr0.05.ucsc.merged.bed

bedtools intersect -u -s -split  -a ${BASEDIR}/profiles/METTL3_polyA/data/nanocompore/bed_files/sig_sites_GMM_pvalue_thr0.05.ucsc.merged.bed   -b ${ANALYSIS}/Ule_miCLIP/sig_clusters_target_olap.bed  > ${ANALYSIS}/rip_clip_olaps/nanocompore_with_clip.bed

bedtools intersect -u -s -split  -a ${BASEDIR}/profiles/METTL3_polyA/data/nanocompore/bed_files/sig_sites_GMM_pvalue_thr0.05.ucsc.merged.bed   -b $ANALYSIS/rip_coverage/rip_diff_peaks_nanocompore_tx.txt  > ${ANALYSIS}/rip_clip_olaps/nanocompore_with_rip.bed

bedtools intersect -u -s -split  -a ${ANALYSIS}/rip_clip_olaps/nanocompore_with_clip.bed  -b ${ANALYSIS}/rip_clip_olaps/nanocompore_with_rip.bed  > ${ANALYSIS}/rip_clip_olaps/nanocompore_with_clip_and_rip.bed 

(cat <<EOF
Type	N
Total	$(wc -l ${BASEDIR}/profiles/METTL3_polyA/data/nanocompore/bed_files/sig_sites_GMM_pvalue_thr0.05.ucsc.merged.bed | cut -d ' ' -f1)
With_CLIP	$(wc -l ${ANALYSIS}/rip_clip_olaps/nanocompore_with_clip.bed | cut -d ' ' -f1)
With_RIP	$(wc -l ${ANALYSIS}/rip_clip_olaps/nanocompore_with_rip.bed | cut -d ' ' -f1)
With_Both	$(wc -l ${ANALYSIS}/rip_clip_olaps/nanocompore_with_clip_and_rip.bed | cut -d ' ' -f1)
EOF
) | perl -pe 's/ +/\t/g' > ${ANALYSIS}/rip_clip_olaps/nanocompore_counts.txt


