#!/bin/bash
BASEDIR=$(git rev-parse --show-toplevel)
CONDA_ENVS="$BASEDIR/conda_environments/"
ANALYSIS="${BASEDIR}/METTL3_polyA/analysis"
RESULTS="${BASEDIR}/METTL3_polyA/results"
TRANSCRIPTOME="${BASEDIR}/METTL3_polyA/data/references/reference_transcriptome.bed"
RIP_PEAKS="${BASEDIR}/METTL3_polyA/data/MOLM13_m6A_RIP/diff_peak.tsv"

mkdir -p $BASEDIR
mkdir -p $ANALYSIS/rip_coverage

# Create the conda environment
if [[ ! -d "$CONDA_ENVS/homer" ]]; then
	conda env create -p "$CONDA_ENVS/homer" --file "$CONDA_ENVS/homer.yml"
	source activate "$CONDA_ENVS/homer"
else
	source activate "$CONDA_ENVS/homer"
fi

tail -n +2 ${BASEDIR}/METTL3_polyA/data/nanocompore/bed_files/sig_sites_GMM_pvalue_thr0.05.bed | bedparse convertChr --assembly hg38 --target ucsc > ${BASEDIR}/METTL3_polyA/data/nanocompore/bed_files/sig_sites_GMM_pvalue_thr0.05.ucsc.bed
tail -n +2 ${BASEDIR}/METTL3_polyA/data/nanocompore/bed_files/sig_sites_GMM_pvalue_context_2_thr0.05.bed | bedparse convertChr --assembly hg38 --target ucsc > ${BASEDIR}/METTL3_polyA/data/nanocompore/bed_files/sig_sites_GMM_pvalue_context_2_thr0.05.ucsc.bed
tail -n +2 ${BASEDIR}/METTL3_polyA/data/nanocompore/bed_files/sig_sites_KS_intensity_pvalue_thr0.05.bed | bedparse convertChr --assembly hg38 --target ucsc > ${BASEDIR}/METTL3_polyA/data/nanocompore/bed_files/sig_sites_KS_intensity_pvalue_thr0.05.ucsc.bed


# Filter RIP sites that overlap Nanocompore transcripts
bedparse filter --annotation <(awk 'BEGIN{OFS=FS="\t"}NR>1{gsub(/_.+$/, "", $4); print $4}' ${BASEDIR}/METTL3_polyA/data/nanocompore/GMM_pvalue_thr1.bed) -c 1 $TRANSCRIPTOME | bedparse convertChr --assembly hg38 --target ucsc > $ANALYSIS/rip_coverage/nanocompore_targeted_tx.bed
bedparse bed12tobed6  $ANALYSIS/rip_coverage/nanocompore_targeted_tx.bed > $ANALYSIS/rip_coverage/nanocompore_targeted_tx.bed6

cat <(head -n 1 $RIP_PEAKS) <(bedtools intersect -u -s -a $RIP_PEAKS -b $ANALYSIS/rip_coverage/nanocompore_targeted_tx.bed6) > $ANALYSIS/rip_coverage/rip_diff_peaks_nanocompore_tx.txt

# Filter METTL3 dep RIP sites that overlap Nanocompore transcripts



# Filter miCLIP sites that overlap Nanocompore transcripts
mkdir -p ${ANALYSIS}/Vu_m6A_miCLIP/
awk 'BEGIN{OFS=FS="\t"}{print $1, $2-2, $3+2, $4, $5, $6}'  ${BASEDIR}/METTL3_polyA/data/Vu_m6A_miCLIP/Vu_m6A_miCLIP_Rep1.bed > ${ANALYSIS}/Vu_m6A_miCLIP/Vu_m6A_miCLIP_Rep1_5nt.bed
awk 'BEGIN{OFS=FS="\t"}{print $1, $2-2, $3+2, $4, $5, $6}'  ${BASEDIR}/METTL3_polyA/data/Vu_m6A_miCLIP/Vu_m6A_miCLIP_Rep2.bed > ${ANALYSIS}/Vu_m6A_miCLIP/Vu_m6A_miCLIP_Rep2_5nt.bed
awk 'BEGIN{OFS=FS="\t"}{print $1, $2-2, $3+2, $4, $5, $6}'  ${BASEDIR}/METTL3_polyA/data/Vu_m6A_miCLIP/Vu_m6A_miCLIP_Rep3.bed > ${ANALYSIS}/Vu_m6A_miCLIP/Vu_m6A_miCLIP_Rep3_5nt.bed

bedtools intersect -s -u -a ${ANALYSIS}/Vu_m6A_miCLIP/Vu_m6A_miCLIP_Rep1_5nt.bed -b ${ANALYSIS}/Vu_m6A_miCLIP/Vu_m6A_miCLIP_Rep2_5nt.bed > ${ANALYSIS}/Vu_m6A_miCLIP/Vu_m6A_miCLIP_Rep1_Rep2_5nt.bed
bedtools intersect -s -u -a ${ANALYSIS}/Vu_m6A_miCLIP/Vu_m6A_miCLIP_Rep3_5nt.bed -b ${ANALYSIS}/Vu_m6A_miCLIP/Vu_m6A_miCLIP_Rep1_Rep2_5nt.bed > ${ANALYSIS}/Vu_m6A_miCLIP/Vu_m6A_miCLIP_allReps_5nt.bed
bedtools intersect -u -s -a ${ANALYSIS}/Vu_m6A_miCLIP/Vu_m6A_miCLIP_allReps_5nt.bed -b $ANALYSIS/rip_coverage/nanocompore_targeted_tx.bed6 > ${ANALYSIS}/Vu_m6A_miCLIP/Vu_m6A_miCLIP_allReps_targeted_5nt.bed

cat ${ANALYSIS}/Vu_m6A_miCLIP/Vu_m6A_miCLIP_Rep1_5nt.bed ${ANALYSIS}/Vu_m6A_miCLIP/Vu_m6A_miCLIP_Rep2_5nt.bed ${ANALYSIS}/Vu_m6A_miCLIP/Vu_m6A_miCLIP_Rep3_5nt.bed | sort -k1,1 -k2,2n | bedtools merge -s -c 6 -o distinct | awk 'BEGIN{OFS=FS="\t"}{print $1, $2, $3, ".", ".", $4}' > ${ANALYSIS}/Vu_m6A_miCLIP/Vu_m6A_miCLIP_merged_rep.bed

bedtools intersect -u -s -a ${ANALYSIS}/Vu_m6A_miCLIP/Vu_m6A_miCLIP_merged_rep.bed -b $ANALYSIS/rip_coverage/nanocompore_targeted_tx.bed6 > ${ANALYSIS}/Vu_m6A_miCLIP/Vu_m6A_miCLIP_merged_rep_target_olap.bed

cat ${BASEDIR}/METTL3_polyA/data/Vu_m6A_miCLIP/Vu_m6A_miCLIP_Rep1.bed ${BASEDIR}/METTL3_polyA/data/Vu_m6A_miCLIP/Vu_m6A_miCLIP_Rep2.bed ${BASEDIR}/METTL3_polyA/data/Vu_m6A_miCLIP/Vu_m6A_miCLIP_Rep3.bed | sort -k1,1 -k2,2n | bedtools merge -s -c 6 -o distinct | awk 'BEGIN{OFS=FS="\t"}{print $1, $2, $3, ".", ".", $4}' > ${ANALYSIS}/Vu_m6A_miCLIP/Vu_m6A_miCLIP_merged_rep_noExt.bed


mkdir ${ANALYSIS}/rip_clip_olaps
# CLIP 
bedtools intersect -u -s -split \
	-a ${ANALYSIS}/Vu_m6A_miCLIP/Vu_m6A_miCLIP_merged_rep_target_olap.bed \
	-b $ANALYSIS/rip_coverage/rip_diff_peaks_nanocompore_tx.txt \
	> ${ANALYSIS}/rip_clip_olaps/clip_w_rip.bed

bedtools intersect -u -s -split \
	-a ${ANALYSIS}/Vu_m6A_miCLIP/Vu_m6A_miCLIP_merged_rep_target_olap.bed \
	-b ${BASEDIR}/METTL3_polyA/data/nanocompore/bed_files/sig_sites_GMM_pvalue_thr0.05.ucsc.bed \
	> ${ANALYSIS}/rip_clip_olaps/clip_w_nanocompore.bed

bedtools intersect -u -s -split \
	-a ${ANALYSIS}/rip_clip_olaps/clip_w_rip.bed \
	-b ${BASEDIR}/METTL3_polyA/data/nanocompore/bed_files/sig_sites_GMM_pvalue_thr0.05.ucsc.bed  \
	> ${ANALYSIS}/rip_clip_olaps/clip_w_rip_w_nanocompore.bed

# RIP
bedtools intersect -u -s -split \
	-a $ANALYSIS/rip_coverage/rip_diff_peaks_nanocompore_tx.txt \
	-b ${ANALYSIS}/Vu_m6A_miCLIP/Vu_m6A_miCLIP_merged_rep_target_olap.bed \
	> ${ANALYSIS}/rip_clip_olaps/rip_w_clip.bed

bedtools intersect -u -s -split \
	-a $ANALYSIS/rip_coverage/rip_diff_peaks_nanocompore_tx.txt \
	-b ${BASEDIR}/METTL3_polyA/data/nanocompore/bed_files/sig_sites_GMM_pvalue_thr0.05.ucsc.bed \
	> ${ANALYSIS}/rip_clip_olaps/rip_w_nanocompore.bed

bedtools intersect -u -s -split \
	-a ${ANALYSIS}/rip_clip_olaps/rip_w_clip.bed \
	-b ${BASEDIR}/METTL3_polyA/data/nanocompore/bed_files/sig_sites_GMM_pvalue_thr0.05.ucsc.bed  \
	> ${ANALYSIS}/rip_clip_olaps/rip_w_clip_w_nanocompore.bed


(cat <<EOF
Type	Validation	Total	w_nanocompore
RIP	RIP_only	$(wc -l $ANALYSIS/rip_coverage/rip_diff_peaks_nanocompore_tx.txt|cut -d ' ' -f1)	$(wc -l ${ANALYSIS}/rip_clip_olaps/rip_w_nanocompore.bed|cut -d ' ' -f1)
RIP	With_CLIP	$(wc -l ${ANALYSIS}/rip_clip_olaps/rip_w_clip.bed|cut -d ' ' -f1)	$(wc -l ${ANALYSIS}/rip_clip_olaps/rip_w_clip_w_nanocompore.bed|cut -d ' ' -f1)
CLIP	CLIP_only	$(wc -l ${ANALYSIS}/Vu_m6A_miCLIP/Vu_m6A_miCLIP_merged_rep_target_olap.bed|cut -d ' ' -f1)	$(wc -l ${ANALYSIS}/rip_clip_olaps/clip_w_nanocompore.bed|cut -d ' ' -f1)
CLIP	With_RIP	$(wc -l ${ANALYSIS}/rip_clip_olaps/clip_w_rip.bed|cut -d ' ' -f1)	$(wc -l ${ANALYSIS}/rip_clip_olaps/clip_w_rip_w_nanocompore.bed|cut -d ' ' -f1)
EOF
)| perl -pe 's/ +/\t/g' >${ANALYSIS}/rip_clip_olaps/rip_clip_counts.txt



bedtools intersect -u -s -split \
	-a ${BASEDIR}/METTL3_polyA/data/nanocompore/bed_files/sig_sites_GMM_pvalue_thr0.05.ucsc.bed  \
	-b ${ANALYSIS}/Vu_m6A_miCLIP/Vu_m6A_miCLIP_merged_rep_target_olap.bed \
	> ${ANALYSIS}/rip_clip_olaps/nanocompore_with_clip.bed

bedtools intersect -u -s -split \
	-a ${BASEDIR}/METTL3_polyA/data/nanocompore/bed_files/sig_sites_GMM_pvalue_thr0.05.ucsc.bed  \
	-b $ANALYSIS/rip_coverage/rip_diff_peaks_nanocompore_tx.txt \
	> ${ANALYSIS}/rip_clip_olaps/nanocompore_with_rip.bed

bedtools intersect -u -s -split \
	-a ${ANALYSIS}/rip_clip_olaps/nanocompore_with_clip.bed \
	-b ${ANALYSIS}/rip_clip_olaps/nanocompore_with_rip.bed \
	> ${ANALYSIS}/rip_clip_olaps/nanocompore_with_clip_and_rip.bed


(cat <<EOF
Type	N
Total	$(wc -l ${BASEDIR}/METTL3_polyA/data/nanocompore/bed_files/sig_sites_GMM_pvalue_thr0.05.ucsc.bed | cut -d ' ' -f1)
With_CLIP	$(wc -l ${ANALYSIS}/rip_clip_olaps/nanocompore_with_clip.bed | cut -d ' ' -f1)
With_RIP	$(wc -l ${ANALYSIS}/rip_clip_olaps/nanocompore_with_rip.bed | cut -d ' ' -f1)
With_Both	$(wc -l ${ANALYSIS}/rip_clip_olaps/nanocompore_with_clip_and_rip.bed | cut -d ' ' -f1)
EOF
) | perl -pe 's/ +/\t/g' > ${ANALYSIS}/rip_clip_olaps/nanocompore_counts.txt


