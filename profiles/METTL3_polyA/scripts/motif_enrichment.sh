#!/bin/bash
BASEDIR=$(git rev-parse --show-toplevel)
CONDA_ENVS="$BASEDIR/conda_environments/"
ANALYSIS="${BASEDIR}/profiles/METTL3_polyA/analysis"
RESULTS="${BASEDIR}/profiles/METTL3_polyA/results"
BED="${BASEDIR}/nanocompore_pipelines/METTL3_KD_polyA/results/nanocompore/out_sig_sites_GMM_logit_pvalue_thr_0.05.bed"
TRANSCRIPTOME="${BASEDIR}/nanocompore_pipelines/METTL3_KD_polyA/results/references/reference_transcriptome.bed"
NANOCOMPORE_FULL="${BASEDIR}/nanocompore_pipelines/METTL3_KD_polyA/results/nanocompore/out_nanocompore_results.tsv"
GENOME="$HOME/tmp/Homo_sapiens.GRCh38.dna.toplevel.fa"
if [[ ! -f $GENOME ]]; then
	wget -O - "ftp://ftp.ensembl.org/pub/release-96/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.toplevel.fa.gz" | gunzip > $GENOME
fi

mkdir -p $BASEDIR
# Create the conda environment
if [[ ! -d "$CONDA_ENVS/homer" ]]; then
	conda env create -p "$CONDA_ENVS/homer" --file "$CONDA_ENVS/homer.yml"
	conda activate "$CONDA_ENVS/homer"
else
	conda activate "$CONDA_ENVS/homer"
fi

mkdir -p "$ANALYSIS/motifs"


# Extract the top 100 sites

awk 'BEGIN{OFS=FS="\t"}NR>1{print $1,$2-5, $3+5, $4, $5, $6}' $BED | sort -k5,5nr | head -n 100 > $ANALYSIS/motifs/sig_sites_ext.bed

# Remove duplicated genomic coordinates to avoid
# counting twice sites that are a significant hit
# for two isoforms of the same gene.
sort -k1,1 -k2,2n $ANALYSIS/motifs/sig_sites_ext.bed | bedtools merge -s -c 5,6 -o max,distinct -i - | awk 'BEGIN{OFS=FS="\t"}NR>1{print $1,$2, $3, ".", $4, $5}' | sort -k 5,5nr > $ANALYSIS/motifs/sig_sites_ext_unique.bed
bedtools getfasta -s -fullHeader -fi $GENOME -bed $ANALYSIS/motifs/sig_sites_ext_unique.bed > $ANALYSIS/motifs/sig_sites_ext_unique.fa 

# Use nanocompore trancripts as background
tail -n +2 $BED | cut -f4 | perl -pe 's/_.+//' | sort -u > $ANALYSIS/motifs/bg_transcripts.txt

# Get the a BED file of BG transcripts
bedparse filter --annotation $ANALYSIS/motifs/bg_transcripts.txt $TRANSCRIPTOME > $ANALYSIS/motifs/bg.bed

# Convert BG transcripts to BED6 and subtract significant sites
bedparse bed12tobed6 --appendExN $ANALYSIS/motifs/bg.bed | sort -k1,1 -k2,2n > $ANALYSIS/motifs/bg.bed6
bedtools subtract -s -a $ANALYSIS/motifs/bg.bed6 -b $ANALYSIS/motifs/sig_sites_ext_unique.bed | sort -k1,1 -k2,2n | bedtools merge -s -c 6 -o distinct -i - | awk 'BEGIN{OFS=FS="\t"}NR>1{print $1,$2, $3, ".", ".", $4}' > $ANALYSIS/motifs/bg_noOlap.bed

# Get background Fasta
bedtools getfasta -s -name -fi $GENOME -bed $ANALYSIS/motifs/bg_noOlap.bed > $ANALYSIS/motifs/bg_noOlap.fa

# Run Homer
rm -rf $ANALYSIS/motifs/homer
findMotifs.pl $ANALYSIS/motifs/sig_sites_ext_unique.fa fasta $ANALYSIS/motifs/homer -fasta $ANALYSIS/motifs/bg_noOlap.fa -p 5 -len 5 -S 5 -norevopp -rna -basic

#cp -R $ANALYSIS/motifs/homer $RESULTS

