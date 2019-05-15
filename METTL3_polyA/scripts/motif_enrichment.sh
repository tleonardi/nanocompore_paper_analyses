#!/bin/bash
BASEDIR=$(git rev-parse --show-toplevel)
CONDA_ENVS="$BASEDIR/conda_environments/"
ANALYSIS="${BASEDIR}/METTL3_polyA/analysis"
BED="${BASEDIR}/METTL3_polyA/data/nanocompore/bed_files/sig_sites_GMM_pvalue_thr0.05.bed"
TRANSCRIPTOME="${BASEDIR}/METTL3_polyA/data/references/reference_transcriptome.bed"

GENOME="$HOME/tmp/Homo_sapiens.GRCh38.dna.toplevel.fa"
if [[ ! -f $GENOME ]]; then
	wget -O - "ftp://ftp.ensembl.org/pub/release-96/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.toplevel.fa.gz" | gunzip > $GENOME
fi

mkdir -p $BASEDIR
# Create the conda environment
if [[ ! -d "$CONDA_ENVS/homer" ]]; then
	conda env create -p "$CONDA_ENVS/homer" --file "$CONDA_ENVS/homer.yml"
	source activate "$CONDA_ENVS/homer"
else
	source activate "$CONDA_ENVS/homer"
fi

mkdir -p "$ANALYSIS/motifs"


# Extract the top 100 sites
awk 'BEGIN{OFS=FS="\t"}NR>1{print $1,$2-5, $3+5, $4, $5, $6}' $BED | sort -k5,5nr | head -n 100 > $ANALYSIS/motifs/sig_sites_ext.bed

# Remove duplicated genomic coordinates to avoid
# counting twice sites that are a significant hit
# for two isoforms of the same gene.
awk 'BEGIN{OFS=FS="\t"}NR>1{print $1,$2, $3, ".", $5, $6}' $ANALYSIS/motifs/sig_sites_ext.bed > $ANALYSIS/motifs/sig_sites_ext_unique.bed
bedtools getfasta -s -name -fi $GENOME -bed $ANALYSIS/motifs/sig_sites_ext_unique.bed > $ANALYSIS/motifs/sig_sites_ext_unique.fa 

# Use the very same trancripts as background
cut -f4 $ANALYSIS/motifs/sig_sites_ext.bed | perl -pe 's/_.+//' | sort -u > $ANALYSIS/motifs/bg_transcripts.bed

# Get the a BED file of BG transcripts
bedparse filter --annotation $ANALYSIS/motifs/bg_transcripts.bed $TRANSCRIPTOME > $ANALYSIS/motifs/bg.bed

# Convert BG transcripts to BED6 and subtract significant sites
bedparse bed12tobed6 --appendExN $ANALYSIS/motifs/bg.bed > $ANALYSIS/motifs/bg.bed6
bedtools subtract -s -a $ANALYSIS/motifs/bg.bed6 -b $ANALYSIS/motifs/sig_sites_ext.bed > $ANALYSIS/motifs/bg_noOlap.bed

# Get background Fasta
bedtools getfasta -s -split -name -fi $GENOME -bed $ANALYSIS/motifs/bg_noOlap.bed > $ANALYSIS/motifs/bg_noOlap.fa

# Run Homer
findMotifs.pl $ANALYSIS/motifs/sig_sites_ext_unique.fa fasta $ANALYSIS/motifs/homer -fasta $ANALYSIS/motifs/bg_noOlap.fa -p 5 -len 5 -S 5 -norevopp 

