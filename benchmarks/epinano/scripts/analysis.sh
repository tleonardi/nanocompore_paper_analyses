#!/bin/bash
BASEDIR=$(git rev-parse --show-toplevel)
ANALYSIS="${BASEDIR}/benchmarks/epinano/analysis/"
NANOCOMPORE_FULL="${BASEDIR}/nanocompore_pipelines/METTL3_KD_polyA/results/nanocompore/out_nanocompore_results.tsv"
REF_TRANSCRIPTOME="${BASEDIR}/nanocompore_pipelines/METTL3_KD_polyA/results/references/reference_transcriptome.bed"
REF_TRANSCRIPTOME_FA="${BASEDIR}/nanocompore_pipelines/METTL3_KD_polyA/results/references/reference_transcriptome.fa"
GENOME="/home/tleonardi/nobackup/Homo_sapiens.GRCh38.dna.toplevel.fa"
TRANSCRIPTOME="/home/tleonardi/nobackup/Homo_sapiens.GRCh38.94.chr.bed"
SING="$BASEDIR/nanocompore_pipelines/singularity_images/nanocompore_4443e07.img"

source /hps/nobackup/enright/tom/miniconda3/etc/profile.d/conda.sh

if [[ -d $BASEDIR/conda_environments/epinano ]]; then
	conda activate $BASEDIR/conda_environments/epinano
else
	conda env create -f $BASEDIR/conda_environments/epinano.yml -p $BASEDIR/conda_environments/epinano
fi

mkdir -p $ANALYSIS/data/

#wget -O $ANALYSIS/data/unmod_rep1.tgz https://sra-pub-src-1.s3.amazonaws.com/SRR8767350/RNAAB089716.fast5.tar.gz.4
#wget -O $ANALYSIS/data/m6a_rep1.tgz https://sra-pub-src-1.s3.amazonaws.com/SRR8767351/RNAAB090763.fast5.tar.gz.3
#wget -O $ANALYSIS/data/unmod_rep2.tgz https://sra-pub-src-1.s3.amazonaws.com/SRR8767348/RNA081120181.fast5.tar.gz.2
#wget -O $ANALYSIS/data/m6a_rep2.tgz https://sra-pub-src-1.s3.amazonaws.com/SRR8767349/RNA081120182.fast5.tar.gz.2
#
#
#samp="m6a_rep1"
#
#
#
mkdir -p $ANALYSIS/METTL3_KD/
ln -s $REF_TRANSCRIPTOME_FA $ANALYSIS/METTL3_KD/reference_transcriptome.fa
ln -s $REF_TRANSCRIPTOME_FA.fai $ANALYSIS/METTL3_KD/reference_transcriptome.fa.fai
picard CreateSequenceDictionary R=$ANALYSIS/METTL3_KD/reference_transcriptome.fa

# Albacore basecalling (v2.1.7)
bsub -q standard -oo $BASEDIR/benchmarks/epinano/logs/%J.sout.log -e $BASEDIR/benchmarks/epinano/logs/%J.serr.log -M 20000 -n 10 -R 'rusage[mem=20000]'  "source /nfs/software/enright/virtualenvs/Albacore_2.1.7/bin/activate && read_fast5_basecaller.py -r -i /nfs/leia/research/enright/nanopore/datasets/RNA_RUN_04_METTL3_KD/20180522_1533_69/ -t 12 -s $ANALYSIS/METTL3_KD/albacore/KD_1 -f "FLO-MIN106" -k "SQK-RNA001" -o fast5 -q 0 --disable_pings --disable_filtering"
bsub -q standard -oo $BASEDIR/benchmarks/epinano/logs/%J.sout.log -e $BASEDIR/benchmarks/epinano/logs/%J.serr.log -M 20000 -n 10 -R 'rusage[mem=20000]'  "source /nfs/software/enright/virtualenvs/Albacore_2.1.7/bin/activate && read_fast5_basecaller.py -r -i /nfs/leia/research/enright/nanopore/datasets/RNA_RUN_06_METTL3_KD_2/20181206_1707_KD_2/ -t 12 -s $ANALYSIS/METTL3_KD/albacore/KD_2 -f "FLO-MIN106" -k "SQK-RNA001" -o fast5 -q 0 --disable_pings --disable_filtering"
bsub -q standard -oo $BASEDIR/benchmarks/epinano/logs/%J.sout.log -e $BASEDIR/benchmarks/epinano/logs/%J.serr.log -M 20000 -n 10 -R 'rusage[mem=20000]'  "source /nfs/software/enright/virtualenvs/Albacore_2.1.7/bin/activate && read_fast5_basecaller.py -r -i /nfs/leia/research/enright/nanopore/datasets/RNA_RUN_04_METTL3_KD/20180523_1703_CT_2/ -t 12 -s $ANALYSIS/METTL3_KD/albacore/WT_1 -f "FLO-MIN106" -k "SQK-RNA001" -o fast5 -q 0 --disable_pings --disable_filtering"
bsub -q standard -oo $BASEDIR/benchmarks/epinano/logs/%J.sout.log -e $BASEDIR/benchmarks/epinano/logs/%J.serr.log -M 20000 -n 10 -R 'rusage[mem=20000]'  "source /nfs/software/enright/virtualenvs/Albacore_2.1.7/bin/activate && read_fast5_basecaller.py -r -i /nfs/leia/research/enright/nanopore/datasets/RNA_RUN_06_METTL3_KD_2/20181206_1707_WT_2/ -t 12 -s $ANALYSIS/METTL3_KD/albacore/WT_2 -f "FLO-MIN106" -k "SQK-RNA001" -o fast5 -q 0 --disable_pings --disable_filtering"


bsub -q standard -oo $BASEDIR/benchmarks/epinano/logs/%J.sout.log -e $BASEDIR/benchmarks/epinano/logs/%J.serr.log -M 20000 -n 10 -R 'rusage[mem=20000]' "singularity exec $SING nanopolish extract -r -q -o $ANALYSIS/METTL3_KD/albacore/KD_1/KD_1.fastq $ANALYSIS/METTL3_KD/albacore/KD_1/workspace/"
bsub -q standard -oo $BASEDIR/benchmarks/epinano/logs/%J.sout.log -e $BASEDIR/benchmarks/epinano/logs/%J.serr.log -M 20000 -n 10 -R 'rusage[mem=20000]' "singularity exec $SING nanopolish extract -r -q -o $ANALYSIS/METTL3_KD/albacore/KD_2/KD_2.fastq $ANALYSIS/METTL3_KD/albacore/KD_2/workspace/"
bsub -q standard -oo $BASEDIR/benchmarks/epinano/logs/%J.sout.log -e $BASEDIR/benchmarks/epinano/logs/%J.serr.log -M 20000 -n 10 -R 'rusage[mem=20000]' "singularity exec $SING nanopolish extract -r -q -o $ANALYSIS/METTL3_KD/albacore/WT_1/WT_1.fastq $ANALYSIS/METTL3_KD/albacore/WT_1/workspace/"
bsub -q standard -oo $BASEDIR/benchmarks/epinano/logs/%J.sout.log -e $BASEDIR/benchmarks/epinano/logs/%J.serr.log -M 20000 -n 10 -R 'rusage[mem=20000]' "singularity exec $SING nanopolish extract -r -q -o $ANALYSIS/METTL3_KD/albacore/WT_2/WT_2.fastq $ANALYSIS/METTL3_KD/albacore/WT_2/workspace/"


# Mapping
mkdir -p $ANALYSIS/METTL3_KD/minimap/

bsub -q standard -oo $BASEDIR/benchmarks/epinano/logs/minimap_%J.sout.log -e $BASEDIR/benchmarks/epinano/logs/minimap_%J.serr.log -M 20000 -n 10 -R 'rusage[mem=20000]' "singularity exec $SING minimap2 -ax map-ont -t 10 -uf -k14 --for-only $ANALYSIS/METTL3_KD/reference_transcriptome.fa $ANALYSIS/METTL3_KD/albacore/KD_1/KD_1.fastq | samtools view -hSb - | samtools sort -@ 2 -m 5G -o $ANALYSIS/METTL3_KD/minimap/KD_1.bam"
bsub -q standard -oo $BASEDIR/benchmarks/epinano/logs/minimap_%J.sout.log -e $BASEDIR/benchmarks/epinano/logs/minimap_%J.serr.log -M 20000 -n 10 -R 'rusage[mem=20000]' "singularity exec $SING minimap2 -ax map-ont -t 10 -uf -k14 --for-only $ANALYSIS/METTL3_KD/reference_transcriptome.fa $ANALYSIS/METTL3_KD/albacore/KD_2/KD_2.fastq | samtools view -hSb - | samtools sort -@ 2 -m 5G -o $ANALYSIS/METTL3_KD/minimap/KD_2.bam"
bsub -q standard -oo $BASEDIR/benchmarks/epinano/logs/minimap_%J.sout.log -e $BASEDIR/benchmarks/epinano/logs/minimap_%J.serr.log -M 20000 -n 10 -R 'rusage[mem=20000]' "singularity exec $SING minimap2 -ax map-ont -t 10 -uf -k14 --for-only $ANALYSIS/METTL3_KD/reference_transcriptome.fa $ANALYSIS/METTL3_KD/albacore/WT_1/WT_1.fastq | samtools view -hSb - | samtools sort -@ 2 -m 5G -o $ANALYSIS/METTL3_KD/minimap/WT_1.bam"
bsub -q standard -oo $BASEDIR/benchmarks/epinano/logs/minimap_%J.sout.log -e $BASEDIR/benchmarks/epinano/logs/minimap_%J.serr.log -M 20000 -n 10 -R 'rusage[mem=20000]' "singularity exec $SING minimap2 -ax map-ont -t 10 -uf -k14 --for-only $ANALYSIS/METTL3_KD/reference_transcriptome.fa $ANALYSIS/METTL3_KD/albacore/WT_2/WT_2.fastq | samtools view -hSb - | samtools sort -@ 2 -m 5G -o $ANALYSIS/METTL3_KD/minimap/WT_2.bam"


for bam in $ANALYSIS/METTL3_KD/minimap/*bam; do
	name=$(basename $bam .bam)
	python bam_filter.py $bam 500 $ANALYSIS/METTL3_KD/minimap/${name}_filtered.bam
	samtools index $ANALYSIS/METTL3_KD/minimap/${name}_filtered.bam
done


# SAM2TSV
bsub -q standard -oo $BASEDIR/benchmarks/epinano/logs/sam2tsv_%J.sout.log -e $BASEDIR/benchmarks/epinano/logs/minimap_%J.serr.log -M 20000 -n 1 -R 'rusage[mem=20000]' "samtools view -h -F 3844 $ANALYSIS/METTL3_KD/minimap/WT_1_filtered.bam | sam2tsv -r $ANALYSIS/METTL3_KD/reference_transcriptome.fa > $ANALYSIS/METTL3_KD/sam2tsv/WT_1.tsv"
bsub -q standard -oo $BASEDIR/benchmarks/epinano/logs/sam2tsv_%J.sout.log -e $BASEDIR/benchmarks/epinano/logs/minimap_%J.serr.log -M 20000 -n 1 -R 'rusage[mem=20000]' "samtools view -h -F 3844 $ANALYSIS/METTL3_KD/minimap/WT_2_filtered.bam | sam2tsv -r $ANALYSIS/METTL3_KD/reference_transcriptome.fa > $ANALYSIS/METTL3_KD/sam2tsv/WT_2.tsv"
bsub -q standard -oo $BASEDIR/benchmarks/epinano/logs/sam2tsv_%J.sout.log -e $BASEDIR/benchmarks/epinano/logs/minimap_%J.serr.log -M 20000 -n 1 -R 'rusage[mem=20000]' "samtools view -h -F 3844 $ANALYSIS/METTL3_KD/minimap/KD_1_filtered.bam | sam2tsv -r $ANALYSIS/METTL3_KD/reference_transcriptome.fa > $ANALYSIS/METTL3_KD/sam2tsv/KD_1.tsv"
bsub -q standard -oo $BASEDIR/benchmarks/epinano/logs/sam2tsv_%J.sout.log -e $BASEDIR/benchmarks/epinano/logs/minimap_%J.serr.log -M 20000 -n 1 -R 'rusage[mem=20000]' "samtools view -h -F 3844 $ANALYSIS/METTL3_KD/minimap/KD_2_filtered.bam | sam2tsv -r $ANALYSIS/METTL3_KD/reference_transcriptome.fa > $ANALYSIS/METTL3_KD/sam2tsv/KD_2.tsv"



mkdir -p $BASEDIR/benchmarks/epinano/logs/
mkdir -p $ANALYSIS/METTL3_KD/variants/
bsub -q standard -oo $BASEDIR/benchmarks/epinano/logs/tsv_to_var_%J.sout.log -e $BASEDIR/benchmarks/epinano/logs/%J.serr.log -M 120000 -n 2 -R 'rusage[mem=120000]'  "source /hps/nobackup/enright/tom/miniconda3/etc/profile.d/conda.sh && conda activate $BASEDIR/conda_environments/epinano && python $BASEDIR/benchmarks/epinano/bin/EpiNano-python3/scripts/TSV_to_Variants_Freq.py3 -f $ANALYSIS/METTL3_KD/sam2tsv/WT_1.tsv -t 2 > $ANALYSIS/METTL3_KD/variants/var_freq_WT_1.tsv"
bsub -q standard -oo $BASEDIR/benchmarks/epinano/logs/tsv_to_var_%J.sout.log -e $BASEDIR/benchmarks/epinano/logs/%J.serr.log -M 120000 -n 2 -R 'rusage[mem=120000]'  "source /hps/nobackup/enright/tom/miniconda3/etc/profile.d/conda.sh && conda activate $BASEDIR/conda_environments/epinano && python $BASEDIR/benchmarks/epinano/bin/EpiNano-python3/scripts/TSV_to_Variants_Freq.py3 -f $ANALYSIS/METTL3_KD/sam2tsv/WT_2.tsv -t 2 > $ANALYSIS/METTL3_KD/variants/var_freq_WT_2.tsv"
bsub -q standard -oo $BASEDIR/benchmarks/epinano/logs/tsv_to_var_%J.sout.log -e $BASEDIR/benchmarks/epinano/logs/%J.serr.log -M 120000 -n 2 -R 'rusage[mem=120000]'  "source /hps/nobackup/enright/tom/miniconda3/etc/profile.d/conda.sh && conda activate $BASEDIR/conda_environments/epinano && python $BASEDIR/benchmarks/epinano/bin/EpiNano-python3/scripts/TSV_to_Variants_Freq.py3 -f $ANALYSIS/METTL3_KD/sam2tsv/KD_1.tsv -t 2 > $ANALYSIS/METTL3_KD/variants/var_freq_KD_1.tsv"
bsub -q standard -oo $BASEDIR/benchmarks/epinano/logs/tsv_to_var_%J.sout.log -e $BASEDIR/benchmarks/epinano/logs/%J.serr.log -M 120000 -n 2 -R 'rusage[mem=120000]'  "source /hps/nobackup/enright/tom/miniconda3/etc/profile.d/conda.sh && conda activate $BASEDIR/conda_environments/epinano && python $BASEDIR/benchmarks/epinano/bin/EpiNano-python3/scripts/TSV_to_Variants_Freq.py3 -f $ANALYSIS/METTL3_KD/sam2tsv/KD_2.tsv -t 2 > $ANALYSIS/METTL3_KD/variants/var_freq_KD_2.tsv"


bsub -q standard -oo $BASEDIR/benchmarks/epinano/logs/svm_%J.sout.log -e $BASEDIR/benchmarks/epinano/logs/svm_%J.serr.log -M 50000 -n 1 -R 'rusage[mem=120000]'  "source /hps/nobackup/enright/tom/miniconda3/etc/profile.d/conda.sh && conda activate $BASEDIR/conda_environments/epinano && python  $BASEDIR/benchmarks/epinano/bin/EpiNano-python3/scripts/SVM.py -M $BASEDIR/benchmarks/epinano/bin/EpiNano-python3/models/SVM20/model2.1-q3.mis1.mis2.mis3.mis4.mis5.del1.del2.del3.del4.del5.poly.dump -p $ANALYSIS/METTL3_KD/sam2tsv/WT_1.tsv.per.site.var.per_site_var.5mer.csv -cl 7,10,11,12,13,14,20,21,22,23,24 -o $ANALYSIS/METTL3_KD/epinano_predictions/WT_1.txt"

bsub -q standard -oo $BASEDIR/benchmarks/epinano/logs/svm_%J.sout.log -e $BASEDIR/benchmarks/epinano/logs/svm_%J.serr.log -M 50000 -n 1 -R 'rusage[mem=120000]'  "source /hps/nobackup/enright/tom/miniconda3/etc/profile.d/conda.sh && conda activate $BASEDIR/conda_environments/epinano && python  $BASEDIR/benchmarks/epinano/bin/EpiNano-python3/scripts/SVM.py -M $BASEDIR/benchmarks/epinano/bin/EpiNano-python3/models/SVM20/model2.1-q3.mis1.mis2.mis3.mis4.mis5.del1.del2.del3.del4.del5.poly.dump -p $ANALYSIS/METTL3_KD/sam2tsv/WT_2.tsv.per.site.var.per_site_var.5mer.csv -cl 7,10,11,12,13,14,20,21,22,23,24 -o $ANALYSIS/METTL3_KD/epinano_predictions/WT_2.txt"

bsub -q standard -oo $BASEDIR/benchmarks/epinano/logs/svm_%J.sout.log -e $BASEDIR/benchmarks/epinano/logs/svm_%J.serr.log -M 50000 -n 1 -R 'rusage[mem=120000]'  "source /hps/nobackup/enright/tom/miniconda3/etc/profile.d/conda.sh && conda activate $BASEDIR/conda_environments/epinano && python  $BASEDIR/benchmarks/epinano/bin/EpiNano-python3/scripts/SVM.py -M $BASEDIR/benchmarks/epinano/bin/EpiNano-python3/models/SVM20/model2.1-q3.mis1.mis2.mis3.mis4.mis5.del1.del2.del3.del4.del5.poly.dump -p $ANALYSIS/METTL3_KD/sam2tsv/KD_1.tsv.per.site.var.per_site_var.5mer.csv -cl 7,10,11,12,13,14,20,21,22,23,24 -o $ANALYSIS/METTL3_KD/epinano_predictions/KD_1.txt"

bsub -q standard -oo $BASEDIR/benchmarks/epinano/logs/svm_%J.sout.log -e $BASEDIR/benchmarks/epinano/logs/svm_%J.serr.log -M 50000 -n 1 -R 'rusage[mem=120000]'  "source /hps/nobackup/enright/tom/miniconda3/etc/profile.d/conda.sh && conda activate $BASEDIR/conda_environments/epinano && python  $BASEDIR/benchmarks/epinano/bin/EpiNano-python3/scripts/SVM.py -M $BASEDIR/benchmarks/epinano/bin/EpiNano-python3/models/SVM20/model2.1-q3.mis1.mis2.mis3.mis4.mis5.del1.del2.del3.del4.del5.poly.dump -p $ANALYSIS/METTL3_KD/sam2tsv/KD_2.tsv.per.site.var.per_site_var.5mer.csv -cl 7,10,11,12,13,14,20,21,22,23,24 -o $ANALYSIS/METTL3_KD/epinano_predictions/KD_2.txt"


join -a 1 -1 1 -2 1 -t $'\t' \
	<(tail -n +2 $ANALYSIS/METTL3_KD/epinano_predictions/WT_2*.csv | awk 'BEGIN{FS=","; OFS="\t"}{split($2, A, ":"); split($4, C, ":"); print $3"_"A[1]-1,$3,$1,A[1]-1,C[1],$27}' | sort -k1,1) \
	<(awk 'BEGIN{OFS=FS="\t"}NR>1{print $4"_"$1,$1,$3,$4,$6}' $BASEDIR/nanocompore_pipelines/METTL3_KD_polyA/results/nanocompore/out_nanocompore_results.tsv | sort -k1,1 | grep -v '!!!') | \
	cut -f2,3,4,5,6,7,8,10 > $ANALYSIS/METTL3_KD/epinano_predictions/WT_2_parsed.txt

join -a 1 -1 1 -2 1 -t $'\t' \
	<(tail -n +2 $ANALYSIS/METTL3_KD/epinano_predictions/WT_1*.csv | awk 'BEGIN{FS=","; OFS="\t"}{split($2, A, ":"); split($4, C, ":"); print $3"_"A[1]-1,$3,$1,A[1]-1,C[1],$27}' | sort -k1,1) \
	<(awk 'BEGIN{OFS=FS="\t"}NR>1{print $4"_"$1,$1,$3,$4,$6}' $BASEDIR/nanocompore_pipelines/METTL3_KD_polyA/results/nanocompore/out_nanocompore_results.tsv | sort -k1,1 | grep -v '!!!') | \
	cut -f2,3,4,5,6,7,8,10 > $ANALYSIS/METTL3_KD/epinano_predictions/WT_1_parsed.txt

join -a 1 -1 1 -2 1 -t $'\t' \
	<(tail -n +2 $ANALYSIS/METTL3_KD/epinano_predictions/KD_2*.csv | awk 'BEGIN{FS=","; OFS="\t"}{split($2, A, ":"); split($4, C, ":"); print $3"_"A[1]-1,$3,$1,A[1]-1,C[1],$27}' | sort -k1,1) \
	<(awk 'BEGIN{OFS=FS="\t"}NR>1{print $4"_"$1,$1,$3,$4,$6}' $BASEDIR/nanocompore_pipelines/METTL3_KD_polyA/results/nanocompore/out_nanocompore_results.tsv | sort -k1,1 | grep -v '!!!') | \
	cut -f2,3,4,5,6,7,8,10 > $ANALYSIS/METTL3_KD/epinano_predictions/KD_2_parsed.txt

join -a 1 -1 1 -2 1 -t $'\t' \
	<(tail -n +2 $ANALYSIS/METTL3_KD/epinano_predictions/KD_1*.csv | awk 'BEGIN{FS=","; OFS="\t"}{split($2, A, ":"); split($4, C, ":"); print $3"_"A[1]-1,$3,$1,A[1]-1,C[1],$27}' | sort -k1,1) \
	<(awk 'BEGIN{OFS=FS="\t"}NR>1{print $4"_"$1,$1,$3,$4,$6}' $BASEDIR/nanocompore_pipelines/METTL3_KD_polyA/results/nanocompore/out_nanocompore_results.tsv | sort -k1,1 | grep -v '!!!') | \
	cut -f2,3,4,5,6,7,8,10 > $ANALYSIS/METTL3_KD/epinano_predictions/KD_1_parsed.txt
