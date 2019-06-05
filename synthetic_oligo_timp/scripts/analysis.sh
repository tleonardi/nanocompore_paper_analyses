#/bin/bash
BASEDIR=$(git rev-parse --show-toplevel)
CONDA_ENVS="$BASEDIR/conda_environments/"
ANALYSIS="${BASEDIR}/synthetic_oligo_timp/analysis"
IMAGE="$BASEDIR/synthetic_oligo_timp/nanocompore_af97cf1.img"
OLDIMAGE="/hps/nobackup/enright/tom/new_nanocompore/7SK_DKC1/nanocompore_pipeline/singularity/nanocompore_b93c6c2.img"

# Download the sequencing data
mkdir -p $ANALYSIS/data
wget -O $ANALYSIS/data/3primem6aoligo.tar.gz https://s3.amazonaws.com/timpnanopore/oxford/171208_m6A/3primem6aoligo.tar.gz

mkdir -p $ANALYSIS/data/3primem6aoligo
tar -zxf $ANALYSIS/data/3primem6aoligo.tar.gz -C $ANALYSIS/data/3primem6aoligo


# BASECALL WITH GUPPY
export cpus=30
export flowcell=FLO-MIN106
export kit=SQK-RNA001
mkdir -p $ANALYSIS/guppy

GUPPY_COMMAND="guppy_basecaller -i $ANALYSIS/data/3primem6aoligo/ -s $ANALYSIS/guppy --fast5_out --recursive --num_callers ${cpus} --disable_pings --reverse_sequence true --u_substitution true --trim_strategy rna --flowcell $flowcell --kit $kit"

singularity exec $IMAGE $GUPPY_COMMAND


# MAP TO THE REFERENCE
mkdir -p $ANALYSIS/minimap/
MINIMAP_COMMAND="minimap2 -x map-ont -t ${cpus} -a $ANALYSIS/data/ligated_oligo.fa $ANALYSIS/guppy/*.fastq"
singularity exec $IMAGE $MINIMAP_COMMAND > $ANALYSIS/minimap/minimap.sam

awk 'BEGIN{OFS=FS="\t"}{print $1, $2-70, $2, ".", 0, "+"}' $ANALYSIS/data/ligated_oligo.fa.fai > $ANALYSIS/data/oligo_ROI.bed
samtools view $ANALYSIS/minimap/minimap.sam -bh -t $ANALYSIS/data/ligated_oligo.fa.fai -F 2324 -L $ANALYSIS/data/oligo_ROI.bed | samtools sort -@ ${cpus} -o $ANALYSIS/minimap/minimap_filtered.bam
samtools index $ANALYSIS/minimap/minimap_filtered.bam $ANALYSIS/minimap/minimap_filtered.bam.bai

mkdir -p $ANALYSIS/nanopolish
cat $ANALYSIS/guppy/*.fastq > $ANALYSIS/nanopolish/basecalled.fastq
NANOPOLISH_COMMAND="nanopolish index -s $ANALYSIS/guppy/sequencing_summary.txt -d $ANALYSIS/data/3primem6aoligo/ $ANALYSIS/nanopolish/basecalled.fastq"
singularity exec $IMAGE $NANOPOLISH_COMMAND

mkdir -p $ANALYSIS/nanopolishcomp/
EVENTALIGN_COMMAND="nanopolish eventalign -t 15 --reads $ANALYSIS/nanopolish/basecalled.fastq --bam $ANALYSIS/minimap/minimap_filtered.bam --genome $ANALYSIS/data/ligated_oligo.fa --samples --print-read-names --scale-events > $ANALYSIS/nanopolish/eventalign.txt"
singularity exec $IMAGE bash -c "$EVENTALIGN_COMMAND"

# reformat nanopolish results
singularity exec $IMAGE python3 $BASEDIR/synthetic_oligo_timp/scripts/reformat_eventalign.py $ANALYSIS/nanopolish/eventalign.txt $ANALYSIS/nanopolishcomp/

NANOPOLISHCOMP1_COMMAND="NanopolishComp Eventalign_collapse -i $ANALYSIS/nanopolish/mod_eventalign.txt -t 15 -o $ANALYSIS/nanopolishcomp/ -p mod"
NANOPOLISHCOMP2_COMMAND="NanopolishComp Eventalign_collapse -i $ANALYSIS/nanopolish/unmod_eventalign.txt -t 15 -o $ANALYSIS/nanopolishcomp/ -p unmod"
singularity exec $IMAGE bash -c "$NANOPOLISHCOMP1_COMMAND"
singularity exec $IMAGE bash -c "$NANOPOLISHCOMP2_COMMAND"

mkdir $ANALYSIS/nanocompore/
echo -e ">FLuc_Control_Plasmid:88-1805\nTGAGGACTGTA" > $ANALYSIS/nanocompore/ref.fa

NANOCOMPORE_COMMAND="nanocompore sampcomp --file_list1 $ANALYSIS/nanopolishcomp/mod_eventalign_collapse.tsv --file_list2 $ANALYSIS/nanopolishcomp/unmod_eventalign_collapse.tsv --label1 mod --label2 unmod --fasta $ANALYSIS/nanocompore/ref.fa --outpath $ANALYSIS/nanocompore/ --min_ref_length 1 --sequence_context 1 --pvalue_thr 0 --overwrite"
singularity exec $OLDIMAGE bash -c "$NANOCOMPORE_COMMAND"

