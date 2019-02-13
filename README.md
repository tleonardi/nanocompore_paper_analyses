# Nanocompore paper analyses

This repository contains all the analyses performed for the [Nanocompore](https://github.com/tleonardi/nanocompore) paper

## Draft paper

https://www.authorea.com/users/34664/articles/334998-comparative-detection-of-rna-modifications-with-nanopore-direct-rna-sequencing

## In silico benchmarking

### In silico reference

To assess the accuracy of Nanocompore, we generated a set of in silico reference sequences. In order to cover all the possible sequence context, we first generated a De Bruijn sequence with all possible 7 bases long kmers.
The Bruijn sequence is the shortest circular sequence of length such that every string of length on the alphabet of size occurs as a contiguous substring of the sequence. However it lacks the complexity and randomness of real biological data. In addition we also wanted to create a high level of redundancy to have the same kmers in various larger sequence contexts. Thus, the De Bruijn sequence was sliced into fragments of 15 bases with a 1 base step sliding window. A set of 2000 sequences 500 bases long were then created by joining randomly selected fragments together. The selected set covers all of the 7-mers at least 35 times and up to 93, vith a median coverage of 60 occurences. 

* Analysis python notebook: [De_Bruijn_shuffle_seq_gen.ipynb](https://github.com/a-slide/nanocompore_paper_analyses/blob/master/in_silico_dataset/De_Bruijn_shuffle_seq_gen.ipynb)
* Fasta reference : [shuffle_de_bruijn_ref.fa](https://raw.githubusercontent.com/a-slide/nanocompore_paper_analyses/master/in_silico_dataset/shuffle_de_bruijn_ref.fa?token=AFb-SIR5A2x5Q0Ak4qhYcfrajIHz5e5tks5cbReIwA%3D%3D)

### kmer intensity and dwell time model

In order to generate in silico reads, we first needed to obtain a model distribution for each kmers that could subsequently be altered to mimick the effect of  modifications on the signal. To do so, we used an in vitro synthesized non-modified RNA dataset obtained from the [Human Nanopore Consortium data repository](https://github.com/nanopore-wgs-consortium/NA12878/blob/master/nanopore-human-transcriptome/fastq_fast5_bulk.md). After mapping and signal realignment with Nanopolish, we collected median intensity and dwell time for each 5-mers. We then fitted a large number of distributions and selected the one minimizing the overall sum of squared errors between the real data and the model distribution. The median intensity is best modeled by a [logistic distribution](https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.logistic.html) and the dwell time by a [Wald distribution](https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.wald.html).  

- Analysis python notebook: [Kmer_Model_Generation.ipynb](https://github.com/a-slide/nanocompore_paper_analyses/blob/master/in_silico_dataset/Kmer_Model_Generation.ipynb)
- Kmers model: [kmer_model_RNA.tsv](https://raw.githubusercontent.com/a-slide/nanocompore_paper_analyses/master/in_silico_dataset/kmer_model_RNA.tsv?token=AFb-SClDNjGaXPwUmSsM9X9EgGWf-aatks5cbTpPwA%3D%3D)

### Simulated datasets generation

Nanocompore comes with a companion tool called `simulate_reads_from_fasta` which can generate simulated read data based on a fasta reference and a kmer model file (previously described). In addition one can also use this tool to generate signal alterations in reads by tweaking the kmer model to mimmick the effect of RNA modifications. We generated a large number of simulated datasets with various amplitude of modification of the median signal intensity and dwell time.

* Analysis python notebook: []()
* Dataset repository...


### Nanocompore


### Post processing analyses (threshold of detection + accuracy)


## Targeted sequencing samples

Description: ...

* Basecalling, alignment and resquigling
* Nanocompore
* Post processing analyses

## METTL3 polya+ RNAs

Description: ...

* Basecalling, alignment and resquigling
* Nanocompore
* Post processing analyses
