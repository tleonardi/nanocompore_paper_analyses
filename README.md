# Nanocompore paper analyses

This repository contains all the analyses performed for the [Nanocompore](https://github.com/tleonardi/nanocompore) paper

## Draft paper

https://www.authorea.com/users/34664/articles/334998-comparative-detection-of-rna-modifications-with-nanopore-direct-rna-sequencing

## In silico benchmarking

### In silico reference

To assess the accuracy of Nanocompore, we generated a set of in silico reference sequences. In order to maximise the sequence diversity and kmer coverage we wrote a "guided" random sequence generator. In brief, the sequences are generated base per base using a randon function, but the program keeps a track of the number of times each kmers was already used. The sequence is extended, based on a random function with a weighted probability for each kmer invertly proportional to their occurence in the sequences already generated. Essentially, this ensure that all kmers are represented as uniformly as possible, but it leaves some limited space to randomness. We generated a set of 2000 sequences 500 bases long each that maximises the 9-mer coverage. We also excluded any homopolymer longer than 5 bases, as they are likely to be miscalled in nanopore data. All sequence were aligned in a pairwise fashion to ensure they are sufficiently different from each other. Kmer coverage in the final sequence set are summarised in the table below.  

| Kmer length | % kmer found | Median  occurences |
| ----------- | ------------ | ------------------ |
| 5           | 100.00%      | 970                |
| 7           | 99.83%       | 60                 |
| 9           | 99.66%       | 4                  |

* Analysis python notebook: [Random_guided_seq_gen.ipynb](https://github.com/a-slide/nanocompore_paper_analyses/blob/master/in_silico_dataset/Random_guided_seq_gen.ipynb)
* Fasta reference : [random weighted function.fa](https://raw.githubusercontent.com/a-slide/nanocompore_paper_analyses/master/in_silico_dataset/random weighted function.fa)

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
