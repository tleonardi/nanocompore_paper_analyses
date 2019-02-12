# Nanocompore paper analyses

Analyses performed for the nanocompore paper

## Draft paper

https://www.authorea.com/users/34664/articles/334998-comparative-detection-of-rna-modifications-with-nanopore-direct-rna-sequencing#

## In  silico Benchmarking

### In silico reference generation

To assess the accuracy of Nanocompore, we generated a set of in silico reference transcripts. In order to cover all the possible sequence context, we first generated a De Bruijn sequence with all possible 7 based long kmers.
The Bruijn sequence is the shortest circular sequence of length such that every string of length on the alphabet of size occurs as a contiguous subrange of the sequence. However it lacks the complexity and randomness of real biologica data. In addition we also wanted to create a high level of redundancy to have the same kmers in various larger sequence contexts. Thus, the De Bruijn sequence was sliced into fragments of 20 bases with a 1 base step sliding window. A set of 2000 sequences 500 bases long were then created by joining randomly selected fragments together. The selected set covers all of the 7-mers at least XXXXXXXXXX times

The full analysis is available in the following python notebook: [De_Bruijn_shuffle_seq_gen]() as well as the [fasta_reference file]()

### Generate artificial modification data

Nanocompore comes with a companion tool called `simulate_reads_from_fasta`. This tool can generate a NanopolishComp Eventalign like file corresponding to a set of simulated reads matching each sequences in a reference fasta file. The intensities and dwell time per 5-mers are based on a model file provided in Nanocompore.

The model was generated based on a non modified RNA standard available ... 

This analysis is available at in the following python notebook: [De_Bruijn_shuffle_seq_gen]() 

...


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
