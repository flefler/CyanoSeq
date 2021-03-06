# CyanoSeq

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.6745183.svg)](https://doi.org/10.5281/zenodo.6745183)

Version 1.0.0 available on Zenodo at https://zenodo.org/record/6745183

Current version: 0.9.6

CyanoSeq is a curated database of cyanobacterial 16S rRNA sequences for taxonomic assignment of metagenomic/metabarcoding/amplicon reads. CyanoSeq is assembled from 16S rRNA sequences found within NCBI, with their taxonomies curated from cyanobacterial taxonomic literature as well as a systematic assessment of uncharacterized cyanobacterial sequences. When possible, the full length 16S rRNA sequences are provided, allowing use for several 16S rRNA primer sets to be used for metabarcoding, as well as full length 16S rRNA sequences for taxonomic assignment. The taxonomy of CyanoSeq is meant to reflect the current state of cyanobacterial taxonomy with curated clades of described and undescribed taxa. A provisional rank was given to those taxa that fell outside of the sensu stricto clade in an attempt to resolve polyphletic ranks. CyanoSeq does not aim to revise cyanobacterial taxonomy nor become a taxonomic authority, rather it serves as a starting point to identify and name monophyletic clades which do not belong to any established taxonomic rank and require revision. CyanoSeq currently contains 5365 cyanobacterial sequences and 123 chloroplast and bacterial sequences for use in classifying reads from metagenomic studies.


Two fastq.qz files are provdied for taxonomic assignment using the "assignTaxonomy" function in [DADA2.](https://benjjneb.github.io/dada2/tutorial.html). These files should work with QIIME2, but have not currently been tested.

[CyanoSeq0.9.6_dada2.fastq.gz](https://github.com/flefler/CyanoSeq/blob/main/CyanoSeq0.9.6_dada2.fastq.gz) is the Cyanobacterial data bases which contains 5365 cyanobacterial sequences and 123 chloroplast and bacterial sequences and can be used with cyanobacterial specific primers (i.e., those described by [Nübel et al., 1997](https://journals.asm.org/doi/10.1128/aem.63.8.3327-3332.1997)) 

[CyanoSeq0.9.6_SILVA138.1_dada2.fastq.gz](https://github.com/flefler/CyanoSeq/blob/6c8f3b669936083a366267d691e222cc6569517f/CyanoSeq0.9.6_SILVA138.1_dada2.fastq.gz) is the cyanobacterial database merged with [SILVA](https://www.arb-silva.de/), current version 138.1, with the cyanobacterial sequences from SILVA removed and replaced with those curated here, and the Class “Cyanobacteriia” replaced with “Cyanophyceae”. This can be used general bacterial primers to fully understand the bacterial and cyanobacterial communities. 

# Questions, comments, concerns?

Leave a request or start a discussion here, or email me at flefler(at)ufl(dot)edu
