# CyanoSeq V1.3

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.7864137.svg)](https://doi.org/10.5281/zenodo.7864137)

[Current version: 1.3](https://zenodo.org/record/)

CyanoSeq is published in the Journal of Phycology: https://doi.org/10.1111/jpy.13335

CyanoSeq is a curated database of cyanobacterial 16S rRNA sequences for taxonomic assignment of metagenomic/metabarcoding/amplicon reads. CyanoSeq is assembled from 16S rRNA sequences found within NCBI, with their taxonomies curated from cyanobacterial taxonomic literature as well as a systematic assessment of uncharacterized cyanobacterial sequences. When possible, the full length 16S rRNA sequences are provided, allowing for several 16S rRNA primer sets to be used for taxonomic assignment of metabarcoding data, as well as de novo phylogenetic reconstruction. The taxonomy of CyanoSeq is meant to reflect the current state of cyanobacterial taxonomy with curated clades of described and undescribed taxa. A provisional rank was given to those taxa that fell outside of the sensu stricto clade in an attempt to resolve polyphyletic ranks. CyanoSeq does not aim to revise cyanobacterial taxonomy nor become a taxonomic authority, rather it serves as a starting point to identify and name monophyletic clades which do not belong to any established taxonomic rank and require revision. CyanoSeq currently contains 4174 cyanobacterial sequences and 123 chloroplast and bacterial sequences for use in classifying reads from metabarcoding studies.

This update was done in conjuction with an update to the cyanobacterial taxonomy of [ITIS](https://itis.gov/) using resources such as [AlgaeBase](AlgaeBase.org), [CyanoDB](www.cyanodb.cz/), as well as recent literature. 

## Key updates V1.3
A few changes have been made since the last version, which are noted in the [change log](https://github.com/flefler/CyanoSeq/blob/main/ChangeLog.md). A few key points listed below.

1: We have switched from using SILVA 138.1 as the bacterial database, we have now incorperated the bacteria, excluding cyanobacteria, from [GSR-DB](https://manichanh.vhir.org/gsrdb/index.php). Please cite ```Leidy-Alejandra G. Molano, Sara Vega-Abellaneda, Chaysavanh Manichanh. GSR-DB: a manually curated and optimized taxonomical database for 16S rRNA amplicon analysis. mSystems (2024) https://doi.org/10.1128/msystems.00950-23``` in addition to CyanoSeq if you use this version

2: The family Prochlorococcaceae has proven to be a headache for curation and now has its taxonomy modeled after GTDB R220 with some minor changes, see the Taxonomy_V1.3.xlsx files for more in depth information.

3: Greatly reduced the number of sequences, primarily in over represented genera (e.g., <I>Dolichospermum</I>, <I>Prochlorococcus</I>, etc) to reduce redundancy. This <I>should</I> help with classification, especially with the difficult groups (e.g., ADA-<I>Aphanizomenon</I>/<I>Dolichospermum</I>/<I>Anabaena</I>)

4: To facilite de novo phylogenetic reconstruction, the file ```NCBI_ClassifiedSeqs.tsv``` is provided on Zenodo. [Example ```R``` script for manipulation can be found here](https://github.com/flefler/CyanoSeq/blob/main/RetrivingSeqs.md). Script on how this was done can be found [here](https://github.com/flefler/CyanoSeq/blob/main/ScrapeNCBI_Classify.md).

## Files

Two fastq.qz files are provdied for taxonomic assignment using the "assignTaxonomy" function in [DADA2](https://benjjneb.github.io/dada2/tutorial.html) and [IdTaxa](https://microbiomejournal.biomedcentral.com/articles/10.1186/s40168-018-0521-5) classifiers. IDTAXA files have not been tested, please let me know if these work or not.

Scripts are now provided to create [QIIME2](https://docs.qiime2.org/2022.8/) classifiers. Thanks to Lucija Kranjer. Necessary files for are provided.

CyanoSeq_1.3_dada2.fastq.gz is the Cyanobacterial data bases which contains 4174 cyanobacterial sequences with 123 chloroplast and bacterial sequences. This should only be used with cyanobacterial specific primers (i.e., those described by [Nübel et al., 1997](https://journals.asm.org/doi/10.1128/aem.63.8.3327-3332.1997)) 

CyanoSeq_1.3_GSRDB_dada2.fastq.gz is the cyanobacterial database merged with [GSR-DB](https://manichanh.vhir.org/gsrdb/index.php), with the cyanobacterial sequences from GSR-DB removed and replaced with those curated here. This can be used general bacterial primers to characterize the total bacterial community. 

Fasta and nwk files of each order are provided to facilitate de novo phylogenetic tree reconstruction for novel sequences and use of tools such as [epa-ng](https://github.com/pierrebarbera/epa-ng) for placement of your ASVs

# Questions, comments, concerns?

Leave a request or start a discussion here
