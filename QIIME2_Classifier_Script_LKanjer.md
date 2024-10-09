# The below script was provided by Lucija Kanjer, ORCID 0000-0003-3869-8742. Tell her thanks!
```
##### CyanoSeq QIIME2 taxonomy files and classifiers #####
# 0n 26 June 2023 by L. Kanjer

# activate QIIME2 environment
conda activate qiime2-2023.2

# Import sequence data into QZA format
qiime tools import \
  --type 'FeatureData[Sequence]' \
  --input-path CyanoSeq_1.2.1_SILVA138.1_QIIME2.fasta \
  --output-path CyanoSeq_1.2.1_SILVA138.1_QIIME2_sequences.qza

# Import taxonomy data into QZA format
qiime tools import \
  --type 'FeatureData[Taxonomy]' \
  --input-format HeaderlessTSVTaxonomyFormat \
  --input-path CyanoSeq_1.2.1_SILVA138.1_QIIME2.tsv \
  --output-path CyanoSeq_1.2.1_SILVA138.1_QIIME2_taxonomy.qza

# Extract selected V4 reads (515F/806R) from the CyanoSeq database as reference reads
# parameters are the same as in QIIME2 tutorial
qiime feature-classifier extract-reads \
  --i-sequences CyanoSeq_1.2.1_SILVA138.1_QIIME2_sequences.qza \
  --p-f-primer GTGCCAGCMGCCGCGGTAA \
  --p-r-primer GGACTACHVGGGTWTCTAAT \
  --p-trunc-len 120 \
  --p-min-length 100 \
  --p-max-length 400 \
  --o-reads ref-seqs-CyanoSeq-515F-806R.qza

# Train classifier - adding taxonomy assignment to extracted V4 reads (515F/806R)
qiime feature-classifier fit-classifier-naive-bayes \
  --i-reference-reads ref-seqs-CyanoSeq-515F-806R.qza \
  --i-reference-taxonomy CyanoSeq_1.2.1_SILVA138.1_QIIME2_taxonomy.qza \
  --o-classifier classifier-CyanoSeq-515F-806R.qza

# Extract selected V3-V4 reads (341F/805R) from the CyanoSeq database as reference reads
qiime feature-classifier extract-reads \
  --i-sequences CyanoSeq_1.2.1_SILVA138.1_QIIME2_sequences.qza \
  --p-f-primer CCTACGGGNGGCWGCAG \
  --p-r-primer GACTACHVGGGTATCTAATCC \
  --p-min-length 250 \
  --p-max-length 500 \
  --o-reads ref-seqs-CyanoSeq-341F-805R.qza

# Train classifier - adding taxonomy assignment to extracted V3-V4 reads (341F/805R)
qiime feature-classifier fit-classifier-naive-bayes \
  --i-reference-reads ref-seqs-CyanoSeq-341F-805R.qza \
  --i-reference-taxonomy CyanoSeq_1.2.1_SILVA138.1_QIIME2_taxonomy.qza \
  --o-classifier classifier-CyanoSeq-341F-805R.qza
```
# This is how I made the basic classifiers

## CyanoSeq SILVA
```
qiime tools import \
  --type 'FeatureData[Sequence]' \
  --input-path CyanoSeqV1.3_SILVA138.2_QIIME2.fasta \
  --output-path CyanoSeqV1.3_SILVA138.2_QIIME2.qza

qiime tools import \
  --type 'FeatureData[Taxonomy]' \
  --input-format HeaderlessTSVTaxonomyFormat \
  --input-path CyanoSeqV1.3_SILVA138.2_QIIME2_tax.tsv \
  --output-path ref-taxonomy.qza

qiime feature-classifier fit-classifier-naive-bayes \
  --i-reference-reads CyanoSeqV1.3_SILVA138.2_QIIME2.qza \
  --i-reference-taxonomy ref-taxonomy.qza \
  --o-classifier CyanoSeqV1.3_SILVA138.2_QIIME2_classifier.qza
```
## CyanoSeq GSR-DB
```
qiime tools import \
  --type 'FeatureData[Sequence]' \
  --input-path CyanoSeqV1.3_GSRDB_QIIME2.fasta \
  --output-path CyanoSeqV1.3_GSRDB_QIIME2.qza

qiime tools import \
  --type 'FeatureData[Taxonomy]' \
  --input-format HeaderlessTSVTaxonomyFormat \
  --input-path CyanoSeqV1.3_GSRDB_QIIME2_tax.tsv \
  --output-path ref-taxonomy.qza

qiime feature-classifier fit-classifier-naive-bayes \
  --i-reference-reads CyanoSeqV1.3_GSRDB_QIIME2.qza \
  --i-reference-taxonomy ref-taxonomy.qza \
  --o-classifier CyanoSeqV1.3_GSRDB_QIIME2_classifier.qza
```
