# This script is to scrape cyanobacterial 16S rRNA sequences off NCBI

## Load Libraries
```
library(rentrez)
library(Biostrings)
```
## Set query and other fun stuff
### Last done 2024/09/17:3000[PDAT])
```
query <- "(((txid1117[ORGN] AND 16S[Title]) NOT Genome[All Fields]) NOT Genomic[All Fields] AND ddbj_embl_genbank[filter])"
search_results <- entrez_search(db = "nucleotide", term = query, retmax = 90000)
id_list <- search_results$ids
batch_size <- 100
id_batches <- split(id_list, ceiling(seq_along(id_list)/batch_size))

## Download sequences in batches
sequences <- character()
for (batch in id_batches) {
  sequences_batch <- entrez_fetch(db = "nucleotide", id = batch, rettype = "fasta")
  sequences <- paste(sequences, sequences_batch)
  Sys.sleep(2) # add a small delay between batches to avoid overloading the server
}
```
## Write the sequences to a fasta file, which needs the ```"``` characters removed, how annoying.
```
write.table(sequences, file = "output.txt")
```
# Classify all sequences from NCBI

## Load in sequences, keep only sequences 1000\<n\<2000

```{r message=FALSE}
library(Biostrings)
NCBI_Seqs = readDNAStringSet("~/UFL Dropbox/Forrest Lefler/Laughinghouse_Lab/CyanoSeq/External_databases/NCBI_20240917_cyano_16S.fa")

long_seqs <- NCBI_Seqs[width(NCBI_Seqs) >= 1000]
long_seqs2 <- long_seqs[width(long_seqs) <= 2000]

head(long_seqs2)

Biostrings::writeXStringSet(long_seqs2, file = "/Users/flefler/UFL Dropbox/Forrest Lefler/Laughinghouse_Lab/CyanoSeq/External_databases/NCBI_20240917_cyano_16S_FILTERED.fa")
```

## Classify and filter the seqs, save the df

```{r message=FALSE}
#This gives dada2 the path to the new sequecnes
NCBI_sequences = dada2::getSequences("/Users/flefler/UFL Dropbox/Forrest Lefler/Laughinghouse_Lab/CyanoSeq/External_databases/NCBI_20240917_cyano_16S_FILTERED.fa", collapse = FALSE, silence = TRUE)

#This assigns taxonomy with a minimum 97% bootstrap confidence
taxa <- dada2::assignTaxonomy(NCBI_sequences, "~/UFL Dropbox/Forrest Lefler/Laughinghouse_Lab/CyanoSeq/V1.3/CyanoSeq_1.3_dada2.fastq.gz", minBoot=80, multithread=2, outputBootstraps = T)

#This creates a dataframe
Seq_DF = as.data.frame(taxa[["tax"]]) %>% rownames_to_column(var = "Sequence")

fasta = seqinr::read.fasta(file = "/Users/flefler/UFL Dropbox/Forrest Lefler/Laughinghouse_Lab/CyanoSeq/External_databases/NCBI_20240917_cyano_16S_FILTERED.fa",
                           as.string = TRUE, forceDNAtolower = F, set.attributes = T, whole.header = T)

#Create data fram from fasta file with only accession number and sequence
Genbank_accession = names(fasta)
Sequence = paste(fasta)
df = data.frame(Genbank_accession,Sequence)


#Join the new sequence files and the taxonomy, select only the good stuff
NCBI_ClassifiedSeqs = full_join(Seq_DF, df, by = "Sequence") %>% 
  filter(Class == "Cyanophyceae") %>% 
  filter(Order != "Chloroplast") %>% 
  filter(!is.na(Genbank_accession))

write.table(NCBI_ClassifiedSeqs, file = "~/UFL Dropbox/Forrest Lefler/Laughinghouse_Lab/CyanoSeq/V1.3/NCBI_ClassifiedSeqs.tsv", sep = "\t", row.names = F, col.names = T, quote = FALSE)
```

