# Example how to pull out sequences from the NCBI_ClassifiedSeqs.tsv

You will need some familiarity with ```R```, and the packages [tidyverse](https://www.tidyverse.org/) and [Biostrings](https://github.com/Bioconductor/Biostrings) are needed.

## Load packages and the file
```
library(tidyverse)
library(Biostrings)

NCBI_ClassifiedSeqs = read_delim("NCBI_ClassifiedSeqs.tsv", 
    delim = "\t", escape_double = FALSE, trim_ws = TRUE)
```
## This will give you a df with ONLY Microcystis
```
Microcystis = NCBI_ClassifiedSeqs %>% 
  filter(Genus == "Microcystis")
```
## This will give you a df with the family Microcystaceae
```
Microcystaceae = NCBI_ClassifiedSeqs %>% 
  filter(Family == "Microcystaceae")
```
## This will give you a df with the order 
```
Chroococcales = NCBI_ClassifiedSeqs %>% 
  filter(Order == "Chroococcales")
```
## Save the file
This saves as a fasta file with the Genbank accession number and name 
```
Microcystis = Biostrings::DNAStringSet(Microcystis$Sequence)
names(Microcystis) = paste(Microcystis$Genbank_accession)
Biostrings::writeXStringSet(Microcystis, "Microcystis.fasta")
```
