# Changes in 1.3

We have switched from using SILVA 138.1 as the bacterial database, we will now use [GSR-DB](https://manichanh.vhir.org/gsrdb/index.php)

The family Prochlorococcaceae has proven to be a headache for curation and now has its taxonomy modeled after [GTDB R220](https://gtdb.ecogenomic.org/) with some minor changes, see below.
There have been several provisional genera (i.e., GTDB names) added. See the Taxonomy_V1.3.xlsx files for more in depth information.

Greatly reduced the number of sequences, primarily in over represented genera (e.g., <i>Dolichospermum, Prochlorococcus</i>, etc) to reduce redundancy. 
This "Should" help with classification, especially with the difficult groups (e.g., ADA-<i>Aphanizomenon/Dolichospermum/Anabaena</i>)

GTDB R220 changes for CyanoSeq V1.3:

WH-5701 and MW101C3 = <i>Regnicoccus</i>

NIES-981 = <i>Ciimarium</i>

RCC307 = <i>Inmanicoccus</i>

Synechococcus_D = <i>Lacustricoccus</i>

Aphanocapsa_A feldmannii, <i>Regnicoccus frigidus</i>, Cyanobium_A sp013205255 were excluded as their placement in the 16S rRNA phylogeny is weird, likely due to these being pulled from MAGs.


## New genera added:
<i>Almyronema, Ciimarium, Cyanoarbor, Limnonema, Egbenema, Komarkovaeasiopsis, Gansulinema, Karukerafilum, Radiculonema,
Copelandiella, Azorthrix, Pseudocalidococcus, Phayaothrix, Trichothermofontia, Monilinema, Floridaenema,
Venetifunis, Ahomia, Parathermosynechococcus, Candidatus Sivonenia, Sphaerothrix, Reofilinostoc, Okeanomitos, Thermospirulina</i>


## Nomenclature notes:
All taxa from the family Aerosakkonemataceae were moved to the novel order "Aerosakkonematales" based on phylogenomic analyses from our accepted paper in J Phycol, DOI pending.
This also resulted in the description of <i>Floridaenema</i>, thus provisional genus "Oscillatoriaceae_X" and what we previously considered "<i>Phormidium</i>" have been with reclassified as <i>Floridaenema</i>.

<i>Funiculus tenuis</i> is now <i>Konicacronema tenuis</i> and <i>Parifilum solicrustae</i> is now <i>Pycnacronema solicrustae</i> (Jusko and Johansen 2023 DOI: 10.1111/jpy.13411)

<i>Leptothermofonsia</i> species were transferred into <i>Kovacikia </i>(Kaštovský et al. https://doi.org/10.3390/d15090975)

Rivulariaceae_XX is now <i>Phayaothrix</i>

Gloeobacteraceae_X is now Candidatus_Sivonenia

Prochlorococcaceae_XX is now NIES-981 based on GTDBR220 addition

Leptolyngbyaceae_XXX is now split between <i>Copelandiella</i> and <i>Khargia</i>

<i>Sphaerothrix</i> and <i>Almyronema</i> are very close on the tree, ignoring for now but deserves more attention

The Class "Candidatus-Melainabacteria" was renamed to "Vampirovibrionophyceae" as per Strunecký and Mareš 2023

## Problematic taxa:
<i>Brunnivagina</i> was not included as the analyses was done on a single partial 16S rRNA fragment. Genomic analyses was not conducted despite availability. Genomic evidence does not support the erection of this genus.

Kaštovský et al. (https://doi.org/10.3390/d15090975) also concluded that <i>Plectolyngbya</i> should be syn with <i>Leptolyngbya</i>, however we retain <i>Plectolyngbya</i> as our analyses do not support their conclusion. 


# Changes in 1.2

The manuscript was accept in the Journal of Phycology (yay), but open access is expensive (not yay). Contact me for a copy!

Taxonomy file was update to include the type/reference strain for each genus. 
In the case of provisional genera or a genus with an unsequenced type, the provided information is what we are considering the reference strain.

<i>Anathece</i> sequences were transfered to <i>Cyanobium</i> due to high similarity between the types strains. <i>Cyanobium</i> takes principle of priority. 
This also increases the number of assigned ASVs for the freshwater Prochlorococcaceae in our tests

Phylum name was changed to Cyanobacteriota

Chrysosporum ovalisporum was changed to <i>Umezakia ovalisporum, C. bergii</i> remain unchanged

Chroococcidiopsidaceae_X was changed to <i>Soropora</i>

Leptolyngbyaceae_X was changed to <i>Monilinema</i>

Leptolyngbyaceae_XX was changed to <i>Apatinema</i>

Leptolyngbyaceae_XXXXX was changed to <i>Romeriopsis</i>

<i>Geitlerinema</i> was incorrectly placed in the Oscillatoriales, this has been resolved

## New taxa added:
<i>Spirirestis, Microlinema, Cyanodorina, Microcrocis, Khargia</i>


## Problematic taxa:
<i>Insularia</i> was found to be monophyletic with, and >97.8% similar to, <i>Thainema</i>. Thus, this was included in the database as <i>Thainema</i>.

<i>Microcoleusiopsis</i> was found to be monophyletic, and >99.8% similar to <i>Arthrospira</i>. Thus, this was included in the database as <i>Arthrospira</i>.


# Changes in 1.1.2

<i>Calothrix</i> reevaluated based on Kumar et al. 2022, seems like no known Calothrix sequences exist.
<i>Fulbrightiella</i> and Sherwoodiella</i> replace Calothrix</i>
<i>Romeriopsis</i> replaces Leptolyngbyaceae_XXXXX
<i>Apatinema</i> replaces Leptolyngbyaceae_XX

## New taxa added:
<i>Radiocystis*, Edaphophycus, Zarconia, Neocylindrospermum, Prochlorotrichaceae_X</i>


<i>Radiocystis</i>: No valid publications, but these three strains form a monophyletic clade. However, these sequences come from trusted cyanobacterial taxonomists. Future works should aim to validate this genus. 
