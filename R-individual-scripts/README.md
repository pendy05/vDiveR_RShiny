[To the main README >](https://github.com/pendy05/DiveR/blob/main/README.md)

# R code for viral diversity dynamics plots and table
These six R scripts are to plot the viral diversity and showcase the patterns of the dynamics of the protein sequence change and the frequency distribution of each of the diversity motifs of viruses as well as completely and/or highly conserved <i>k</i>-mer sequence.<br><br>
## Table of contents
- [Plot 1: Entropy and Incidence of Total Variants](#plot-1-entropy-and-incidence-of-total-variants)
- [Plot 2: Correlation of Entropy](#plot-2-correlation-of-entropy)
- [Plot 3: Dynamics of Diversity Motifs (Proteome)](#plot-3-dynamics-of-diversity-motifs-proteome)
- [Plot 4: Dynamics of Diversity Motifs (Protein)](#plot-4-dynamics-of-diversity-motifs-protein)
- [Plot 5: Frequency Distribution Violin Plots (Protein)](#plot-5-frequency-distribution-violin-plots-protein)
- [Table 1: Identification of Completely (CCS) / Highly Conserved (HCS) Sequences](#table-1)
<br>

## Plot 1: Entropy and Incidence of Total Variants
**R script:** [Entropy-and-Incidence-of-Total-Variants.R](https://github.com/pendy05/DiveR/blob/main/R-individual-scripts/Entropy-and-Incidence-of-Total-Variants.R) 
<p align="center">
<img src="https://user-images.githubusercontent.com/74284859/132232581-ed1a8d49-aefb-4a7c-89bd-09053065e198.jpg" width="1000" height="500">
</p>

**Description:** <br>
Entropy (black) and incidence of total variants (pink) were measured for each aligned *k*-mer position (1-k, 2-k+1, etc.) of the proteins. The entropy values indicate the level of variability at the corresponding nonamer positions, with zero representing completely conserved positions (total variants incidence of 0%). Benchmark reference for values for entropy (black dotted line; 9.2) and total variants (pink dotted line; 98%) that from HIV-1 clade B envelope protein [(Hu *et al.*, 2013)](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0059994) are provided. For both individual protein and across proteome, the minimum entropy value is zero while the maximum entropy value at y-axis is 100.

<!--
Interpertation of result
-->

<br><br>

## Plot 2: Correlation of Entropy
**R script:** [Correlation-of-Entropy.R](https://github.com/pendy05/DiveR/blob/main/R-individual-scripts/Correlation-of-Entropy.R)
<p align="center">
<img src="https://user-images.githubusercontent.com/74284859/132232809-2bc5a305-5570-4528-ba41-9819a6f67a09.jpg" width="500" height="500">
</p>

**Description:** <br>
Relationship between incidence of total variants and entropy for viral proteome nonamer positions.  At y-axis, the minimum entropy value is zero while the maximum entropy value is obtained by rounding the highest entropy encountered up to integer.

<!--
Interpertation of result
-->

<br><br>

## Plot 3: Dynamics of Diversity Motifs (Proteome)
**R script:** [Dynamics-of-Diversity-Motifs(Proteome).R](https://github.com/pendy05/DiveR/blob/main/R-individual-scripts/Dynamics-of-Diversity-Motifs(Proteome).R)
<p align="center">
<img src="https://user-images.githubusercontent.com/74284859/132233240-e63f6f3f-c0cf-4e73-9879-adcf4dac4d44.jpg" width="500" height="600">
</p>

**Description:** <br>
*k*-mers are classified into four different motifs, namely index, major, minor and unique, based on their incidences. *k*-merTypes defines as distinct sequence for a given *k*-mer position. The above dot plot showcases the relationshop between the distribution of four distinct motifs and mutations. The diversity of the position is depicted by the decline of the index incidences (black), the increase of total variant incidences (pink) and corresponding individual patterns of the major, minor, unique and *k*-merTypes motifs. The below violin plot demonstrates the frequency distribution of the motifs. The width of the plot (x-axis) represents the frequency distribution of a given incidence of the indicated motif. The black thick horizontal line of box plot in the middle represents the median incidence value.

<!--
Interpertation of result
-->

<br><br>

## Plot 4: Dynamics of Diversity Motifs (Protein)
**R script:** [Dynamics-of-Diversity-Motifs(Protein).R](https://github.com/pendy05/DiveR/blob/main/R-individual-scripts/Dynamics-of-Diversity-Motifs(Protein).R)
<p align="center">
<img src="https://user-images.githubusercontent.com/74284859/132233808-fa6dd8b4-12c8-4823-ad3f-01271c30550f.jpg" width="500" height="600">
</p>

**Description:** <br>
*k*-mers are classified into four different motifs, namely index, major, minor and unique, based on their incidences. *k*-merTypes defines as distinct sequence for a given *k*-mer position. The above dot plot showcases the relationshop between the distribution of four distinct motifs and mutations. The diversity of the position is depicted by the decline of the index incidences (black), the increase of total variant incidences (pink) and corresponding individual patterns of the major, minor, unique and *k*-merTypes motifs. The below violin plot demonstrates the frequency distribution of the motifs. The width of the plot (x-axis) represents the frequency distribution of a given incidence of the indicated motif. The black thick horizontal line of box plot in the middle represents the median incidence value.

<!--
Interpertation of result
-->

<br><br>

## Plot 5: Frequency Distribution Violin Plots (Protein)
**R script:** [Frequency-Distribution-Violin-Plots(Protein).R](https://github.com/pendy05/DiveR/blob/main/R-individual-scripts/Frequency-Distribution-Violin-Plots(Protein).R)
<p align="center">
<img src="https://user-images.githubusercontent.com/74284859/132236375-0639ab97-1452-46ed-b18e-26540a96ead8.jpg" width="1000" height="600">
</p>

**Description:** <br>
The *k*-mer positions of the proteome and the individual proteins were defined as completely conserved (black) ( 𝑖𝑛𝑑𝑒𝑥 𝑖𝑛𝑐𝑖𝑑𝑒𝑛𝑐𝑒 = 100% ), highly conserved (blue) (90% ≤ 𝑖𝑛𝑑𝑒𝑥 𝑖𝑛𝑐𝑖𝑑𝑒𝑛𝑐𝑒 < 100%), mixed variable (green) ( 20% < 𝑖𝑛𝑑𝑒𝑥 𝑖𝑛𝑐𝑖𝑑𝑒𝑛𝑐𝑒 ≤ 90%), highly diverse (purple) (10% < 𝑖𝑛𝑑𝑒𝑥 𝑖𝑛𝑐𝑖𝑑𝑒𝑛𝑐𝑒 ≤ 20%) and extremely diverse (pink) ( 𝑖𝑛𝑑𝑒𝑥 𝑖𝑛𝑐𝑖𝑑𝑒𝑛𝑐𝑒 ≤ 10% ). 

<!--
Interpertation of result
-->

<br><br>
## Table 1: Identification of Completely (CCS) / Highly Conserved (HCS) Sequences <a id="table-1"></a>
**R script:** [Conserved-seq-concatenation.R](https://github.com/pendy05/DiveR/blob/main/R-individual-scripts/Conserved-seq-concatenation.R) 
<p align="center">
<img src="https://user-images.githubusercontent.com/74284859/159927864-0b1a7c1b-16d2-4aa6-b24a-388bb23dba5a.png" width="1000" height="500">
</p>

**Description:** <br>
The sample image shown above is HCS table produced for two proteins, namely NS1 and VP1. The table has three columns: <br>
1) CCS/HCS: consists of the header in format of "\[conservationlevel\]\_[proteinname]\_[number]" with the numbering starts from one, <br>
2) Position: indicates the start and end position of the concatenated (overlapping / neighbouring) CCS/HCS sequence and <br>
3) Sequence: CCS/HCS sequence. <br><br>
